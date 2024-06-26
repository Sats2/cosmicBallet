import Constants as const
import numpy as np
from typing import Union
import warnings
import math
from CelestialObjects import *


class _Fragments():
    """Internal class that holds the additional fragments generated from the collision of two planets

    Attributes:
        mass (float): Mass of the Fragment
        velocity (array): Velocity of the fragment
        force (array): Force acting on the fragment
        momentum (array): Momentum of the fragment
        radius (float): Radius of the fragment
        position (array): Position of the fragment in space.
        name (str): Generic Name for the Fragment
        material_property (dict): Material Property of the Fragment.
        volume (float): Volume of the Fragment
        density (float): Density of the Fragment
        planet_type (str): Define the object as a Fragment
    """
    def __init__(self, name:str, mass:float, velocity:np.array, radius:float, position:np.array, 
                 material_property:dict, force:np.array=None) -> None:
        """Constructor for the Fragments class

        Args:
            name (str): Name of the fragment (generic name)
            mass (float): Mass of the fragment
            velocity (np.array): Velocity of the fragment
            force (np.array): Force acting on the fragment
            radius (float): Radius of the fragment
            position (np.array): Position of the fragment in space.
            material_property (dict): Material Property of the fragment.
        """
        self.mass = mass
        self.velocity = velocity
        if force is not None:
            self.force = force
        else:
            self.force = np.zeros(3)
        self.momentum = mass * velocity
        self.radius = radius
        self.position = position
        self.name = name
        self.material_property = material_property
        self.volume = np.power(self.radius, 3) * np.pi * (4/3)
        self.density = self.mass / self.volume
        self.planet_tyoe = "fragment"


class Simulator():
    """This is a class that is used to simulate the orbital trajectory of celestial bodies.

    This simulator class performs the N-Body Simulation using Newtonian Mechanics for Orbital Motion using Newton's Second Law
    of Motion where the force is the cummulative of all gravitational interaction between objects. There are two simulation 
    methods available to solve, either with Lagranian Mechanics or with Hamiltonian Mechanics. For Langrangian Mechanics, the 
    available solvers are the forward euler and runge-kutta. As for Hamiltonian Mechanics, the simplectic euler and leapfrog
    solvers are available. High Density objects such as Black Holes and Neutron Stars cannot be solved with this simulator as
    Newtonian Mechanics fails to capture the orbital behaviour, instead use the simulator in MassiveObjectSimulator. Galaxies
    are also not simulated with this simulator, instead use GalaxySimulator.

    Attributes:
        celestial_bodies (list): A list of all celestial bodies used for N-Body Simulation, each item in list being an object
                                available in the CelestialObjects module.
        time_step (float/int): Discretization interval for the time used for simulation
        simulation_time (flaot/int): Total time for which simulation is needed.
        time_unit (str, optional): Unit of both the time values in the class. Defaults to seconds
        positions (array): An array containing the trajectory of all celestial objects.
    
    Methods:
        time_unit_correction(): Performs the conversion of values of the attributes 'simulation_time' and 'time_step' from
                                the entered unit to seconds that is later used for simulation.
        solve(simulation_method, solver): Performs the N-Body Simulation based on the input parameters and attributes of 
                                            the class and the user input choice of formulation and solver.

    """
    def __init__(self, celestial_bodies:list, time_step:Union[float,int], 
                 simulation_time:Union[float,int], time_unit:str="seconds") -> None:
        """Initializes the Simulator class.

        Args:
            celestial_bodies (list): A list of objects from the CelestialObjects module for which the simulation is needed.
            time_step (float,int): Time discretization value for the simulation
            simulation_time (float,int): Total time for simulation (same as end time)
            time_unit (str, optional): Unit of the arguements 'time_step' and 'simulation_time'. Defaults to seconds.

        Raises:
            TypeError: Raised when the input arguements do not match the intended data type.
            ValueError: Raised when parameter values are out-of-bounds.
        """
        try:
            assert isinstance(time_step, (float,int)), "Simulator attribute 'time_step' can only be of type float"
            assert isinstance(simulation_time, (float,int)), "Simulator attribute 'simulation time' can only be of type float/int"
            if time_unit is not None:
                assert isinstance(time_unit, str), "Simulator attribute 'time_unit' can only be of type string"
            for item in celestial_bodies:
                assert isinstance(item, (Stars,Planets,Galaxy,BlackHole)), "All items in the Simulator attribute 'celestial_bodies' must be class objects available in module CelestialObjects"
        except AssertionError:
            raise TypeError
        try:
            assert time_step>0, "Simulator attribute 'time_step' must be a positive value"
            assert simulation_time>0, "Simulator attribute 'simulation_time' must be a positive value"
            assert time_step<=simulation_time, "The chosen 'time_step' must be lesser than or equal to the 'simulation_time'"
            for item in celestial_bodies:
                assert (item.object_type != 'galaxy'), "The NBodySimulator cannot simulate Galaxies. Please use the Simulator in GalaxySimulator"
                assert (item.object_type != 'black_hole'), "The NBodySimulator cannot simulate Black Holes. Please use the Simulator in MassiveObjectSimulator"
                if item.object_type == "star":
                    assert (item.star_type != "Neutron"), "The NBodySimulator cannot simulate Neutron Stars. Please use the Simulator in MassiveObjectSimulator"
        except AssertionError:
            raise ValueError
        self.celestial_bodies = celestial_bodies
        self.time_step = time_step
        self.simulation_time = simulation_time
        self.time_unit = time_unit
        self.positions = None
        self.n_body = len(celestial_bodies)
        self.removed_object_list = []

    @property
    def time_unit_correction(self):
        """Updates the simulation_time and time_step attributes from the given unit to seconds

        Raises:
            ValueError: Raised when the user given unit does not match the available units for conversion
        """
        if self.time_unit is not None:
            if self.time_unit == "days":
                self.simulation_time *= const.DAY_TO_SEC
                self.time_step *= const.DAY_TO_SEC
            elif self.time_unit == "hours":
                self.simulation_time *= const.HOUR_TO_SEC
                self.time_step *= const.HOUR_TO_SEC
            elif self.time_unit == "months":
                self.simulation_time *= const.MON_TO_SEC
                self.time_step = const.MON_TO_SEC
            elif self.time_unit == "years":
                self.simulation_time *= const.YEAR_TO_sEC
                self.time_unit *= const.YEAR_TO_sEC
            else:
                raise ValueError("Time Unit Unrecognized. Select from 'hours/days/months/years'")

    def __merge_objects(self, p1:object, p2:object)->None:
        p1.mass += p2.mass
        new_volume = p1.volume + p2.volume
        new_radius = np.cbrt(new_volume * 0.75 / np.pi)
        p1.radius = new_radius
        self.removed_object_list.append(self.celestial_bodies.index(p2))
        self.celestial_bodies.remove(p2)
        p1.momentum = p1.mass * p1.velocity

    def __compute_fragments(self, planets:list, impact_energy:float)->list:
        p1, p2 = planets
        mass_list = []
        fragments_list = []
        for p in planets:
            obj_mass = p.mass
            obj_radius = p.radius
            if p is p1:
                other = p2
            else:
                other = p1
            num_fragments = np.random.rand(2,5)
            mass_fraction = 0.001 * np.random.rand(1,5)
            for _ in range(num_fragments):
                frag_velocity = np.sqrt(2 * impact_energy / other.mass) * np.sin(np.random(0,0.5*np.pi)) \
                    * other.velocity / np.linalg.norm(other.velocity)
                frag_position = frag_velocity * 10
                frag_mass = mass_fraction * p.mass
                frag_radius = np.cbrt(0.75 * p.mass / (p.density*np.pi))
                frag = _Fragments(name="Fragment", mass=frag_mass, velocity=frag_velocity, radius=frag_radius,
                                  position=frag_position, material_property=p.material_property)
                fragments_list.append(frag)
                obj_mass -= frag_mass
                obj_radius -= frag_radius
            mass_list.append([obj_mass, obj_radius])
        for i,p in enumerate(planets):
            p.mass = mass_list[i][0]
            p.radius = mass_list[i][1]
        return fragments_list
    
    def __elastic_collision(self, p1:object, p2:object)->None:
        velocity1 = p1.velocity
        velocity2 = p2.velocity
        p1.velocity = velocity2
        p2.velocity = velocity1
        p1.momentum = p1.mass * p1.velocity
        p2.momentum = p2.mass * p2.velocity
        return

    def __handle_collisions(self, p1:object, p2:object)->None:
        if p1.planet_type == "fragment" and p2.planet_type == "fragment":
            self.__elastic_collision(p1=p1, p2=p2)
            return
        impact_velocity = np.linalg.norm(p1.velocity - p2.velocity)
        if p1.mass > p2.mass:
            impactor = p2
            impacted = p2
        else:
            impactor = p1
            impacted = p2
        impact_energy = 0.5 * impactor.mass * impact_velocity**2
        if p1.planet_type.lower() == "rocky" and p2.planet_type.lower() == "rocky":
            if impact_energy > impacted.material_property["yield_strength"] * impacted.volume and \
                impacted.volume > 1e4*impactor.volume:
                fragments = self.__compute_fragments(impactor, impacted)
                self.celestial_bodies.append([item for item in fragments])
            else:
                self.__merge_objects(p1=impacted, p2=impactor)
        else:
            self.__merge_objects(p1=impacted, p2=impactor)
        return
    
    def __calculate_forces(self):
        """A private method within the Simulator Class that calculates the total acceleration acting on a celestial body in an N-Body Simulation. Collision checks are also
        performed to prevent numerical singularities in the force calculations.
        """
        for body in self.celestial_bodies:
            body.force[:] = 0.0
        for i, body1 in enumerate(self.celestial_bodies):
            for j, body2 in enumerate(self.celestial_bodies):
                if i != j:
                    r = body2.position - body1.position
                    distance = np.linalg.norm(r)
                    if distance < (body1.radius + body2.radius):
                        print(f"Collision detected between Celestial Objects {body1.name} and {body2.name}")
                        self.__handle_collisions(body1, body2)
                    elif distance > 0:
                        force_magnitude = const.G * body1.mass * body2.mass / np.power(distance, 2)
                        force_direction = r  / distance
                        body1.force += force_magnitude * force_direction
                        
    
    def __forward_euler_update(self):
        """Private method within the Simulator Class that updates the trajectory of the celestial objects calculated with the forward euler solver.
        """
        for body in self.celestial_bodies:
            body.velocity += (body.force * float(self.time_step / body.mass)) 
            body.position += body.velocity * self.time_step
    
    def __forward_euler(self):
        """Private method within the Simulator Class that implements the ODE Solver with Forward Euler method.
        """
        num_steps = math.ceil(self.simulation_time / self.time_step)
        for step in range(num_steps):
            if step == 0:
                for i,body in enumerate(self.celestial_bodies):
                    body.position = body.init_position.astype(np.float64)
                    body.velocity = body.init_velocity.astype(np.float64)
                    body.trajectory.append(body.position.copy())
            self.__calculate_forces()
            self.__forward_euler_update()
            for i, body in enumerate(self.celestial_bodies):
                body.trajectory.append(body.position.copy())
    
    def __rk4_step(self):
        """Private Method within the Simulator Class that performs the time update according to the Runge-Kutta method.
        """
        original_position = np.array([body.position for body in self.celestial_bodies])
        original_velocity = np.array([body.velocity for body in self.celestial_bodies])
        # Compute k1
        self.__calculate_forces()
        k1_r = self.time_step * original_velocity
        k1_v = self.time_step * np.array([body.force / body.mass for body in self.celestial_bodies])
        # Compute k2
        for i,body in enumerate(self.celestial_bodies):
            body.position = original_position[i] + 0.5*k1_r[i]
            body.velocity = original_velocity[i] + 0.5*k1_v[i]
        self.__calculate_forces()
        k2_r = self.time_step * np.array([body.velocity for body in self.celestial_bodies])
        k2_v = self.time_step * np.array([body.force / body.mass for body in self.celestial_bodies])
        # Compute k3
        for i,body in enumerate(self.celestial_bodies):
            body.position = original_position[i] + 0.5*k2_r[i]
            body.velocity = original_velocity[i] + 0.5*k2_v[i]
        self.__calculate_forces()
        k3_r = self.time_step * np.array([body.velocity for body in self.celestial_bodies])
        k3_v = self.time_step * np.array([body.force / body.mass for body in self.celestial_bodies])
        # Compute k4
        for i,body in enumerate(self.celestial_bodies):
            body.position = original_position[i] + k3_r[i]
            body.velocity = original_velocity[i] + k3_v[i]
        self.__calculate_forces()
        k4_r = self.time_step * np.array([body.velocity for body in self.celestial_bodies])
        k4_v = self.time_step * np.array([body.force / body.mass for body in self.celestial_bodies])
        # Position and Velocity Update
        for i,body in enumerate(self.celestial_bodies):
            body.position = original_position[i] + (k1_r[i] + 2*k2_r[i] + 2*k3_r[i] + k4_r[i])/6
            body.velocity = original_velocity[i] + (k1_v[i] + 2*k2_v[i] + 2*k3_v[i] + k4_v[i])/6

    def __runge_kutta4(self):
        """Private method of the Simulator Class that implements the 4th Order Runge-Kutta solver to solve the ODE System.
        """
        num_steps = math.ceil(self.simulation_time / self.time_step)
        for step in range(num_steps):
            if step == 0:
                for i,body in enumerate(self.celestial_bodies):
                    body.position = body.init_position.astype(np.float64)
                    body.velocity = body.init_velocity.astype(np.float64)
                    body.trajectory.append(body.position.copy())
            self.__rk4_step()
            for i,body in enumerate(self.celestial_bodies):
                body.trajectory.append(body.position.copy())
    
    def __calculate_potential_gradient(self):
        """Private function of the Simulator Class that generates the gradient of the potential energy for each object
        in a time slice that is used to solve the ODE System generated with Hamiltonian Dynamics.
        """
        num_body = len(self.celestial_bodies)
        self.potential_gradient = np.zeros((num_body,3))
        for i in range(num_body):
            for j in range(i+1,num_body):
                r = self.celestial_bodies[j].position - self.celestial_bodies[i].position
                distance = np.linalg.norm(r)
                if distance <= (self.celestial_bodies[j].radius + self.celestial_bodies[i].radius):
                    print(f"Collision detected between Celestial Objects {self.celestial_bodies[i].name} and \
                          {self.celestial_bodies[j].name}")
                    fragments = self.__handle_collisions(self.celestial_bodies[i], self.celestial_bodies[j])
                    if len(fragments) != 0:
                        for item in fragments:
                            self.celestial_bodies.append(item)
                else:
                    force_magnitude = const.G * self.celestial_bodies[i].mass * self.celestial_bodies[j].mass / distance**2
                    force = force_magnitude * r / distance
                    self.potential_gradient[i] += force
                    self.potential_gradient[j] -= force

    def __leapfrog_step(self):
        """Private method of the Simulator class that updates the momentum and position of all bodies being simulated 
        in a time slice for the Leapfrog Integrator.
        """
        self.__calculate_potential_gradient()
        for i,body in enumerate(self.celestial_bodies):
            body.momentum += 0.5 * self.time_step * self.potential_gradient[i]
        for body in self.celestial_bodies:
            body.position += self.time_step * body.momentum / body.mass
        self.__calculate_potential_gradient()
        for i,body in enumerate(self.celestial_bodies):
            body.momentum += 0.5 * self.time_step * self.potential_gradient[i]


    def __leapfrog(self):
        """Private method of the Simulator class that implements the Leapfrog Integrator to solve the ODE System that
        is generated from Hamiltonian Dynamics for the N-Body Problem using the Leapfrog Integrator scheme.
        """
        num_steps = math.ceil(self.simulation_time / self.time_step)
        for step in range(num_steps):
            if step == 0:
                for i,body in enumerate(self.celestial_bodies):
                    body.position = body.init_position.astype(np.float64)
                    body.momentum = body.init_velocity.astype(np.float64) * body.mass
                    body.trajectory.append(body.position.copy())
            self.__leapfrog_step()
            for i,body in enumerate(self.celestial_bodies):
                body.trajectory.append(body.position.copy())

    def __forest_ruth_step(self):
        """Private method of the Simulator class that updates the momentum and position of all bodies being simulated
        in a time slice for the Forest-Ruth Integrator.
        """
        gamma = 1 / (2 - np.cbrt(2))
        w1 = gamma / 2
        w2 = (1 - gamma) / 2
        w3 = w2
        w4 = w1
        steps = [w1, w2, w3, w4]
        for w in steps:
            for body in self.celestial_bodies:
                body.position += (w * self.time_step * body.momentum /body.mass)
            self.__calculate_potential_gradient()
            for i,body in enumerate(self.celestial_bodies):
                body.momentum += (w * self.time_step * self.potential_gradient[i])

    def __forest_ruth(self):
        """Private method of the Simulator class that implements the Forest-Ruth Integrator to solve the ODE System that
        is generated from Hamiltonian Dynamics for the N-Body Problem using the Forest-Ruth Integrator scheme.
        """
        num_steps =math.ceil(self.simulation_time / self.time_step)
        for stps in range(num_steps):
            if stps == 0:
                for i,body in enumerate(self.celestial_bodies):
                    body.position = body.init_position.astype(np.float64)
                    body.momentum = body.init_velocity.astype(np.float64) * body.mass
                    body.trajectory.append(body.position.copy())
            self.__forest_ruth_step()
            for i,body in enumerate(self.celestial_bodies):
                body.trajectory.append(body.position.copy())

    def solve(self, simulation_method:str="Lagrangian", solver:str="RK4")->None:
        """Method of the Simulator class that solves the ODE System of the N-Body Problem to generate results for the 
        N-Body Problem Simulation.

        The method contains two solvers for Lagrangian and Hamiltonian Mechanics each. The Forward Euler and 4th Order 
        Runge-Kutta method for Lagrangian Mechanics, and, the Leapfrog and Forest-Ruth methods for the Hamiltonian Mechanics
        formulation.

        Args:
            simulation_method (str, optional): The formulation of the N-Body Problem the method must follow to perform the
                                                simulation. Defaults to "Lagrangian".
            solver (str, optional): The integration scheme the method needs to use to perform the simulation. Defaults to "RK4".

        Raises:
            ValueError: Raised when an unrecognized solver or simulation scheme is entered.
        """ 
        if simulation_method.lower() == "lagrangian":
            if solver.lower() == "euler":
                if math.ceil(self.simulation_time / self.time_step) > 1e6:
                    warnings.warn("Number of Time Steps Too Large. Error accumulation may impact accuracy of results. Consider rk4 or Hamiltonian Mechanics to solve!")
                self.__forward_euler()
            elif solver.lower() == "rk4":
                self.__runge_kutta4()
            else:
                raise ValueError("Unidentified Solver. Select either Euler or RK4 (Runge-Kutta 4th Order)")
        elif simulation_method.lower() == "hamiltonian":
            if solver.lower() == "leapfrog":
                self.__leapfrog()
            elif solver.lower() == "forest_ruth":
                self.__forest_ruth()
            else:
                raise ValueError("Unidentified Solver. Select either Leapfrog or Forest_Ruth")
        else:
            raise ValueError("Unidentified Simulation Method. Select either Lagrangian or Hamiltonian")