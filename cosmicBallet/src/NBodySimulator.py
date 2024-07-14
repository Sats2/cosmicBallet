import Constants as const
from typing import Union
import warnings
import math
from CelestialObjects import *


def _calculate_PN1(p1:object, p2:object)->np.array:
    """Private function of the script that calculates the first term of the Post-Newtonian Expansion

    The 1PN term accounts for relativistic correction accounting for speed of the object and the speed of light.

    Args:
        p1 (object): One of the bodies in motion.
        p2 (object): The other body in motion.

    Returns:
        F_PN1 (np.array): The first term of the Post-Newtonian Expansion.
    """
    mu = p1.mass * p2.mass
    M = p1.mass + p2.mass
    eta = mu / M**2
    r = np.linalg.norm(p2.position - p1.position)
    n = (p1.position - p2.position) / r
    v = p1.velocity - p2.velocity
    r_dot = np.dot(v,n)
    factor = - const.G * M / np.power(r*const.C, 2)
    term1 = (1 + 3*eta) * np.dot(v,v) - 1.5*eta*r_dot**2 - 2*(2+eta)*const.G*M/r
    term2 = 2 * (2 - eta) * r_dot * v
    F_PN1 = p1.mass * factor * (term1 * n + term2)
    return F_PN1


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
        removed_object_list (list): A list containing all the celestial object that merged into other objects.
        formulation (str): The type of ODE Formulation used -> Newtonian or Hamiltonian
        solver (str): The solver used to solve the system of ODEs
        post_newton_correction (bool): User input whether to use Post Newton Correction of 1st order for the calculated forces.
        total_collisions (int): Holds the total number of collisions detected during the simulation
        fragment_collisions (int): Holds the total number of collisions between fragments that are detected during the simulation.
        time_val (float): Holds the current time during the simulation. Begins at 0 and ends at simulation_time.
        fragmentation_occured (bool): Holds whether fragmentation has occured during the simulation or not at a tim slice.
    
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
                assert isinstance(item, (Stars,Planets,BlackHole)), "All items in the Simulator attribute 'celestial_bodies' must be class objects available in module CelestialObjects"
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
        self.total_collisions = 0
        self.fragment_collisions = 0
        self.removed_object_list = []
        self.time_val = 0.0
        self.fragmentation_occured = False

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
        """Private method of the Simulation class that merges two objects on collision

        The total volume and mass of the colliding objects are preserved. The mass of the impactor gets added
        to the mass of the impacted and the radius changes accordingly to preserve the total volume of both 
        colliding objects.

        Args:
            p1 (object): The impacted planet that will hold the colliding object.
            p2 (object): The impactor planet that gets merged into the colliding planet.
        """
        p1.mass += p2.mass
        new_volume = p1.volume + p2.volume
        new_radius = np.cbrt(new_volume * 0.75 / np.pi)
        p1.radius = new_radius
        self.removed_object_list.append(p2)
        self.celestial_bodies.remove(p2)
        p1.momentum = p1.mass * p1.velocity

    def __compute_fragments(self, planets:list, fragmentation:list)->list:
        """Private method of the Simulator class that computes the fragmentation of planets on collision and
        generates the fragments. The mass of the fragments are based on the ratio of the kinetic energy of the
        collision to the gravitational binding energy of the planet (max. value of 0.8).

        Args:
            planets (list): List of planets that collide with each other of length 2.
            fragmentation (list): List containing whether the planets have fragmented or not.

        Returns:
            fragment_list (list): The total fragments generated that gets added to the class's celestial_objects
                                attribute for further computation of the N-body problem.
        """
        p1, p2 = planets
        fragments_list = []
        v = p2.velocity - p1.velocity
        kinetic_energy = 0.5 * (p1.mass * p2.mass) * np.linalg.norm(v)**2 / (p1.mass + p2.mass)
        for i,p in enumerate([p1,p2]):
            if fragmentation[i]:
                planet_mass = p.mass
                gravitational_binding_energy = 3 * const.G * p.mass**2 / (5 * p.radius)
                if kinetic_energy > 2*gravitational_binding_energy:
                    mass_fraction = 0.5
                    num_fragments = 1
                    break_apart = True
                else:
                    mass_fraction = 0.05
                    num_fragments = 5
                    break_apart = False
                for j in range(num_fragments):
                    if break_apart:
                        mass_frac = np.random.normal(loc=mass_fraction, scale=0.005)
                    else:
                        mass_frac = mass_fraction
                    mass = mass_frac * p.mass
                    volume = mass / p.density
                    radius = np.cbrt(volume * 0.75 / np.pi)
                    frag_velocity = np.sqrt(2 * p.mass * np.dot(p.velocity,p.velocity) / mass) * p.position \
                                    * np.random.normal(loc=1, scale=0.5) / np.linalg.norm(p.position)
                    frag_position = p.position + self.time_step * frag_velocity
                    fragment = Fragments(name=f"{p.name}_fragment_{j}", mass=mass, radius=radius, 
                                         material_property=p.material_property, position=frag_position,
                                         velocity=frag_velocity)
                    fragments_list.append(fragment)
                    planet_mass -= mass
                    if planet_mass < 0.25 * p.mass:
                        break
                p.radius = np.cbrt(planet_mass * 0.75 / (np.pi * p.density))
                p.mass = planet_mass
                print(p.name, p.mass, p.radius)
            else:
                pass
        return fragments_list
    
    def __elastic_collision(self, p1:object, p2:object)->None:
        """Private method of the Simulation class that handles perfect elastic collisions between fragments

        Args:
            p1 (object): Object for one of the fragments
            p2 (object): Object for the other fragment
        """
        velocity1 = p1.velocity
        velocity2 = p2.velocity
        p1.velocity = velocity2
        p2.velocity = velocity1
        p1.momentum = p1.mass * p1.velocity
        p2.momentum = p2.mass * p2.velocity
        return
    
    def __impact_fragments(self, p1: object, p2: object) -> list:
        """Calculates the impact fragments resulting from the collision between two objects.

        Args:
            p1 (object): The first object involved in the collision.
            p2 (object): The second object involved in the collision.

        Returns:
            fragments (list): A list indicating whether each object has fragmented after the collision.
                            The first element corresponds to p1 and the second element corresponds to p2.
                            True indicates fragmentation, while False indicates no fragmentation.
        """
        v = p2.velocity - p1.velocity
        kinetic_energy = 0.5 * (p1.mass * p2.mass) * np.linalg.norm(v)**2 / (p1.mass + p2.mass)
        impact_force_p1 = kinetic_energy / p1.radius
        impact_force_p2 = kinetic_energy / p2.radius
        impact_area = np.pi * np.sqrt(p1.radius*p2.radius / (p1.radius + p2.radius))**2
        impact_stress_p1 = impact_force_p1 / impact_area
        impact_stress_p2 = impact_force_p2 / impact_area
        fragments = [False, False]
        if impact_stress_p1 > p1.material_property["yield_strength"]:
            fragments[0] = True
        if impact_stress_p2 > p2.material_property["yield_strength"]:
            fragments[1] = True
        return fragments

    def __handle_collisions(self, p1:object, p2:object)->None:
        """Private method of the Simulator class that handles any detected collisions

        The collision handling is performed with the following criteria:
            - Star and Star collision: The stars merge into a single star
            - Star and Planet/Fragment collision: The planet/fragments gets merged into the star
            - Planet and Planet collision: Merge of planets or fragmentation of planets depending on the impact energy
            - Planet and Fragment collision: Fragment gets merged into the planet
            - Fragment and Fragment collision: Perfect elastic collision between the fragments

        Args:
            p1 (object): One of the colliding objects
            p2 (object): The other colliding object
        """
        self.total_collisions += 1
        collision_line = f"Collision detected between Celestial Objects {p1.name} and {p2.name} at position {p1.position} \
              at time {self.time_val} s"
        if p1.object_type == "star" and p2.object_type == "star":
            if p1.mass > p2.mass:
                impactor = p2
                impacted = p1
            else:
                impactor = p1
                impacted = p2
            self.__merge_objects(p1=impacted, p2=impactor)
        elif p1.object_type == "star" or p2.object_type == "star":
            if p1.object_type == "star":
                impactor = p2
                impacted = p1
            else:
                impactor = p1
                impacted = p2
            self.__merge_objects(p1=impacted, p2=impactor)
        else:
            if p1.planet_type == "fragment" and p2.planet_type == "fragment":
                self.fragment_collisions += 1
                collision_line = ""
                self.__elastic_collision(p1=p1, p2=p2)
                return
            if p1.mass > p2.mass:
                impactor = p2
                impacted = p2
            else:
                impactor = p1
                impacted = p2
            fragments_happen = self.__impact_fragments(p1=impactor, p2=impacted)
            if p1.planet_type.lower() == "rocky" and p2.planet_type.lower() == "rocky":
                if any(fragments_happen):
                    fragments = self.__compute_fragments([impactor, impacted], fragmentation=fragments_happen)
                    for item in fragments:
                        self.celestial_bodies.append(item)
                    self.fragmentation_occured = True
                else:
                    self.__merge_objects(p1=impacted, p2=impactor)
            else:
                self.__merge_objects(p1=impacted, p2=impactor)
        if collision_line == "":
            pass
        else:
            print(collision_line)
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
                        self.__handle_collisions(body1, body2)
                    elif distance > 0:
                        if self.post_newton_correction:
                            n = r / np.linalg.norm(r)
                            force = const.G * body1.mass * body2.mass * n / np.power(distance,2)
                            force_PN1 = _calculate_PN1(body1, body2)
                            body1.force += (force + force_PN1)
                        else:
                            force_magnitude = const.G * body1.mass * body2.mass / np.power(distance,3)
                            body1.force += force_magnitude*r
                        
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
            self.time_val += self.time_step
            if step == 0:
                for i,body in enumerate(self.celestial_bodies):
                    body.position = body.init_position.astype(np.float64)
                    body.velocity = body.init_velocity.astype(np.float64)
                    body.trajectory.append(np.concatenate(([self.time_val], body.position.copy())))
            self.__calculate_forces()
            self.__forward_euler_update()
            for i, body in enumerate(self.celestial_bodies):
                body.trajectory.append(np.concatenate(([self.time_val], body.position.copy())))
    
    def __rk4_step(self):
        """Private Method within the Simulator Class that performs the time update according to the Runge-Kutta method.
        """
        original_position = np.array([body.position for body in self.celestial_bodies])
        original_velocity = np.array([body.velocity for body in self.celestial_bodies])
        # Compute k1
        self.__calculate_forces()
        if self.fragmentation_occured:
            return
        k1_r = self.time_step * original_velocity
        k1_v = self.time_step * np.array([body.force / body.mass for body in self.celestial_bodies])
        # Compute k2
        for i,body in enumerate(self.celestial_bodies):
            body.position = original_position[i] + 0.5*k1_r[i]
            body.velocity = original_velocity[i] + 0.5*k1_v[i]
        self.__calculate_forces()
        if self.fragmentation_occured:
            return
        k2_r = self.time_step * np.array([body.velocity for body in self.celestial_bodies])
        k2_v = self.time_step * np.array([body.force / body.mass for body in self.celestial_bodies])
        # Compute k3
        for i,body in enumerate(self.celestial_bodies):
            body.position = original_position[i] + 0.5*k2_r[i]
            body.velocity = original_velocity[i] + 0.5*k2_v[i]
        self.__calculate_forces()
        if self.fragmentation_occured:
            return
        k3_r = self.time_step * np.array([body.velocity for body in self.celestial_bodies])
        k3_v = self.time_step * np.array([body.force / body.mass for body in self.celestial_bodies])
        # Compute k4
        for i,body in enumerate(self.celestial_bodies):
            body.position = original_position[i] + k3_r[i]
            body.velocity = original_velocity[i] + k3_v[i]
        self.__calculate_forces()
        if self.fragmentation_occured:
            return
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
            self.time_val += self.time_step
            if step == 0:
                for i,body in enumerate(self.celestial_bodies):
                    body.position = body.init_position.astype(np.float64)
                    body.velocity = body.init_velocity.astype(np.float64)
                    body.trajectory.append(np.concatenate(([self.time_val], body.position.copy())))
            self.__rk4_step()
            for i,body in enumerate(self.celestial_bodies):
                body.trajectory.append(np.concatenate(([self.time_val], body.position.copy())))
    
    def __calculate_potential_gradient(self):
        """Private function of the Simulator Class that generates the gradient of the potential energy for each object
        in a time slice that is used to solve the ODE System generated with Hamiltonian Dynamics.
        """
        num_body = len(self.celestial_bodies)
        self.potential_gradient = np.zeros((num_body,3))
        for i in range(num_body):
            for j in range(i+1,num_body):
                try:
                    r = self.celestial_bodies[j].position - self.celestial_bodies[i].position
                except IndexError:
                    break
                distance = np.linalg.norm(r)
                if distance <= (self.celestial_bodies[j].radius + self.celestial_bodies[i].radius):
                    self.__handle_collisions(self.celestial_bodies[i], self.celestial_bodies[j])
                    if len(self.celestial_bodies) != num_body:
                        num_body = len(self.celestial_bodies)
                        self.__calculate_potential_gradient()
                else:
                    force_magnitude = const.G * self.celestial_bodies[i].mass * self.celestial_bodies[j].mass / distance**2
                    force = force_magnitude * r / distance
                    if self.post_newton_correction:
                        force_pn1 = _calculate_PN1(self.celestial_bodies[i], self.celestial_bodies[j])
                        force += force_pn1
                    self.potential_gradient[i] += force
                    self.potential_gradient[j] -= force

    def __leapfrog_step(self):
        """Private method of the Simulator class that updates the momentum and position of all bodies being simulated 
        in a time slice for the Leapfrog Integrator.
        """
        self.__calculate_potential_gradient()
        for i,body in enumerate(self.celestial_bodies):
            body.momentum += 0.5 * self.time_step * self.potential_gradient[i,:]
            body.velocity = body.momentum / body.mass
        for body in self.celestial_bodies:
            body.position += self.time_step * body.velocity
        self.__calculate_potential_gradient()
        if self.fragmentation_occured:
            return
        for i,body in enumerate(self.celestial_bodies):
            body.momentum += 0.5 * self.time_step * self.potential_gradient[i,:]
            body.velocity = body.momentum / body.mass


    def __leapfrog(self):
        """Private method of the Simulator class that implements the Leapfrog Integrator to solve the ODE System that
        is generated from Hamiltonian Dynamics for the N-Body Problem using the Leapfrog Integrator scheme.
        """
        num_steps = math.ceil(self.simulation_time / self.time_step)
        for step in range(num_steps):
            self.time_val += self.time_step
            if step == 0:
                for i,body in enumerate(self.celestial_bodies):
                    body.position = body.init_position.astype(np.float64)
                    body.momentum = body.init_velocity.astype(np.float64) * body.mass
                    body.velocity = body.momentum / body.mass
                    body.trajectory.append(np.concatenate(([self.time_val], body.position.copy())))
            self.__leapfrog_step()
            for i,body in enumerate(self.celestial_bodies):
                body.trajectory.append(np.concatenate(([self.time_val], body.position.copy())))

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
                body.position += w * body.velocity * self.time_step
            self.__calculate_potential_gradient()
            for i,body in enumerate(self.celestial_bodies):
                body.momentum += (w * self.time_step * self.potential_gradient[i])
                body.velocity = body.momentum / body.mass

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
                    body.velocity = body.momentum / body.mass
                    body.trajectory.append(np.concatenate(([self.time_val], body.position.copy())))
            self.__forest_ruth_step()
            for i,body in enumerate(self.celestial_bodies):
                body.trajectory.append(np.concatenate(([self.time_val], body.position.copy())))

    def solve(self, formulation:str="Hamiltonian", solver:str="forest_ruth", correction:bool=False)->None:
        """Method of the Simulator class that solves the ODE System of the N-Body Problem to generate results for the 
        N-Body Problem Simulation.

        The method contains two solvers for Newtonian and Hamiltonian Mechanics each. The Forward Euler and 4th Order 
        Runge-Kutta method for Newtonian Mechanics, and, the Leapfrog and Forest-Ruth methods for the Hamiltonian Mechanics
        formulation.

        Args:
            simulation_method (str, optional): The formulation of the N-Body Problem the method must follow to perform the
                                                simulation. Defaults to "Newtonian".
            solver (str, optional): The integration scheme the method needs to use to perform the simulation. Defaults to "RK4".
            correction (bool, optional): User input to whether the Post-Newton Correction term is applied or not. This is 
                                        applicable to cases where heavy objects are relatively close to each other. Applied 
                                        only to the Newtonian Formulation. Defaults to False.

        Raises:
            TypeError: Raised when the input arguements are not of the defined datatype.
            ValueError: Raised when an unrecognized solver or simulation scheme is entered.
        """ 
        try:
            for body in self.celestial_bodies:
                body.trajectory = []
            self.time_val = 0.0
        except:
            pass
        try:
            assert isinstance(correction, bool), "The arguement correction needs to be of type bool"
            assert isinstance(formulation, str), "The arguement formulation needs to be a string object"
            assert isinstance(solver, str), "The arguement solver needs to be a string object"
        except AssertionError:
            raise TypeError
        self.post_newton_correction = correction
        if formulation.lower() == "newtonian":
            if solver.lower() == "euler":
                if math.ceil(self.simulation_time / self.time_step) > 1e6:
                    warnings.warn("Number of Time Steps Too Large. Error accumulation may impact accuracy of results. Consider rk4 or Hamiltonian Mechanics to solve!")
                self.__forward_euler()
            elif solver.lower() == "rk4":
                self.__runge_kutta4()
            else:
                raise ValueError("Unidentified Solver. Select either Euler or RK4 (Runge-Kutta 4th Order)")
        elif formulation.lower() == "hamiltonian":
            if solver.lower() == "leapfrog":
                self.__leapfrog()
            elif solver.lower() == "forest_ruth":
                self.__forest_ruth()
            else:
                raise ValueError("Unidentified Solver. Select either Leapfrog or Forest_Ruth")
        else:
            raise ValueError("Unidentified Simulation Method. Select either Newtonian or Hamiltonian")
        if self.total_collisions == 0:
            return
        print(f"Total Collisions Detected: {self.total_collisions} out of which {self.fragment_collisions} are fragment collision.")