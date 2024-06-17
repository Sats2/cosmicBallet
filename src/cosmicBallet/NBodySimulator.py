import Constants as const
import numpy as np
from typing import Union
import warnings
from numba import njit
from CelestialObjects import *


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
        simulation_method (str): A string indicating user choise between Langrangian and Hamiltonian Mechanics for simulation
        sovler (str): The type of solver to be used for solving the system of Ordinary Differential Equations in the N-Body
                     Simulations.
        time_step (float/int): Discretization interval for the time used for simulation
        simulation_time (flaot/int): Total time for which simulation is needed.
        time_unit (str, optional): Unit of both the time values in the class. Defaults to seconds
    
    Methods:
        time_unit_correction(): Performs the conversion of values of the attributes 'simulation_time' and 'time_step' from
                                the entered unit to seconds that is later used for simulation.
        solve(): Performs the N-Body Simulation based on the input parameters and attributes of the class

    """
    def __init__(self, celestial_bodies:list, method:str, solver:str, time_step:Union[float,int], 
                 simulation_time:Union[float,int], time_unit:str="seconds") -> None:
        """Initializes the Simulator class.

        Args:
            celestial_bodies (list): A list of objects from the CelestialObjects module for which the simulation is needed.
            method (str): Selection between Langranian and Hamiltonian Mechanics
            solver (str): Choice of solver for simulation
            time_step (float,int): Time discretization value for the simulation
            simulation_time (float,int): Total time for simulation (same as end time)
            time_unit (str, optional): Unit of the arguements 'time_step' and 'simulation_time'. Defaults to seconds.

        Raises:
            TypeError: Raised when the input arguements do not match the intended data type.
            ValueError: Raised when parameter values are out-of-bounds.
        """
        try:
            assert isinstance(method, str), "Simulator attribute 'method' can only be of type string"
            assert isinstance(solver, str), "Simulator attribute 'solver' can only be of type string"
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
        self.simulation_method = method
        self.time_step = time_step
        self.simulation_time = simulation_time 
        self.solver = solver
        self.time_unit = time_unit

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
            
    def __collision_check(self):
        pass
    
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
                    if distance <= (body1.radius + body2.radius):
                        self.__collision_check()
                    else:
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
        self.positions = [[] for _ in self.celestial_bodies]
        num_steps = math.ceil(self.simulation_time / self.time_step)
        for step in range(num_steps):
            if step == 0:
                for body in self.celestial_bodies:
                    body.position = body.init_position.astype(np.float64)
                    body.velocity = body.init_velocity.astype(np.float64)
            self.__calculate_forces()
            self.__forward_euler_update()
            for i, body in enumerate(self.celestial_bodies):
                self.positions[i].append(body.position.copy())
    
    def __rk4_step(self):
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
        num_steps = math.ceil(self.simulation_time / self.time_step)
        self.positions = [[] for _ in self.celestial_bodies]
        for step in range(num_steps):
            if step == 0:
                for body in self.celestial_bodies:
                    body.position = body.init_position.astype(np.float64)
                    body.velocity = body.init_velocity.astype(np.float64)
            self.__rk4_step()
            for i,body in enumerate(self.celestial_bodies):
                self.positions[i].append(body.position.copy())
    
    def __leapfrog(self):
        pass

    def __simplectic_euler(self):
        pass

    #TODO: Complete the function calls to simulate the trajectories
    def solve(self)->None:
        if self.simulation_method.lower() == "lagrangian":
            if self.solver.lower() == "euler":
                if math.ceil(self.simulation_time / self.time_step) > 1e6:
                    warnings.warn("Number of Time Steps Too Large. Error accumulation may impact accuracy of results. Consider rk4 or Hamiltonian Mechanics to solve!")
                self.__forward_euler()
            elif self.solver.lower() == "rk4":
                self.__runge_kutta4()
            else:
                raise ValueError("Unidentified Solver. Select either Euler or RK4 (Runge-Kutta 4th Order)")
        elif self.simulation_method.lower() == "hamiltonian":
            if self.solver.lower() == "leapfrog":
                self.__leapfrog()
            elif self.solver.lower() == "simp_euler":
                self.__simplectic_euler()
            else:
                raise ValueError("Unidentified Solver. Select either Leapfrog or Simp_Euler (Simplectic Euler)")
        else:
            raise ValueError("Unidentified Simulation Method. Select either Lagrangian or Hamiltonian")