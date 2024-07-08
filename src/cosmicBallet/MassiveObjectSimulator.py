import numpy as np
from CelestialObjects import BlackHole, Stars
from typing import Union
import math
import Constants as const


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
    v = p2.velocity - p1.velocity
    r = np.linalg.norm(p2.position - p1.position)
    r_hat = (p2.position - p1.position)/r
    F_PN1 = const.G * mu * ((1 + 1.5*(v/const.C)**2 - 5*const.G*M/(r*const.C**2) + 0.5*(np.dot(v,r_hat)**2)/const.C**2)*r_hat \
                            - 4*np.dot(v,r_hat)*v/const.C**2) / r**2
    return F_PN1


def _calculate_PN2(p1:object, p2:object)->np.array:
    pass


def _calculate_PN3(p1:object, p2:object)->np.array:
    pass

class DenseObjectSimulator():
    """A class that simulates the trajectory of stars around a dense object in space.

    The class simulates the trajectory of stars around a dense object such as a black hole or neutron star using Newtonian Gravity
    with Parameterized Post-Newtonian Approximation for approximating general relativity. The time integration for the resulting
    Ordinary Differential Equation (ODE) is performed using 4th Order Runge-Kutta (RK4) Integration. 

    Attributes:
        celestial_bodies (list): A list containing all the celestial objects that are simulated. The first item in the list 
                                contains the dense object followed by the stars whose trajectory is simulated.
        time_step (int/float): The time step for RK4 Integration.
        simulation_time (int/float): The total amount of time for the trajectory simulation (also t_end in simulations).
        n (int): The total number of time steps in the simulation.
    
    Methods:
        solve(): A method that solves the N-Body problem using the RK4 Integration.
    """
    def __init__(self, dense_body:object, celestial_object_list:list, time_step:Union[int,float],
                 simulation_time:Union[int,float]) -> None:
        """Constructor for the class DenseObjectSimulator.

        Args:
            dense_body (object): The dense object being simulated.
            celestial_object_list (list): A list containing the total number of stars whose trajectory is simulated for.
            time_step (Union[int,float]): The time step for RK4 Integration
            simulation_time (Union[int,float]): The total time for trajectory simulation.

        Raises:
            TypeError: Raised when dense_body is not a Black Hole or Neutron Star, the time step is not a numerical data type,
                        the simulation_time is not a numerical data type, or the celestial_object_list is not a list type.
            ValueError: Raised when the time_step or simulation_time is not positive, the time_step is too large, the objects
                        in celestial_object_list are not Stars or are Neutron Stars.
        """ 
        try:
            assert isinstance(dense_body, (Stars, BlackHole)), "dense_body must be a Neutron Star or Black Hole"
            if dense_body.object_type == "Star":
                assert (dense_body.star_type == "Neutron"), "dense_body must be a Neutron Star or Black Hole"
            assert isinstance(time_step, (int,float)), "time_step must be an integer or floating point value"
            assert isinstance(simulation_time, (int,float)), "simulation_time must be an integer or floating point value"
            assert isinstance(celestial_object_list, list), "celestial_object_list must be a list of stars"
        except AssertionError:
            raise TypeError
        try:
            assert time_step>0, "time_step cannot be lesser than 0"
            assert simulation_time>0, "simulation_time cannot be lesser than 0"
            assert time_step<0.5*simulation_time, "time step cannot be larger than 50% of simulation_time"
            for body in celestial_object_list:
                assert isinstance(body, Stars), f"{body.name} is not a star. All objects in celestial_object_list needs to be a non-Neutron Star"
                assert (body.star_type!="Neutron"), f"{body.name} is a Neutron Star. All objects in celestial_object_list must be non-Neutron Stars"
        except AssertionError:
            raise ValueError
        self.celestial_objects = celestial_object_list
        self.celestial_objects.insert(0, dense_body)
        self.time_step = time_step
        self.simulation_time = simulation_time
        self.n = math.ceil(self.simulation_time / self.time_step)
    
    def __force_correction(self, b1:object, b2:object)->np.array:
        """Private Method of the DenseObjectSimulator class that computes the total force correction required

        Args:
            b1 (object): Object 1
            b2 (object): Object 2

        Returns:
            force_correction_term (np.array): An array containing the total correction value to the force based on PN expansion.
        """
        pn1_term = _calculate_PN1(b1, b2)
        pn2_term = _calculate_PN2(b1, b2)
        pn3_term = _calculate_PN3(b1, b2)
        force_correction_term = pn1_term + pn2_term + pn3_term
        return force_correction_term
    
    def __compute_force(self)->None:
        """Private method of the DenseObjectSimulator class that computes the total force acting on all objects.
        """
        for body in self.celestial_objects:
            body.force[:] = 0.0
        for i, body1 in enumerate(self.celestial_objects):
            for j,body2 in enumerate(self.celestial_objects):
                r = body2.position - body1.position
                distance = np.linalg.norm(r)
                if distance <= (body1.radius + body2.radius):
                    self.__merge_body(body1, body2)
                    print(f"Collision between {body1.name} and {body2.name} detected! Objects will merge.")
                else:
                    force = (const.G * body1.mass * body2.mass / np.power(distance,2)) * r/distance
                    force_correction_term = self.__force_correction(body1, body2)
                    body1.force += force + force_correction_term 

    def __rk4_step(self)->None:
        """A private method of the DenseObjectSimulator class that performs the RK4 integration in a single time step.
        """
        original_position = np.array([body.position for body in self.celestial_objects])
        original_velocity = np.array([body.velocity for body in self.celestial_objects])
        self.__compute_force()
        k1_r = self.time_step * original_velocity
        k1_v = self.time_step * np.array([body.force / body.mass for body in self.celestial_bodies])
        for i,body in enumerate(self.celestial_bodies):
            body.position = original_position[i] + 0.5*k1_r[i]
            body.velocity = original_velocity[i] + 0.5*k1_v[i]
        self.__compute_force()
        k2_r = self.time_step * np.array([body.velocity for body in self.celestial_bodies])
        k2_v = self.time_step * np.array([body.force / body.mass for body in self.celestial_bodies])
        for i,body in enumerate(self.celestial_bodies):
            body.position = original_position[i] + 0.5*k2_r[i]
            body.velocity = original_velocity[i] + 0.5*k2_v[i]
        self.__compute_force()
        k3_r = self.time_step * np.array([body.velocity for body in self.celestial_bodies])
        k3_v = self.time_step * np.array([body.force / body.mass for body in self.celestial_bodies])
        for i,body in enumerate(self.celestial_bodies):
            body.position = original_position[i] + k3_r[i]
            body.velocity = original_velocity[i] + k3_v[i]
        self.__compute_force()
        k4_r = self.time_step * np.array([body.velocity for body in self.celestial_bodies])
        k4_v = self.time_step * np.array([body.force / body.mass for body in self.celestial_bodies])
        for i,body in enumerate(self.celestial_bodies):
            body.position = original_position[i] + (k1_r[i] + 2*k2_r[i] + 2*k3_r[i] + k4_r[i])/6
            body.velocity = original_velocity[i] + (k1_v[i] + 2*k2_v[i] + 2*k3_v[i] + k4_v[i])/6

    def solve(self)->None:
        """Method of the DenseObjectSimulator class that performs the RK4 Integration and updates trajectory of objects.
        """
        for i in range(self.n):
            if i == 0:
                for body in enumerate(self.celestial_objects):
                    body.position = body.init_position
                    body.velocity = body.init_velocity
                    body.trajectory.append(body.position.copy())
            self.__rk4_step()
            for body in enumerate(self.celestial_objects):
                body.trajectory.append(body.position.copy())
            self.dense_body.trajectory.append(self.dense_body.position.copy())