import numpy as np
from CelestialObjects import BlackHole, Stars
from typing import Union
import math
import Constants as const


class SchwarzschildSimulator():
    """Class that simulates the trajectory of stars orbiting a dense object.

    The class allows for the orbital simulation of stars around a dense object such as a (stationary) Black Hole or a
    Neutron Star by utilizing the Schwarzschild solution to Einstein's theory of General Relativity. The simulation is
    performed using the Schwarzschild metric for space-time. The simulation is performed under the assumption that the
    mass of the dense body is negligible compared to the mass of the stars, hence the force acting on the dense body can
    be ignored.

    For more information on Schwarzschild metric and equations of motion, please refer to the article:
    - https://www.physics.usu.edu/Wheeler/GenRel/Lectures/GRNotesDecSchwarzschildGeodesicsPost.pdf

    Newtonian equations of motion are used to simulate the motion of the stars. The simulation is performed using the
    assumption that the stars are point masses and the dense object is a stationary Black Hole or Neutron Star, with all
    calculations performed in natural units.

    Attributes:
        dense_body (object): The Neutron Star or Black Hole object.
        stars (list): A list containing the objects for stars that are orbiting the dense body.
        dt (int/float): The time step value for numerical integration.
        t_end (int/float): The end time for numerical integration.
        n (int): The total number of time steps in the numerical integration.
        G (int): Newton's Gravitational constant in natural units (=1).
        bh_mass (float): Mass of the Dense Body in SI Units.

    Methods:
        solve(): Performs the numerical integration for orbital simulation.
    """
    def __init__(self, dense_body:object, stars:list, time_step:Union[float,int],
                 simulation_time:Union[float,int]) -> None:
        """Constructor for the SchwarzschildSimulator class.

        Args:
            dense_body (object): The Neutron Star or Black Hole as an initialize object.
            stars (list): A list containing the objects for stars that are orbiting the dense body.
            time_step (Union[float,int]): The time step for numerical integration.
            simulation_time (Union[float,int]): The end time for numerical integration.

        Raises:
            TypeError: Raised when the time_step or simulation_time are not numerical values, or if the objects in the
                        list of stars does not belong to the Stars class or is Neutron Star, or the dense_body is not
                        a Black Hole or Neutron Star.
            ValueError: Raised when the time_step or simulation_time are lesser than or equal to zero, or when the 
                        time_step is too large when compared to the simulation_time, or when the dense body is not stationary.
        """
        try:
            assert isinstance(time_step, (int,float)), "time_step can only be a float/int value"
            assert isinstance(simulation_time, (int,float)), "simulation_time can only be a float/int value"
            for star in stars:
                assert (star.star_type != "Neutron"), f"Neutron Star orbit cannot be simulated. {star.name} is a Neutron Star"
                assert isinstance(star, Stars), f"{star.name} must be a star"
            if isinstance(dense_body, Stars):
                assert (star.star_type == "Neutron"), "dense_body must be a Neutron Star or Black Hole"
            else:
                assert isinstance(dense_body, BlackHole), "dense_body must be a Neutron Star or Black Hole"
            assert dense_body.init_velocity.all() == 0, "dense_body must be stationary"
        except AssertionError:
            raise TypeError
        try:
            assert time_step>0, "time_step must be a positive value"
            assert simulation_time>0, "simulation_time must be a positive value"
            assert time_step<0.25*simulation_time, "time_step is too large compared to simulation_time"
        except AssertionError:
            raise ValueError
        self.dense_body = dense_body
        self.stars = stars
        self.dt = time_step
        self.t_end = simulation_time
        self.n = math.ceil(simulation_time / time_step)
    
    def __convert_to_natural_units(self)->None:
        """Private Method of the SchwarzschildSimulator class that converts the input values of the celestial objects
        to natural units from SI units.
        """
        self.G = 1
        self.bh_mass = self.dense_body.mass
        self.bh_radius = self.dense_body.radius
        self.dense_body.mass = 1
        for star in self.stars:
            star.mass = star.mass / self.bh_mass
            star.radius *= (2 / self.bh_radius)
            star.init_position *= (2 / self.bh_radius)
            star.init_velocity *= 1 / const.C
    
    def __convert_to_SI_units(self)->None:
        """Private Method of the SchwarzschildSimulator class that converts the input values of the celestial objects
        to SI units from natural units.
        """
        self.dense_body.mass = self.bh_mass
        for star in self.stars:
            star.mass *= self.bh_mass
            star.radius *= 0.5*self.dense_body.radius
            for i in range(len(star.trajectory)):
                star.trajectory[i] *= 0.5*self.dense_body.radius

    def __calculate_force(self)->None:
        """Private method of the SchwarzschildSimulator class that computes the total force acting on the stars as they
        orbit the central dense object.

        Raises:
            RuntimeWarning: Raised when a star is at the center of the Black Hole or when two stars collide.
        """
        for star in self.stars:
            star.force[:] = 0
            r = np.linalg.norm(star.position)
            if r == 0:
                print(f"{star.name} is at the center of the Black Hole")
                raise RuntimeWarning
            star.force = -self.G * star.mass * self.dense_body.mass * star.position / r**3
        if len(self.stars) > 1:
    	        for i,star1 in enumerate(self.stars): 
                    for j,star2 in enumerate(self.stars):
                        if i != j:
                            r = np.linalg.norm(star1.position - star2.position)
                            if r <= 0:
                                print(f"{star1.name} and {star2.name} are at the same position")
                                raise RuntimeWarning
                            n = (star1.position - star2.position) / r
                            star1.force -= self.G * star1.mass * star2.mass * n / r**2

    def __rk4_step(self)->None:
        """Private method of the SchwarzschildSimulator class that performs the numerical integration for the system of ODEs
        based on the Runge-Kutta method for one time step
        """
        original_position = np.array([star.position for star in self.stars])
        original_velocity = np.array([star.velocity for star in self.stars])
        self.__calculate_force()
        k1_v = self.dt * np.array([star.force/star.mass for star in self.stars])
        k1_r = self.dt * original_velocity
        for i,star in enumerate(self.stars):
            star.position = original_position[i] + 0.5 * k1_r[i]
            star.velocity = original_velocity[i] + 0.5 * k1_v[i]
        self.__calculate_force()
        k2_v = self.dt * np.array([star.force/star.mass for star in self.stars])
        k2_r = self.dt * np.array([star.velocity for star in self.stars])
        for i,star in enumerate(self.stars):
            star.position = original_position[i] + 0.5 * k2_r[i]
            star.velocity = original_velocity[i] + 0.5 * k2_v[i]
        self.__calculate_force()
        k3_v = self.dt * np.array([star.force/star.mass for star in self.stars])
        k3_r = self.dt * np.array([star.velocity for star in self.stars])
        for i,star in enumerate(self.stars):
            star.position = original_position[i] + k3_r[i]
            star.velocity = original_velocity[i] + k3_v[i]
        self.__calculate_force()
        k4_v = self.dt * np.array([star.force/star.mass for star in self.stars])
        k4_r = self.dt * np.array([star.velocity for star in self.stars])
        for i,star in enumerate(self.stars):
            star.position = original_position[i] + (k1_r[i] + 2*k2_r[i] + 2*k3_r[i] + k4_r[i])/6
            star.velocity = original_velocity[i] + (k1_v[i] + 2*k2_v[i] + 2*k3_v[i] + k4_v[i])/6

    def solve(self)->None:
        """Private method of the SchwarzschildSimulator class that performs the numerical integration for
        the ODE System using the Runge-Kutta method.
        """
        self.__convert_to_natural_units()
        for i in range(self.n):
            if i == 0:
                for star in self.stars:
                    star.position = star.init_position.astype(np.float64)
                    star.velocity = star.init_velocity.astype(np.float64)
                    star.trajectory.append(np.concatenate(([(i+1)*self.dt], star.position.copy())))
            self.__rk4_step()
            for star in self.stars:
                star.trajectory.append(np.concatenate(([(i+1)*self.dt], star.position.copy())))
        self.__convert_to_SI_units()
            