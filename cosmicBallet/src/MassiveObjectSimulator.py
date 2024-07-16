import numpy as np
from CelestialObjects import BlackHole, Stars
from typing import Union
import math
import platform
import utils.Constants as const


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
    factor = - 1 * M / np.power(r*1, 2)
    term1 = (1 + 3*eta) * np.dot(v,v) - 1.5*eta*r_dot**2 - 2*(2+eta)*1*M/r
    term2 = 2 * (2 - eta) * r_dot * v
    F_PN1 = p1.mass * factor * (term1 * n - term2)
    return F_PN1


def _time_step_condition(stars:list, dense_object:object)->float:
    """Private Function of the script that calculates the minimum time step for the numerical integration.

    Args:
        stars (list): The list of all stars being simulated.
        dense_object (object): The dense object around which the stars are orbiting.

    Returns:
        float: The minimum time step for the numerical integration.
    """
    time_steps = []
    for star in stars:
        r = np.linalg.norm(star.init_position)
        t = np.sqrt(np.power(r,3) / (dense_object.mass * const.G))
        time_steps.append(t)
    return min(time_steps)/100


def _set_origin(stars:list, dense_object:object)->np.array:
    """Private Function of the script that sets the origin of the simulation to the dense object.

    Args:
        stars (list): The list of all stars being simulated.
        dense_object (object): The dense object around which the stars are orbiting.
    """
    position = dense_object.init_position
    for star in stars:
        star.init_position -= position
    dense_object.init_position = np.array([0,0,0]).astype(np.float64)
    return position


def _revert_origin(stars:list, dense_object:object, position:np.array)->None:
    """Private Function of the script that reverts the origin of the simulation to the original position.

    Args:
        stars (list): The list of all stars being simulated.
        dense_object (object): The dense object around which the stars are orbiting.
        position (np.array): The original position of the dense object.
    """
    for star in stars:
        for i in range(len(star.trajectory)):
            star.trajectory[i][1:] += position
    dense_object.position = position


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
        self.dense_body_position = _set_origin(stars, dense_body)
        try:
            assert time_step>0, "time_step must be a positive value"
            assert simulation_time>0, "simulation_time must be a positive value"
            min_timestep = _time_step_condition(stars, dense_body)
            assert time_step<=min_timestep, f"time_step must be lesser than {min_timestep}"
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
            star.init_position *= 0.5*self.dense_body.radius
            star.init_velocity *= const.C
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
            star.force -= self.G * star.mass * self.dense_body.mass * star.position / r**3
            correction_force = _calculate_PN1(star, self.dense_body)
            star.force += correction_force
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
                            correction_force = _calculate_PN1(star1, star2)
                            star1.force += correction_force

    def __forest_ruth_step(self):
        """Private method of the SchwarzschildSimulator class that updates the velocity and position of all bodies being simulated
        in a time slice for the Forest-Ruth Integrator.
        """
        gamma = 1 / (2 - np.cbrt(2))
        w1 = gamma / 2
        w2 = (1 - gamma) / 2
        w3 = w2
        w4 = w1
        steps = [w1, w2, w3, w4]
        for w in steps:
            for body in self.stars:
                body.position += w * body.velocity * self.dt
            self.__calculate_force()
            for i,body in enumerate(self.stars):
                body.velocity += w *self.dt * body.force / body.mass

    def solve(self)->None:
        """Private method of the SchwarzschildSimulator class that performs the numerical integration for
        the ODE System using the Forest-Ruth method.
        """
        try:
            for star in self.stars:
                star.trajectory = []
        except:
            pass
        self.__convert_to_natural_units()
        for i in range(self.n):
            if i == 0:
                for star in self.stars:
                    star.position = star.init_position.astype(np.float64)
                    star.velocity = star.init_velocity.astype(np.float64)
                    star.trajectory.append(np.concatenate(([(i+1)*self.dt], star.position.copy())))
            self.__forest_ruth_step()
            for star in self.stars:
                star.trajectory.append(np.concatenate(([(i+1)*self.dt], star.position.copy())))
        self.__convert_to_SI_units()
        _revert_origin(self.stars, self.dense_body, self.dense_body_position)


class BinaryMerger():
    """Class that simulates the merger of two celestial objects.

    The class allows for the simulation of the merger of two celestial objects such as Neutron Stars or Black Holes. The
    simulation is performed using the NRSur7dq2 and SurfinBH libraries for the gravitational waveforms and the final
    remnant mass and spin calculations. The simulation is performed under the assumption that the two objects are in a
    binary system and are in the process of merging. The simulation is performed in natural units.

    Attributes:
        binary_system (list): List of two objects that are to be simulated for Binary Merger.
        q (float): The mass ratio of the two objects in the binary system.
        chi1 (np.array): The dimensionless spin of the first object in the binary system.
        chi2 (np.array): The dimensionless spin of the second object in the binary system.

    Methods:
        simulate(): Calculates the trajectory of the two objects as they merge.
    """
    def __init__(self, binary_system:list)->None:
        """Constructor for the BinaryMerger class.

        Args:
            binary_system (list): List of two objects that are to be simulated for Binary Merger.

        Raises:
            OSError: Raised when the Binary Merger simulation is attempted on a non-Linux system.
            ImportError: Raised when the NRSur7dq2 and SurfinBH libraries are not installed.
            TypeError: Raised when the binary_system is not a list, or the objects in the list are not of the same type,
                        or the objects are not Neutron Stars or BlackHole class.
            ValueError: Raised when the binary_system does not have exactly two objects.
        """
        try:
            assert platform.system()=="Linux", "Binary Merger simulation is only supported on Linux"
        except AssertionError:
            raise OSError
        try:
            import NRSur7dq2 as sur
            import SurfinBH as sBH
        except ImportError:
            raise ImportError("NRSur7dq2 and SurfinBH libraries are required for Binary Merger simulation")
        try:
            assert isinstance(binary_system, list), "binary_system must be a list"
            assert type(binary_system[0])==type(binary_system[1]), "Objects in the binary_system must be of the same type"
            for obj in binary_system:
                assert isinstance(obj, (Stars, BlackHole)), "Objects in the binary_system must be Neutron Stars or Black Holes"
                if isinstance(obj, Stars):
                    assert obj.star_type=="Neutron", f"{obj.name} is not a Neutron Star"
        except AssertionError:
            raise TypeError
        try:
            assert len(binary_system)==2, "Binary System must have only two objects"
        except AssertionError:
            raise ValueError
        self.binary_system = binary_system
        self.q = binary_system[0].mass / binary_system[1].mass
    
    @property
    def calculate_chi(self)->None:
        """Property of the BinaryMerger class that calculates the dimensionless spins of the two objects in the binary
        system.
        """
        self.chi1 = const.C * self.binary_system[0].angular_momentum / (const.G * self.binary_system[0].mass**2)
        self.chi2 = const.C * self.binary_system[1].angular_momentum / (const.G * self.binary_system[1].mass**2)
    
    def simulate(self):
        pass

    def animate(self):
        pass