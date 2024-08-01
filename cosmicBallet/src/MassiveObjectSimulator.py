import numpy as np
from CelestialObjects import BlackHole, Stars
from typing import Union
import math
import platform
from utils.Visualization import Visualize
import mayavi.mlab as mlab
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


def _animate_merger(times:np.array, traj:list, radius_list:list)->None:
    """Private Function of the script that animates the merger of two celestial objects.

    Args:
        times (np.array): The time values for the merger.
        traj (list): The trajectory of the two objects as they merge.
        radius_list (list): The radius of the two objects.
    """
    bht1, bht2 = traj
    trajectory1 = bht1.T
    trajectory2 = bht2.T
    rad1 = np.min(np.array(radius_list))
    for i,radius in enumerate(radius_list):
        radius_list[i] = radius / rad1

    # Create a figure for the animation
    mlab.figure(size=(1000, 1000), bgcolor=(0.35, 0.35, 0.35))

    # Function to create a black hole
    def create_black_hole(position, radius=1.0):
        sphere = mlab.points3d(position[0], position[1], position[2], scale_factor=radius, color=(0, 0, 0), resolution=100)
        return sphere

    # Function to create an accretion disk
    def create_accretion_disk(center, inner_radius=1.5, outer_radius=3.0, height=0.1, resolution=100):
        r = np.linspace(inner_radius, outer_radius, resolution)
        theta = np.linspace(0, 2*np.pi, resolution)
        r, theta = np.meshgrid(r, theta)
        x = (r * np.cos(theta)).ravel()
        y = (r * np.sin(theta)).ravel()
        z = (np.zeros_like(r) + center[2]).ravel()
        accretion_disk = mlab.mesh(x.reshape(resolution, resolution) + center[0],
                                y.reshape(resolution, resolution) + center[1],
                                z.reshape(resolution, resolution), color=(0.5, 0.5, 0.5), opacity=0.9)
        return accretion_disk

    # Initial positions
    bh1 = create_black_hole(trajectory1[0],radius_list[0])
    bh2 = create_black_hole(trajectory2[0],radius_list[1])
    disk1 = create_accretion_disk(trajectory1[0],inner_radius=radius_list[0]*1.5, outer_radius=radius_list[0]*2.0)
    disk2 = create_accretion_disk(trajectory2[0], inner_radius=radius_list[1]*1.5, outer_radius=radius_list[1]*2.0)

    # Add text to display the time
    time_text = mlab.text(0.7, 0.9, f'Time: {times[0]:.1f} M', color=(1, 1, 1), width=0.2)

    # Function to update black holes, accretion disks, and time text
    def update_scene(i):
        # Update positions based on trajectories
        bh1.mlab_source.set(x=trajectory1[i, 0], y=trajectory1[i, 1], z=trajectory1[i, 2])
        bh2.mlab_source.set(x=trajectory2[i, 0], y=trajectory2[i, 1], z=trajectory2[i, 2])
        
        disk1.mlab_source.set(x=trajectory1[i, 0] + disk1.mlab_source.x - disk1.mlab_source.x.mean(),
                            y=trajectory1[i, 1] + disk1.mlab_source.y - disk1.mlab_source.y.mean(),
                            z=trajectory1[i, 2] + disk1.mlab_source.z - disk1.mlab_source.z.mean())
        disk2.mlab_source.set(x=trajectory2[i, 0] + disk2.mlab_source.x - disk2.mlab_source.x.mean(),
                            y=trajectory2[i, 1] + disk2.mlab_source.y - disk2.mlab_source.y.mean(),
                            z=trajectory2[i, 2] + disk2.mlab_source.z - disk2.mlab_source.z.mean())
        
        # Update time text
        time_text.text = f'Time: {times[i]:.1f} M'

        # Merge black holes and disks at time=0
        if times[i] == 0:
            bh1.mlab_source.set(x=0, y=0, z=0, scale_factor=4)
            disk1.mlab_source.set(x=disk1.mlab_source.x - disk1.mlab_source.x.mean(), 
                                y=disk1.mlab_source.y - disk1.mlab_source.y.mean(), 
                                z=disk1.mlab_source.z - disk1.mlab_source.z.mean())
            bh2.remove()
            disk2.remove()
            time_text.text = 'Time: 0.0 (Merged)'

    # Set the camera view to an isometric perspective
    mlab.view(azimuth=45, elevation=45, distance=50)

    # Animate the black holes
    for i in range(len(times)):
        update_scene(i)
        mlab.process_ui_events()
        mlab.savefig(f'temp/frame_{i:03d}.png')

    mlab.close()



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
        dense_body_position (np.array): The initial position of the dense body.

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
            TypeError: When dense_body is not a Black Hole or Neutron Star 
            TypeError: When stars is not a list of stars
            TypeError: When object in stars list is not a Star object or is a Neutron Star
            TypeError: When time_step is not a float or int value
            TypeError: When simulation_time is not a float or int value
            ValueError: When time_step is not a positive value
            ValueError: When simulation_time is not a positive value
            ValueError: When time_step is greater than the minimum time step for the simulation
            ValueError: When all star objects in Stars list have common names.
        """
        try:
            assert isinstance(dense_body, (Stars, BlackHole))
        except AssertionError:
            raise TypeError("dense_body must be a Neutron Star or Black Hole")
        try:
            assert isinstance(stars, list)
        except AssertionError:
            raise TypeError("stars must be a list of stars")
        try:
            for star in stars:
                assert isinstance(star, Stars), "Object in stars must be a Star object"
                assert star.star_type != "Neutron", "Neutron Star orbit cannot be simulated"
        except AssertionError:
            raise TypeError
        try:
            assert isinstance(time_step, (int,float))
        except AssertionError:
            raise TypeError("time_step can only be a float/int value")
        try:
            assert isinstance(simulation_time, (int,float))
        except AssertionError:
            raise TypeError("simulation_time can only be a float/int value")
        try:
            assert time_step>0
        except AssertionError:
            raise ValueError("time_step must be a positive value")
        try:
            assert simulation_time>0
        except AssertionError:
            raise ValueError("simulation_time must be a positive value")
        min_timestep = _time_step_condition(stars, dense_body)
        try:
            assert time_step<=min_timestep
        except AssertionError:
            raise ValueError(f"time_step must be lesser than {min_timestep}")
        try:
            name_list = []
            for star in stars:
                name_list.append(star.name)
            assert len(name_list) == len(set(name_list))
        except AssertionError:
            raise ValueError("All star objects in Stars list must have unique names")
        self.dense_body = dense_body
        self.stars = stars
        self.dt = time_step
        self.t_end = simulation_time
        self.n = math.ceil(simulation_time / time_step)
        self.dense_body_position = _set_origin(stars, dense_body)
    
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

    def __compile_results(self)->list:
        """Private method of the SchwarzschildSimulator class that compiles the results of the simulation into a list.

        Returns:
            celestial_object_list (list): A list of all celestial objects in the simulation.
        """
        self.dense_body.trajectory = []
        for i in range(self.n+1):
            self.dense_body.trajectory.append(np.concatenate(([(i+1)*self.dt], self.dense_body.init_position.copy())))
        celestial_object_list = [self.dense_body]
        for star in self.stars:
            celestial_object_list.append(star)
        return celestial_object_list

    def visualize(self, visualization_type:str="scientific", save_figure:bool=False, figure_name:str=None,
                  animate:bool=False)->None:
        """Method of the SchwarzschildSimulator class that visualizes the trajectory of the stars.

        The method uses the Visualize class from the Visualization module to generate the visualization of the trajectory of the 
        stars around a central dense body.

        Args:
            visualization_type (str, optional): The type of visualization the method needs to generate. Defaults to "scientific".
            save_figure (bool, optional): User input to whether the generated visualization needs to be saved or not. Defaults to False.
            figure_name (str, optional): The name of the file that holds the visualization if save_figure is True. Defaults to None.
            animate (bool, optional): User input to whether the visualization needs to be animated or not. Defaults to False.
        
        Raises:
            ValueError: Raised when an unrecognized visualization type is entered.
        """
        try:
            for star in self.stars:
                assert len(star.trajectory)>0
        except AssertionError:
            raise ValueError("The trajectory of the celestial objects is empty. Please run the simulation first")
        celestial_objects = self.__compile_results()
        vis = Visualize(celestial_objects=celestial_objects, visualization_type=visualization_type, 
                        save_figure=save_figure, figure_name=figure_name)
        vis.visualize(animate=animate)


class BinaryMerger():
    """Class that simulates the merger of two celestial objects.

    The class allows for the simulation of the merger of two celestial objects such as Neutron Stars or Black Holes. The
    simulation is performed using the NRSur7dq2 and SurfinBH libraries for the gravitational waveforms and the final
    remnant mass and spin calculations. The simulation is performed under the assumption that the two objects are in a
    binary system and are in the process of merging. The simulation is performed in natural units.

    Attributes:
        binaryBH (module): The module that contains the functions for the Binary Merger simulation.
        binary_system (list): List of two objects that are to be simulated for Binary Merger.
        q (float): The mass ratio of the two objects in the binary system.
        chi1 (np.array): The dimensionless spin of the first object in the binary system.
        chi2 (np.array): The dimensionless spin of the second object in the binary system.
        trajectory (list): A list that contains the trajectory of the binary systems as they merge.
        waveforms (dict): The gravitational waveforms generated by the binary system (contains multiple modes).
        time_vals (np.array): The time values for the gravitational waveforms and trajectories.

    Methods:
        simulate(): Calculates the trajectory of the two objects as they merge.
        visualize(): Visualizes the trajectory of the two objects as they merge.
        gw_plot(mode): Plots the gravitational waveforms generated by the binary system.
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
            assert platform.system()=="Linux"
        except AssertionError:
            raise OSError("Binary Merger simulation is only supported on Linux")
        try:
            import utils.binaryBH as binaryBH
            self.binaryBH = binaryBH
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
        self.calculate_chi
    
    @property
    def calculate_chi(self)->None:
        """Property of the BinaryMerger class that calculates the dimensionless spins of the two objects in the binary
        system.
        """
        self.chi1 = const.C * self.binary_system[0].angular_momentum / (const.G * self.binary_system[0].mass**2)
        self.chi2 = const.C * self.binary_system[1].angular_momentum / (const.G * self.binary_system[1].mass**2)
    
    def simulate(self):
        """ Method of the BinaryMerger class that calculates the trajectory of the two objects as they merge.
        """
        vals = self.binaryBH.get_binary_data(self.q, self.chi1, self.chi2, None, None, None)
        self.trajectory = [vals[5], vals[6]]
        self.waveforms = vals[4]
        self.time_vals = vals[0]

    def animate(self):
        radius_list = [self.binary_system[0].radius, self.binary_system[1].radius]
        _animate_merger(self.time_vals, self.trajectory, radius_list)

    def gw_plot(self, mode:tuple)->None:
        """Method of the BinaryMerger class that plots the gravitational waveforms generated by the binary system.

        Args:
            mode (tuple): The mode of the gravitational waveforms to be plotted.
        """
        import matplotlib.pyplot as plt
        try:
            assert mode in self.waveforms.keys()
        except AssertionError:
            raise ValueError(f"Mode {mode} not found in the waveforms")
        hp = self.waveforms[mode].real
        hx = self.waveforms[mode].imag
        plt.plot(self.time_vals, hp, label="hp")
        plt.plot(self.time_vals, hx, label="hx")
        plt.xlabel("Time (M)")
        plt.ylabel("Strain (hr/M)")
        plt.legend()
        plt.show()