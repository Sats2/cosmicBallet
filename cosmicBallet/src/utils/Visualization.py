import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from typing import Union
import numpy as np
import mayavi.mlab as mlab


class Visualize():
    """A class for Visualizing the trajectory of the orbits of the N-Bodies Simulated

    Attributes:
        celestial_objects (object): An object belonging to Simulation Classes containing the simulation results for 
                                    the orbital trajectory.
        visual_type (str): A user input to determine whether to render an animation, or a simple scientific plot
        save_fig (bool): A user input to determine whether to save the trajectory visualization. Defaults to False
        figure_name (str): The name for the figure to be saved. Defaults to None
        ani_name (str): The name for the animation to be saved. Defaults to None
    """
    def __init__(self, celestial_objects:list, visualization_type:str="scientific", save_figure:bool=False,
                 figure_name:str=None) -> None:
        """Constructor for the Visualization class

        Args:
            celestial_objects (object): An object of the simulation class containing the trajectory of the bodies
            visualization_type (str, optional): Type of Visualization (scientific or animation). Defaults to "scientific".
            save_figure (bool, optional): A value to determine whether to save the visualization or not. Defaults to False.
            figure_name (str, optional): Name of the file that will hold the visualization if save_figure is True. Defaults to None.

        Raises:
            TypeError: Raised when input values do not belong to the required datatypes
            ValueError: Raised when visualization_type does not belong to the allowed value.
        """
        try:
            assert isinstance(celestial_objects, object), "The celestial_objects needs to be a list of objects for visualization"
            assert isinstance(visualization_type, str), "The visualization_type needs to be a string"
            assert isinstance(save_figure, bool), "The attribute save_figure must be a boolean"
            assert isinstance(figure_name, str), "The figure_name must be a string"
        except AssertionError:
            raise TypeError
        try:
            accepted_choices = ["scientific", "animation"]
            assert (visualization_type.lower() in accepted_choices), "The visualization type must be either scientific/animation"
        except AssertionError:
            raise ValueError
        self.celestial_objects = celestial_objects
        self.visual_type = visualization_type
        self.save_fig = save_figure
        self.ani_name = None
        if figure_name is not None and visualization_type == "scientific":
            self.figure_name = figure_name + ".jpeg"
            self.ani_name = figure_name + ".gif"
        else:
            if self.figure_name is not None:
                self.ani_name = figure_name + ".gif"

    def __create_matplotlib_animation(self):
        """Private method of the Visualization class that generates an animation of the trajectory of the bodies.	
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        points = []
        traj = []
        all_points = []
        for body in self.celestial_objects:
            spacetime_point = body.trajectory[0]
            if spacetime_point[0] == 0:
                points.append(ax.plot(spacetime_point[1], spacetime_point[2], spacetime_point[3], "o", label=body.name)[0])
                traj.append(ax.plot(spacetime_point[1], spacetime_point[2], spacetime_point[3])[0])
                all_points.append(spacetime_point)
            else:
                points.append(None)
                traj.append(None)
        all_points = np.array(all_points)
        ax.set_xlim(all_points[:,1].min(), all_points[:,1].max())
        ax.set_ylim(all_points[:,2].min(), all_points[:,2].max())
        ax.set_zlim(all_points[:,3].min(), all_points[:,3].max())
        ax.set_xlabel("X [m]")
        ax.set_ylabel("Y [m]")
        ax.set_zlabel("Z [m]")

        def update(frame):
            all_points = []
            for body in self.celestial_objects:
                spacetime_point = body.trajectory[frame]
                if spacetime_point[0] == frame:
                    i = self.celestial_objects.index(body)
                    points[i].set_data(spacetime_point[1], spacetime_point[2])
                    points[i].set_3d_properties(spacetime_point[3])
                    spacetime_trajectory = np.array(body.trajectory)[:frame]
                    traj[i].set_data(spacetime_trajectory[:, 1], spacetime_trajectory[:, 2])
                    traj[i].set_3d_properties(spacetime_trajectory[:, 3])
                    all_points.append(spacetime_trajectory)
            all_points = np.array(all_points)
            ax.set_xlim(all_points[:,1].min(), all_points[:,1].max())
            ax.set_ylim(all_points[:,2].min(), all_points[:,2].max())
            ax.set_zlim(all_points[:,3].min(), all_points[:,3].max())
            return traj+points
        
        ani = FuncAnimation(fig, update, frames=len(self.celestial_objects[0].trajectory), blit=False)
        if self.save_fig:
            ani.save(self.ani_name, dpi=300)
        plt.show()
            
    
    def __scientific_plot(self, animate:bool)->None:
        """Private method of the Visualization class that generates a maplotlib plot of the trajectory of the bodies.

        Raises:
            ValueError: When no trajectory for the celestial objects are found
        """
        try:
            for body in self.celestial_objects:
                assert body.trajectory!=[], f"Empty trajectory found for {body.name}. Run Simulations first"
        except AssertionError:
            raise ValueError
        if animate:
            self.__create_matplotlib_animation()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        i = 1
        for body in self.celestial_objects:
            trajectory = np.array(body.trajectory)
            if body.object_type != "fragment":
                ax.plot(trajectory[:,1], trajectory[:,2], trajectory[:,3], label=body.name)
            else:
                ax.plot(trajectory[:,1], trajectory[:,2], trajectory[:,3], label=f"fragment{i}")
                i += 1
        plt.legend()
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        if self.save_fig:
            plt.savefig(self.figure_name, dpi=300)
        plt.show()
    
    def __animation(self)->None:
        fig = mlab.figure(size=(1280, 720))
        points = []
        traj = []
        for body in self.celestial_objects:
            spacetime_point = body.trajectory[0]
            if spacetime_point[0] == 0:
                points.append(mlab.points3d(spacetime_point[1], spacetime_point[2], spacetime_point[3], 
                                            scale_factor=body.radius, color=body.color))
                traj.append(mlab.plot3d(spacetime_point[:1,1], spacetime_point[:1,2], spacetime_point[:1,3],
                                        color=body.color))
            else:
                points.append(None)
                traj.append(None)
        
        def update(frame):
            spacetime_point = body.trajectory[frame]
            all_points = []
            for i,body in self.celestial_objects:
                if spacetime_point[0] == frame:
                    spacetime_trajectory = np.array(body.trajectory)[:frame]
                    points[i] = mlab.points3d(spacetime_point[1], spacetime_point[2], spacetime_point[3],
                                            scale_factor=body.radius, color=body.color)
                    traj[i] = mlab.plot3d(spacetime_trajectory[:frame,1], spacetime_trajectory[:frame,2], 
                                          spacetime_trajectory[:frame,3], color=body.color)
                    all_points.append(spacetime_trajectory)
            all_points = np.array(all_points)
            mlab.view(focalpoint=(all_points[:,1].mean(), all_points[:,2].mean(), all_points[:,3].mean()))
        
        def animate():
            for frame in range(self.n_total):
                update(frame)
                mlab.process_ui_events()
                mlab.draw()
                mlab.pause(0.05)
        
        animate()
        mlab.show()
    
    def visualize(self, animate:bool=False, time_interval:Union[float,int]=None)->None:
        """Method of the visualization class that generates the visualization for the trajectories based on the visualization type

        Args:
            animate (bool): Used for scientific visualization and control whether an animation is generated.
            time_interval (float/int, optional): Required on for animations to determine the time interval between each frame. 
                                                Defaults to None.

        Raises:
            ValueError: When the time_interval is not specified for the animation render.
        """
        try:
            assert isinstance(animate, bool), "animate can only be set to True/False"
            #assert (isinstance(time_interval, (float,int)) and self.visual_type == "animation"), "Time Interval cannot be None for Animation Renders"
        except AssertionError:
            raise ValueError
        if self.visual_type == "scientific":
            self.__scientific_plot(animate=animate)
        else:
            self.__animation()