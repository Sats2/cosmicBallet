import matplotlib.pyplot as plt
import os
import shutil
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
            TypeError: When celestial_objects is not a list of objects
            TypeError: When visualization_type is not a string
            TypeError: When save_figure is not a boolean
            TypeError: When figure_name is not a string
            ValueError: When visualization_type is not scientific/animation.
        """
        try:
            assert isinstance(celestial_objects, object)
        except AssertionError:
            raise TypeError("celestial_objects needs to be a list of objects for visualization")
        try:
            assert isinstance(visualization_type, str)
        except AssertionError:
            raise TypeError("visualization_type needs to be a string")
        try:
            assert isinstance(save_figure, bool)
        except AssertionError:
            raise TypeError("The attribute save_figure must be a boolean")
        if figure_name is not None:
            try:
                assert isinstance(figure_name, str)
            except AssertionError:
                raise TypeError("The figure_name must be a string")
        try:
            accepted_choices = ["scientific", "animation"]
            assert (visualization_type.lower() in accepted_choices)
        except AssertionError:
            raise ValueError("The visualization type must be either scientific/animation")
        self.celestial_objects = celestial_objects
        self.visual_type = visualization_type
        self.save_fig = save_figure
        self.ani_name = None
        self.n_total = int(0.1*len(self.celestial_objects[0].trajectory))
        if figure_name is not None and visualization_type == "scientific":
            self.figure_name = figure_name + ".jpeg"
            self.ani_name = figure_name + ".mp4"
        else:
            if figure_name is not None:
                self.ani_name = figure_name + ".mp4"

    def _create_figure(self,celestial_objects:list, var:int, zero_padding:str)->None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        for i in range(len(celestial_objects)):
            body = celestial_objects[i]
            trajectory = np.array(body.trajectory[:int(10*var) + 1])
            ax.plot(trajectory[:, 1], trajectory[:, 2], trajectory[:, 3], label=body.name, color=body.color)
            ax.plot(trajectory[-1, 1], trajectory[-1, 2], trajectory[-1, 3], "o", color=body.color)
        ax.legend(loc='center left', bbox_to_anchor=(1.09, 0.35))
        plt.savefig(f"temp/data_{var:0{zero_padding}d}.png", dpi=300)
        #plt.show()
        plt.clf()
        plt.close(fig)

    def __create_matplotlib_animation(self):
        """Private method of the Visualization class that generates an animation of the trajectory of the bodies.	
        """
        try:
            os.mkdir("temp")
        except:
            pass
        zero_padding = len(str(self.n_total))
        for j in range(0, self.n_total):
            var = j
            self._create_figure(self.celestial_objects, var, zero_padding)
        ffmpeg_command = (
            f"ffmpeg -r {15} -i temp/data_%0{zero_padding}d.png "
            f"-vcodec mpeg4 -qscale:v 2 -filter:v 'setpts={1/1}*PTS' -y {self.ani_name}"
        )
        os.system(ffmpeg_command)
        shutil.rmtree("temp")
            
    
    def __scientific_plot(self, animate:bool)->None:
        """Private method of the Visualization class that generates a maplotlib plot of the trajectory of the bodies.

        Raises:
            ValueError: When no trajectory for the celestial objects are found
            ValueError: When no figure name is provided to save the animation
        """
        try:
            for body in self.celestial_objects:
                assert body.trajectory!=[]
        except AssertionError:
            raise ValueError(f"Empty trajectory found for {body.name}. Run Simulations first")
        if animate:
            if self.ani_name is None:
                raise ValueError("Please provide a figure name to save the animation. Animations are saved as videos" 
                                 " and require a name")
            self.__create_matplotlib_animation()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        i = 0
        for body in self.celestial_objects:
            trajectory = np.array(body.trajectory)
            ax.plot(trajectory[:,1], trajectory[:,2], trajectory[:,3], label=body.name, color=body.color)
            ax.scatter(trajectory[-1,1], trajectory[-1,2], trajectory[-1,3], color=body.color)
            i += 1
        plt.legend()
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        if self.save_fig:
            plt.savefig(self.figure_name, dpi=300)
        plt.show()
    
    def __animation(self)->None:
        print(self.celestial_objects)
        fig = mlab.figure(size=(1280, 720))
        points = []
        traj = []
        for body in self.celestial_objects:
            spacetime_point = body.trajectory[0]
            spacetime_trajectory = np.array(body.trajectory[:1])
            if spacetime_point[0] == 1:
                points.append(mlab.points3d(spacetime_point[1], spacetime_point[2], spacetime_point[3], 
                                            scale_factor=body.radius, color=body.color_myv, resolution=20))
                traj.append(mlab.plot3d(spacetime_trajectory[0,1], spacetime_trajectory[0,2], 
                                        spacetime_trajectory[0,3], color=body.color_myv))
            else:
                points.append(None)
                traj.append(None)
        
        def update(frame):
            all_points = []
            for point in points:
                if point is not None:
                    point.remove()
            for i,body in enumerate(self.celestial_objects):
                spacetime_point = body.trajectory[frame]
                if int(spacetime_point[0]) == frame+1:
                    spacetime_trajectory = np.array(body.trajectory[:frame+1])
                    points[i] = mlab.points3d(spacetime_point[1], spacetime_point[2], spacetime_point[3],
                                            scale_factor=body.radius, color=body.color_myv, resolution=20)
                    traj[i] = mlab.plot3d(spacetime_trajectory[:,1], spacetime_trajectory[:,2], 
                                           spacetime_trajectory[:,3], color=body.color_myv, tube_radius=50)
                    all_points.append(spacetime_trajectory)
            all_points = np.array(all_points)
            mlab.view(focalpoint=(all_points[:,0,1].mean(), all_points[:,0,2].mean(), all_points[:,0,3].mean()))

        if self.save_fig:
            try:
                os.mkdir("temp")
            except:
                pass
        
        zero_padding = len(str(self.n_total))

        for frame in range(self.n_total):
            update(frame*10)
            mlab.process_ui_events()
            mlab.draw()
            if self.save_fig:
                mlab.savefig(f"temp/data_{frame:0{zero_padding}d}.png")
        
        mlab.show()

        if self.save_fig:
            ffmpeg_command = (
                            f"ffmpeg -r {15} -i temp/data_%0{zero_padding}d.png "
                            f"-vcodec mpeg4 -qscale:v 2 -filter:v 'setpts={1/1}*PTS' -y {self.ani_name}"
                            )
            os.system(ffmpeg_command)
            shutil.rmtree("temp")
    
    def visualize(self, animate:bool=False)->None:
        """Method of the visualization class that generates the visualization for the trajectories based on the visualization type

        Args:
            animate (bool): Used for scientific visualization and control whether an animation is generated. Defaults to False
            dense_body_sim (bool): Used for scientific visualization and control whether the visualization is for a 
                                dense body simulation. Defaults to False

        Raises:
            TypeError: When the animate not a boolean value.
        """
        try:
            assert isinstance(animate, bool)
        except AssertionError:
            raise TypeError("animate must be a boolean value")
        if self.visual_type == "scientific":
            self.__scientific_plot(animate=animate)
        else:
            self.__animation()    