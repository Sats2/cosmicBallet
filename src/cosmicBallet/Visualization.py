import Constants as const
import matplotlib.pyplot as plt
from typing import Union
import numpy as np


class Visualize():
    """A class for Visualizing the trajectory of the orbits of the N-Bodies Simulated

    Attributes:
        results (object): An object belonging to one of the Simulation Classes containing the simulation results for the orbital trajectory
        visual_type (str): A user input to determine whether to render an animation, or a simple scientific plot
        save_fig (bool): A user input to determine whether to save the trajectory visualization. Defaults to False
        figure_name (str): The name for the figure to be saved. Defaults to None
    """
    def __init__(self, simulation_result:object, visualization_type:str="scientific", save_figure:bool=False,
                 figure_name:str=None) -> None:
        """Constructor for the Visualization class

        Args:
            simulation_result (object): An object of the simulation class containing the trajectory of the bodies
            visualization_type (str, optional): Type of Visualization (scientific or animation). Defaults to "scientific".
            save_figure (bool, optional): A value to determine whether to save the visualization or not. Defaults to False.
            figure_name (str, optional): Name of the file that will hold the visualization if save_figure is True. Defaults to None.

        Raises:
            TypeError: Raised when input values do not belong to the required datatypes
            ValueError: Raised when visualization_type does not belong to the allowed value.
        """
        try:
            assert isinstance(simulation_result, object), "The simulation result needs to be an object of one of the Simulation classes"
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
        self.results = simulation_result
        self.visual_type = visualization_type
        self.save_fig = save_figure
        self.figure_name = figure_name
    
    def __scientific_plot(self)->None:
        """Private method of the Visualization class that generates a maplotlib plot of the trajectory of the bodies.

        Raises:
            ValueError: When no trajectory for the celestial objects are found
        """
        try:
            assert (self.results.positions is not None), "No trajectory found. Please run the simulation before visualization"
        except AssertionError:
            raise ValueError
        solution = np.array(self.results.positions)
        fig = plt.figure(figsize=(20,14))
        ax = fig.add_subplot(111, projection="3d")
        for i in range(len(self.results.celestial_bodies)):
            body_result = solution[i]
            body_name = self.results.celestial_bodies.name
            ax.plot(body_result[:,0], body_result[:,1], body_result[:,2], label=body_name)
        plt.legend()
        plt.xlabel("x [m]")
        plt.ylabel("y [m]")
        plt.zlabel("z [m]")
        plt.title(f"Orbital Trajectory for {len(self.results.celestial_bodies)} Bodies for Time [0, {self.results.simulation_time}]")
        if self.save_fig:
            plt.savefig(self.figure_name, dpi=300)
        plt.show()
    
    def __animation(self)->None:
        pass
    
    def visualize(self, time_interval:Union[float,int]=None)->None:
        """Method of the visualization class that generates the visualization for the trajectories based on the visualization type

        Args:
            time_interval (float/int, optional): Required on for animations to determine the time interval between each frame. Defaults to None.

        Raises:
            ValueError: When the time_interval is not specified for the animation render.
        """
        try:
            assert (isinstance(time_interval, (float,int)) and self.visual_type == "animation"), "Time Interval cannot be None for Animation Renders"
        except AssertionError:
            raise ValueError
        if self.visual_type == "scientific":
            self.__scientific_plot()
        else:
            self.__animation()