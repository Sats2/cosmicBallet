import numpy as np
from typing import Union
import math
from cosmicBallet.utils import Constants as const


class Planets():
    """Class that initializes the object Planet for the N-Body Simulation and Visualization
    
    This is the class that creates a template for Planets that will be used in the N-Body Simulation. The initialized planets are assumed
    to be perfectly spherical for ease of volume/density analysis as well as collisions. The class contains information regarding the planet 
    for simulation, such as the mass and orbital properties, as well as, for visualizations to load the contours of the planets in the animation 
    renders.

    Attributes:
        name (str): Name of the Planet
        mass (float/int): Mass of the Planet in kilograms
        radius (float/int): Radius of the Planet in meters
        volume (float): Volume of the Planet in meter^3
        density (float): Density of the Planet in kilograms/meter^3
        planet_type (str): Type of Planet. Accepted values are "Rocky" or "Gaseous"
        planet_contour (str): Planet contour type ("Earth-like"/"Mars-like" for Rocky Planets and "Jupiter-like"/"Neptune-like" for
                                Gaseous Planets)
        init_position (list): Initial position of the Planet in space as a list of coordinates in meters
        init_velocity (list): Initial orbital velocity of the Planet as a list of directional velocities in meter/second.
        position (array): Position of the planet in a given time slice.
        velocity (array): Velocity of the planet in a given time slice.
        object_type (str): Identifier for type of Celestial Object.
        force (array): Total Force acting on the Planet in a given time slice.
        momentum (array): Momentum of the Planet in a given time slice.
        material_property (dict): Aggregate Material properties of the planet based on abundant materials.
        trajectory (list): Holds the trajectory of the planet as a list of position arrays.
        color (tuple): Color of the planet for visualization purposes.

    Methods:
        radius(): Gets the value of the radius
        radius(value): Sets the value of the radius
        mass(): Gets the value of the mass
        mass(value): Sets the value of the mass
        volume(): Calculates and updates the volume of the planet
        density(): Calculates and updates the density of the planet
    """

    def __init__(self, name:str, mass:Union[float, int], radius:Union[float,int], planet_type:str, 
                 planet_contour:str, init_position:list, init_velocity:list, material_property:dict):
        """Constructor for the Planets class.

        Args:
            name (str): Name of the Planet
            mass (Union[float, int]): Mass of the Planet in kilograms
            radius (Union[float,int]): Radius of the Planet in meters
            planet_type (str): Type of the Planet. Accepted values are "Rocky" or "Gaseous"
            planet_contour (str): Contour of the Planet. Accepted values are "Earth-like", 
                                "Mars-like", "Jupiter-like" or "Neptune-like"
            init_position (list): Initial Position of the Planet as a list of coordinates (in meters)
            init_velocity (list): Initial Velocity of the Planet as a list of directional velocities (in meters/second)
            material_property (dict): Material Properties of the Planet.

        Raises:
            TypeError: name is not of type string
            TypeError: mass is not of type float/int
            TypeError: radius is not of type float/int
            TypeError: planet_type is not of type string
            TypeError: planet_contour is not of type string
            TypeError: init_position is not of type list
            TypeError: init_velocity is not of type list
            TypeError: material_property is not of type dictionary
            ValueError: name is None
            ValueError: mass is zero or negative
            ValueError: radius is zero or negative
            ValueError: planet_type is not "Rocky" or "Gaseous"
            ValueError: planet_contour is not "Earth-like", "Mars-like", "Jupiter-like" or "Neptune-like"
            ValueError: init_position does not contain 3 coordinates
            ValueError: init_velocity does not contain 3 directional velocities
            TypeError: init_position does not contain float/int values
            TypeError: init_velocity does not contain float/int values
        """
        try:
            assert isinstance(name, str)
        except AssertionError:
            raise TypeError("Planet Property 'name' can only be of type string")
        try:
            assert isinstance(mass, (float, int))
        except AssertionError:
            raise TypeError("Planet Property 'mass' can only be of type float/int")
        try:
            assert isinstance(radius, (float, int))
        except AssertionError:
            raise TypeError("Planet Property 'radius' can only be of type float/int")
        try:
            assert isinstance(planet_type, str)
        except AssertionError:
            raise TypeError("Planet Property 'planet_type' can only be of type string")
        try:
            assert isinstance(planet_contour, str)
        except AssertionError:
            raise TypeError("Planet Property 'planet_contour' can only be of type string")
        try:
            assert isinstance(init_position, list)
        except AssertionError:
            raise TypeError("Planet Property 'init_position' can only be of type list")
        try:
            assert isinstance(init_velocity, list)
        except AssertionError:
            raise TypeError("Planet Property 'init_velocity' can only be of type list")
        try:
            assert isinstance(material_property, dict)
        except AssertionError:
            raise TypeError("Planet Property 'material_property' must be of type dictionary")
        try:
            assert mass>0
        except AssertionError:
            raise ValueError("Planet Property 'mass' can not be zero or non-negative")
        try:
            assert radius>0
        except AssertionError:
            raise ValueError("Planet Property 'radius' can not be zero or non-negative")
        try:
            assert (name is not None)
        except AssertionError:
            raise ValueError("Planet Property 'name' cannot be None")
        try:
            assert (planet_type.lower() == "rocky" or planet_type.lower() == "gaseous")
        except AssertionError:
            raise ValueError("Planet Property 'planet_type' can only be 'Rocky' or 'Gaseous'")
        try:
            assert (planet_contour.lower() == "earth-like" or planet_contour.lower() == "mars-like" or \
                    planet_contour.lower() == "jupiter-like" or planet_contour.lower() == "neptune-like")
        except AssertionError:
            raise ValueError("Planet Property 'planet_contour' can only be 'Earth-like', 'Mars-like', 'Jupiter-like' or 'Neptune-like'")
        try:
            assert len(init_position)==3
        except AssertionError:
            raise ValueError("Planet Property 'init_position' must contain 3 coordinates")
        try:
            assert len(init_velocity)==3
        except AssertionError:
            raise ValueError("Planet Property 'init_velocity' must contain 3 directional velocities")
        try:
            assert all(isinstance(i, (float, int)) for i in init_position)
        except AssertionError:
            raise TypeError("Planet Property 'init_position' must contain float/int values")
        try:
            assert all(isinstance(i, (float, int)) for i in init_velocity)
        except AssertionError:
            raise TypeError("Planet Property 'init_velocity' must contain float/int values")
        self.name = name
        self.mass = mass
        self.radius = radius
        self.planet_type = planet_type
        self.planet_contour = planet_contour
        self.init_position = np.array(init_position)
        self.init_velocity = np.array(init_velocity)
        self.material_property = material_property
        self.position = None
        self.velocity = None
        self.momentum = None
        self.trajectory = []
        self.vel_list = []
        self.force = np.zeros(3)
        self.object_type = "planet"
        if self.planet_type.lower() == "rocky":
            if self.planet_contour.lower() == "earth-like":
                self.color_myv = (0,0.2,1)
                self.color = "blue"
            else:
                self.color_myv = (1,0,1)
                self.color = "brown"
        else:
            if self.planet_contour.lower() == "jupiter-like":
                self.color_myv = (0.5,1,0)
                self.color = "orange"
            else:
                self.color_myv = (0,0,1)
                self.color = "blue"
    
    @property
    def radius(self):
        """Gets the radius of the planet as a class property

        Returns:
            float: radius of the planet
        """
        return self._radius
    
    @radius.setter
    def radius(self, value:Union[float,int]):
        """Sets the radius of the Planet.

        Args:
            value (float/int): Radius of the Planet
        """
        self._radius = value

    @property
    def color(self):
        """Gets the color of the planet as a class property

        Returns:
            str: Color of the planet
        """
        return self._color
    
    @color.setter
    def color(self, value:str):
        """Sets the color of the Planet

        Args:
            value (str): Color of the Planet
        """
        self._color = value
    
    @property
    def color_myv(self):
        """Gets the color of the planet as a class property as a tuple
        
        Returns:
            tuple: Color of the planet as a tuple
        """
        return self._color_myv
    
    @color_myv.setter
    def color_myv(self, value:tuple):
        """Sets the color of the Planet as a tuple

        Args:
            value (tuple): Color of the Planet as a tuple
        """
        self._color_myv = value
    
    @property
    def mass(self):
        """Gets the mass of the planet as a class property

        Returns:
            float: mass of the planet
        """
        return self._mass
    
    @mass.setter
    def mass(self, value:float):
        """Sets the mass of the Planet

        Args:
            value (float): Mass of the Planet
        """
        self._mass = value
    
    @property
    def volume(self):
        """Dynamically Calculates the Volume of the Planet as a class property for each update to its radius

        The planet is spherical and calculated using the formula:
            volume = (4/3) * pi * radius^3

        Returns:
            float: Volume of the planet
        """
        return (4/3) * np.pi * np.power(self._radius, 3)
    
    @property
    def density(self):
        """Dynamically Calculates the Density of the Planet as a class property for each update to its mass or radius

        Calculated using the formula: density = mass / volume

        Raises:
            ValueError: In case the volume (or radius) is set as zero, a value error is raised to prevent division by zero

        Returns:
            float: Density of the Planet
        """
        try:
            assert self.volume!=0, "Volume cannot be Zero. Modify Radius"
        except AssertionError:
            raise ValueError
        return self._mass/self.volume


class Fragments():
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
        object_type (str): Defines the object as a free fragment.
        trajectory (list): List containing the trajectory of the frament.
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
        self.planet_type = "fragment"
        self.object_type = "fragment"
        self.trajectory = []
        self.color = "grey"
    

class Stars():
    """Class that initializes the object Star for the N-Body Simulation and Visualization
    
    This is the class that creates a template for Stars that will be used in the N-Body Simulation. The initialized Stars are assumed to be
    perfectly spherical for simplified density/volume analysis and collision cases. The class contains information regarding the star for 
    simulation, such as the mass and orbital properties, as well as, for visualizations to load the color, size and other relevent visual 
    properties of the stars in the animation renders. The class requires either the radius or the density as input (if not both).

    Attributes:
        name (str): Name of the Star
        mass (float/int): Mass of the Star in kilograms
        temperature (float/int): Temperature of the Star in kelvin
        init_position (list): Initial position of the star in space as a list of coordinates in meters
        init_velocity (list): Initial orbital velocity of the star as a list of directional velocities in meter/second
        radius (float/int, optional): Radius of the Star in meters
        density (float, optional): Density of the Star in kilograms/meter^3
        position (array): Position of the star in a given time slice.
        velocity (array): Velocity of the star in a given time slice.
        force (array): Total Force acting on the Planet in a given time slice.
        momentum (array): Momentum of the star in a given time slice.
        object_type (str): Identifier for type of Celestial Object.
        star_type (str): Type of Star based on density
        star_class (str): Star Classification based on temperature
        trajectory (list): Holds the trajectory of the star as a list of position arrays.
        color_myv (tuple): Color of the star for visualization purposes.
        color (str): Color of the star for visualization purposes.

    Methods:
        radius(): Gets the radius of the star
        radius(value): Sets the radius of the star
        mass(): Gets the mass of the star
        mass(value): Sets the mass of the star
        volume(): Calculates and updates the volume of the star
        density(): Calculates and updates the volume of the star
        star_type(): Determines and updates the type of the star
        star_class(): Determines and updates the class of the star
    """

    def __init__(self, name:str, mass:Union[float, int], temperature:Union[float,int], init_position:list, init_velocity:list,
                 radius:Union[float,int]=None, density:float=None, angular_momentum:list=None):
        """Initializes the Star Object

        The class requires either radius or density to be provided. The radius is calculated from the density and mass but for dynamic 
        modifications to the star properties only radius updates are accepted.

        The radius is calculated from the density and mass as follows:
            volume = mass / density
            radius = (3 * volume / (4 * pi))**(1/3)

        Args:
            name (str): Name of the Star
            mass (float/int): Mass of the Star in kilograms
            temperature (float/int): Temperature of the Star in kelvin
            init_position (list): Initial position of the star in space as a list of coordinates in meters
            init_velocity (list): Initial orbital velocity of the star as a list of directional velocities in meter/second
            radius (float/int, optional): Radius of the Star in meters
            density (float, optional): Density of the Star in kilograms/meter^3
            angular_momentum (list, optional): Angular Momentum of the Star as a list of axial momenta in kilograms*radian/second
        
        Raises:
            TypeError: name is not of type string
            TypeError: mass is not of type float/int
            TypeError: temperature is not of type float/int
            TypeError: init_position is not of type list
            TypeError: init_velocity is not of type list
            TypeError: radius is not of type float/int
            TypeError: density is not of type float
            TypeError: angular_momentum is not of type list (if provided)
            TypeError: init_position does not contain float/int values
            TypeError: init_velocity does not contain float/int values
            ValueError: name is None
            ValueError: mass is zero or negative
            ValueError: radius is zero or negative
            ValueError: temperature is zero or negative
            ValueError: density is zero or negative
            ValueError: init_position does not contain 3 coordinates
            ValueError: init_velocity does not contain 3 directional velocities
        """
        try:
            assert isinstance(name, str)
        except AssertionError:
            raise TypeError("Star Property 'name' can only be of type string")
        try:
            assert isinstance(mass, (float, int))
        except AssertionError:
            raise TypeError("Star Property 'mass' can only be of type float/int")
        try:
            assert isinstance(temperature, (float, int))
        except AssertionError:
            raise TypeError("Star Property 'temperature' can only be of type float/int")
        try:
            assert isinstance(init_position, list)
        except AssertionError:
            raise TypeError("Star Property 'init_position' can only be of type list")
        try:
            assert isinstance(init_velocity, list)
        except AssertionError:
            raise TypeError("Star Property 'init_velocity' can only be of type list")
        try:
            assert all(isinstance(i, (float, int)) for i in init_position)
        except AssertionError:
            raise TypeError("Star Property 'init_position' can only contain float/int values")
        try:
            assert all(isinstance(i, (float, int)) for i in init_velocity)
        except AssertionError:
            raise TypeError("Star Property 'init_velocity' can only contain float/int values")
        if radius is not None:
            try:
                assert isinstance(radius, (float, int))
            except AssertionError:
                raise TypeError("Star Property 'radius' can only be of type float/int")
        if density is not None:
            try:
                assert isinstance(density, float)
            except AssertionError:
                raise TypeError("Star Property 'density' must be of type float")
        if angular_momentum is not None:
            try:
                assert isinstance(angular_momentum, list)
            except AssertionError:
                raise TypeError("Star Property 'angular_momentum' must be of type list")
        try:
            assert mass>0
        except AssertionError:
            raise ValueError("Star Property 'mass' must be a positive value")
        try:
            assert temperature>0
        except AssertionError:
            raise ValueError("Star Property 'temperature' must be a positive value")
        try:
            assert (name is not None)
        except AssertionError:
            raise ValueError("Star Property 'name' cannot be None")
        try:
            assert (radius is not None or density is not None)
        except AssertionError:
            raise ValueError("Star Properties 'radius' and 'density' cannot be None")
        if radius is not None:
            try:
                assert radius>0
            except AssertionError:
                raise ValueError("Star Property 'radius' must be a positive value")
        if density is not None:
            try:
                assert density>0
            except AssertionError:
                raise ValueError("Star Property 'density' must be a positive value")
        if radius is not None and density is not None:
            vol = (4/3) * np.pi * np.power(radius, 3)
            calc_density = mass / vol
            try:
                assert math.isclose(calc_density, density, abs_tol=1e-3)
            except AssertionError:
                raise ValueError("Provided Density of Star and Calculated Density of Star do not match")
        self.name = name
        self.mass = mass
        self.temperature = temperature
        self.init_position = np.array(init_position)
        self.init_velocity = np.array(init_velocity)
        if radius is not None:
            self.radius = radius
        if density is not None:
            self.density = density
            vol = self.mass/self.density
            self.radius = np.cbrt((3/4) * vol / np.pi)
        self.position = None
        self.velocity = None
        self.momentum = None
        self.force = np.zeros(3)
        self.object_type = "star"
        self.color = ""
        self.color_myv = (1,1,0)
        self.star_class
        self.trajectory = []
        self.vel_list = []
        if angular_momentum is not None:
            self.angular_momentum = np.array(angular_momentum).astype(np.float64)
        else:
            self.angular_momentum = np.array([0, 0, 0]).astype(np.float64)
        
    @property
    def radius(self):
        """Gets the radius as a dynamic property of the class

        Returns:
            float: Radius of the Star
        """
        return self._radius
    
    @property
    def mass(self):
        """Gets the mass as a dynamic property of the class

        Returns:
            float: Mass of the Star
        """
        return self._mass
    
    @radius.setter
    def radius(self, value:Union[float,int]):
        """Sets the radius of the star

        Args:
            value (float/int): New radius of the star
        """
        self._radius = value

    @mass.setter
    def mass(self, value:Union[float,int]):
        """Sets the mass of the star

        Args:
            value (float/int): New mass of the star
        """
        self._mass = value

    @property
    def volume(self):
        """Dynamically updates the volume of the star

        The star is spherical and the volume is calculated using the formula:
            volume = (4/3) * pi * radius^3

        Returns:
            float: Volume of the star
        """
        return (4/3) * np.pi * np.power(self._radius, 3)
    
    @property
    def density(self):
        """Dynamically calculates the density of the star

        Calculated using the formula: density = mass / volume

        Raises:
            ValueError: In case the volume of the star is zero, in order to prevent division by zero error

        Returns:
            float: Density of the star
        """
        try:
            assert self.volume!=0, "Volume cannot be zero. Modify Radius"
        except AssertionError:
            raise ValueError
        return self._mass/self.volume
    
    @property
    def star_type(self):
        """Determines the type of the star based on its density.

        Classification of the stars follow as (the classification and an example is provided):
            (i)     The star is considered to be a Giant/Hypergiant if its density is less than 1000 kg/m^3 (or < 1 g/cm^3). 
                        Example: Stephanson 2-18, density ~ 70.4 kg/m^3, Red Hypergiant (largest known star)
            (ii)    The star is considered to be a Main Sequence star if its density is between 1000 and 10^8 kg/m^3 (or 1-10^5 g/cm^3)
                        Example: The Sun, density ~ 1408 kg/m^3, Main Sequence Yellow Star
            (iii)   The star is considered to be a White Dwarf (remnant of stellar death) is its density is between 10^8 and 10^16 kg/m^3 (or 10^5 - 10^13 g/cm^3)
                        Example: Sirius B, density ~ 5.7733 * 10^8 kg/m^3, Closest White Dwarf to the Sun
            (iv)    The star is considered a Neutron star if it's density exceeds 10^16 kg/m^3 (or 10^13 g/cm^3)
                        Example: PSR J0952-0607, density ~ 1.1155 * 10^18 kg/m^3, Most massive Neutron Star known as of today.
        
        For more information on the examples, search the star example on Wikipedia, or look up star density classification.
        Infinite Density is not accepted as a star property and the object must be set as a Black Hole! 

        Returns:
            str: Type of the Star
        """
        if self.density > 0 and self.density < 1000:
            return "Giant"
        elif self.density > 1000 and self.density < 1e8:
            return "Main Sequence"
        elif self.density > 1e8 and self.density < 1e16:
            self.color = "white"
            self.color_myv = (1,1,1)
            return "White Dwarf"
        else:
            self.color = "white"
            self.color_myv = (1,1,1)
            return "Neutron"
    
    @property
    def star_class(self):
        """Determines the Star Classification for visualization (purely a visual aspect for model selection)

        Star Classification based on the Harvard spectral classification of stars based on the effective temperature of the surface
        of the star. The classification of the stars work as follows (given -> Star Class, Temperature Range, Chromaticity):
            For Main Sequence Stars:
                (i)     Class M, if 2300K < Temperature < 3900K, light orangish red
                (ii)    Class K, if 3900K < Temperature < 5300K, pale yellowish orange
                (iii)   Class G, if 5300K < Temperature < 6000K, yellowish white
                (iv)    Class F, if 6000K < Temperature < 7300K, white
                (v)     Class A, if 7300K < Temperature < 10000K, bluish white
                (vi)    Class B, if 10000K < Temperature < 33000K, deep bluish white
                (vii)   Class O, if Temperature >= 33000K, blue
            For Giants/Hypergiants:
                (i)     Red, if 3700K < Temperature < 10000K, red giant/hypergiant
                (ii)    Blue, if Temperature >= 10000K, blue giant/hypergiant
            White dwarfs and Neutron Stars do not have spectral classification.
            The color is based on the chromaticity of the star and is used for visualization purposes.

        Raises:
            ValueError: If the surface temperature of the star is too low for luminosity.

        Returns:
            str: Classification of star based on luminosity (or effective surface temperature)
        """
        if self.star_type == "Main Sequence":
            if self.temperature > 2300 and self.temperature <= 3900:
                self.color_myv = (1, 0, 0)
                self.color = "red"
                return "M"
            elif self.temperature > 3900 and self.temperature <= 5300:
                self.color_myv = (1,0.5,0)
                self.color = "orange"
                return "K"
            elif self.temperature > 5300 and self.temperature <= 6000:
                self.color_myv = (1,1,0)
                self.color = "yellow"
                return "G"
            elif self.temperature > 6000 and self.temperature <= 7300:
                self.color_myv = (0.8,0.9,1)
                self.color = "white"
                return "F"
            elif self.temperature > 7300 and self.temperature <= 10000:
                self.color_myv = (1,1,1)
                self.color = "white"
                return "A"
            elif self.temperature > 10000 and self.temperature <= 33000:
                self.color_myv = (0,0.3,1)
                self.color = "purple"
                return "B"
            elif self.temperature > 33000:
                self.color_myv = (0,0,1)
                self.color = "blue"
                return "O"
            else:
                raise ValueError("Temperature of Star Too Low and uncharacteristic of Stars. Modification to Temperature needed")
        elif self.star_type == "Giant":
            if self.temperature > 3700 and self.temperature < 10000:
                self.color_myv = (1,0,0)
                self.color = "red"
                return "Red"
            elif self.temperature > 10000:
                self.color_myv = (0,0,1)
                self.color = "blue"
                return "Blue"
            else:
                raise ValueError("Temperature of Star Too Low and uncharacteristic of Giants. Modification to Temperature needed")


class BlackHole():
    """Class that initializes a Black Hole.

        The characteristics of the black hole is calculated within the class dynamically. If the angular momentum of the black
        hole is not specified, it is assumed that the black is a non-rotating black hole.

    Attributes:
        name (str): Name of the Black Hole
        mass (float/int): Mass of the Black Hole in kilograms
        init_position (np.array): Initial position of the Black Hole in space as a list of coordinates in meters
        init_velocity (np.array): Initial orbital velocity of the Black Hole as a list of directional velocities in meter/seconds
        angular_momentum (list, optional): Angular momentum of a rotating Black Hole as list. Ignore for non-rotating Black Holes.
        radius (float): Schwarzchild radius of the Black Hole
        spin (np.array): Dimensionless spin coefficients of the Black Hole
        position (np.array): Holds the position of the Black Hole at each time step
        velocity (np.array): Holds the velocity of the Black Hole at each time step
        trajectory (list): Contains the trajectory of the Black Hole.
        force (np.array): Attribute that holds the force acting on the Black Hole at each time step.
        color (str): Color of the Black Hole for visualization purposes.
        color_myv (tuple): Color of the Black Hole for visualization purposes.
    
    Methods:
        mass(): Gets the mass of the Black Hole
        mass(value): Sets the mass of the Black Hole
        radius(): Calculates and updates the Schwarzchild radius of the Black Hole
        spin(): Calculates the spin of the Black Hole based on the angular momentum
    """

    def __init__(self, name:str, mass:Union[float,int], init_position:list, init_velocity:list, angular_momentum:list=None):
        """Constructor for the Black Hole Class

        Args:
            name (str): Name of the Black Hole
            mass (Union[float,int]): Mass of the Black Hole in kilograms
            init_position (list): Initial position of the Black Hole in space as a list of coordinates in meters
            init_velocity (list): Initial orbital velocity of the Black Hole as a list of directional velocities in meter/seconds
            angular_momentum (list, optional): Angular Momentum of the Black Hole. Defaults to None.

        Raises:
            TypeError: name is not of type string
            TypeError: mass is not of type float/int
            TypeError: init_position is not of type list
            TypeError: init_velocity is not of type list
            TypeError: angular_momentum is not of type list (if provided)
            TypeError: init_position does not contain float/int values
            TypeError: init_velocity does not contain float/int values
            TypeError: angular_momentum does not contain float/int values (if provided)
            ValueError: name is None
            ValueError: mass is zero or negative
        """
        try:
            assert isinstance(name, str)
        except AssertionError:
            raise TypeError("Black Hole property 'name' must be of type string")
        try:
            assert isinstance(mass, (float, int))
        except AssertionError:
            raise TypeError("Black Hole property 'mass' must be of type float/int")
        try:
            assert isinstance(init_position, list)
        except AssertionError:
            raise TypeError("Black Hole property 'init_position' must be of type list")
        try:
            assert isinstance(init_velocity, list)
        except AssertionError:
            raise TypeError("Black Hole property 'init_velocity' must be of type list")
        try:
            assert all(isinstance(i, (float, int)) for i in init_position)
        except AssertionError:
            raise TypeError("Black Hole Property 'init_position' can only contain float/int values")
        try:
            assert all(isinstance(i, (float, int)) for i in init_velocity)
        except AssertionError:
            raise TypeError("Black Hole Property 'init_velocity' can only contain float/int values")
        if angular_momentum is not None:
            try:
                assert isinstance(angular_momentum, list)
            except AssertionError:
                raise TypeError("Black Hole Property 'angular_momentum' must be of type list")
            try:
                assert all(isinstance(i, (float, int)) for i in angular_momentum)
            except AssertionError:
                raise TypeError("Black Hole Property 'angular_momentum' can only contain float/int values")
        try:
            assert mass>0
        except AssertionError:
            raise ValueError("Black Hole Property 'mass' must be a positive value")
        try:
            assert (name is not None)
        except AssertionError:
            raise ValueError("Black Hole Property 'name' cannot be None")
        self.name = name
        self.mass = mass
        self.init_position = np.array(init_position)
        self.init_velocity = np.array(init_velocity)
        self.position = np.zeros(3)
        self.velocity = np.zeros(3)
        self.trajectory = []
        self.object_type = "black_hole"
        if angular_momentum is not None:
            self.angular_momentum = np.array(angular_momentum)
        else:
            self.angular_momentum = np.array([0, 0, 0]).astype(np.float64)
        self.force = np.zeros(3)
        self.color = "black"
        self.color_myv = (0,0,0)
    
    @property
    def mass(self)->Union[float,int]:
        """Initializes the mass of the black hole as a dynamic property.

        Returns:
            float/int: Mass of the Black Hole
        """
        return self._mass
    
    @mass.setter
    def mass(self, value:Union[float,int])->None:
        """Sets the mass of the Black Hole.

        Args:
            value (float/int): New mass of the black hole
        """
        self._mass = value

    @property
    def radius(self)->float:
        """Calculates the Schwarzchild radius of the Black Hole.

        The Schwarzchild radius of the black hole (in meters) is the distance between the singularity (center) and the edge 
        of the event horizon (visible part of the black hole).
        Calculated using the formula:
                radius_schwarzchild = 2 * G * M / c^2
                    where:
                        G - Universal Gravitational Constant (Newton*meter^2/kilogram^2)
                        M - Mass of the Black Hole (singularity) (kilogram)
                        c - Speed of light (meter/second)

        Returns:
            float: Schwarzchild radius of the Black Hole.
        """
        radius = 2 * const.G * self._mass / (const.C * const.C)
        return radius

    @property
    def spin(self)->np.array:
        """Calculates the dimensionless spin coefficients of the Black Hole based on the angular momentum.

        The spin of the black hole is calculated based on the angular momentum of the black hole. The spin is calculated as:
            spin = angular_momentum * c / (G * m^2)
                where:
                    angular_momentum - Angular Momentum of the Black Hole (kilogram*meter^2/second)
                    mass - Mass of the Black Hole (kilogram)
                    c - Speed of light (meter/second)
                    G - Universal Gravitational Constant (Newton*meter^2/kilogram^2)

        Returns:
            array: Spin of the Black Hole
        """
        if all(self.angular_momentum==0):
            return np.array([0, 0, 0]).astype(np.float64)
        spin = self.angular_momentum * const.C / (const.G * np.power(self._mass, 2))
        return spin
