import numpy as np
from typing import Union
import sys
import math
import Constants as const

# TODO: Add Docstrings.

class Planets():
    def __init__(self, mass:Union[float, int], radius:Union[float,int], planet_type:str, 
                 planet_contour:str, init_position:list, init_velocity:list):
        try:
            assert isinstance(mass, (float, int)), "Planet Property 'mass' can only be of type float/int"
            assert isinstance(radius, (float, int)), "Planet Property 'radius' can only be of type float/int"
            assert isinstance(planet_type, str), "Planet Property 'planet_type' can only be of type string"
            assert isinstance(planet_contour, str), "Planet Property 'planet_contour' can only be of type string"
            assert isinstance(init_position, list), "Planet Property 'init_position' can only be of type list"
            assert isinstance(init_velocity, list), "Planet Property 'init_velocity' can only be of type list"
        except AssertionError:
            raise TypeError
        try:
            assert mass>0, "Planet Property 'mass' can not be zero or non-negative"
            assert radius>0, "Planet Property 'radius' can not be zero or non-negative"
        except AssertionError:
            raise ValueError
        self.mass = mass
        self.radius = radius
        self.planet_type = planet_type
        self.planet_contour = planet_contour
        self.init_position = np.array(init_position)
        self.init_velocity = np.array(init_velocity)
    
    @property
    def radius(self):
        return self._radius
    
    @radius.setter
    def radius(self, value):
        self._radius = value
    
    @property
    def mass(self):
        return self._mass
    
    @mass.setter
    def mass(self, value):
        self._mass = value
    
    @property
    def volume(self):
        return (4/3) * np.pi * np.power(self._radius, 3)
    
    @property
    def density(self):
        try:
            assert self.volume!=0, "Volume cannot be Zero. Modify Radius"
        except AssertionError:
            raise ValueError
        return self._mass/self.volume
    


class Stars():
    def __init__(self, mass:Union[float, int], temperature:Union[float,int], init_position:list, init_velocity:list,
                 radius:Union[float,int]=None, density:float=None):
        try:
            assert isinstance(mass, (float, int)), "Star Property 'mass' can only be of type float/int"
            assert isinstance(temperature, (float, int)), "Star Property 'temperature' can only be of type float/int"
            assert isinstance(init_position, list), "Star Property 'init_position' can only be of type list"
            assert isinstance(init_velocity, list), "Star Property 'init_velocity' can only be of type list"
            if radius is not None:
                assert isinstance(radius, (float, int)), "Star Property 'radius' can only be of type float/int"
            if density is not None:
                assert isinstance(density, float), "Star Property 'density' must be of type float"
            assert (radius is not None and density is not None), "Star Properties 'radius' and 'density' cannot be None"
        except AssertionError:
            raise TypeError
        try:
            assert mass>0, "Star Property 'mass' must be a positive value"
            assert temperature>0, "Star Property 'temperature' must be a positive value"
            if radius is not None:
                assert radius>0, "Star Property 'radius' must be a positive value"
            if density is not None:
                assert density>0, "Star Property 'density' must be a positive value"
            if radius is not None and density is not None:
                vol = (4/3) * np.pi * np.power(radius, 3)
                calc_density = mass / vol
                assert math.isclose(calc_density, density, abs_tol=1e-8), "Provided Density of Star and Calculated Density of Star do not match"
        except AssertionError:
            raise ValueError
        self.mass = mass
        self.temperature = temperature
        self.init_position = np.array(init_position)
        self.init_velocity = np.array(init_velocity)
        if radius is not None:
            self.radius = radius
        if density is not None:
            self.density = density
        
    @property
    def radius(self):
        return self._radius
    
    @property
    def mass(self):
        return self._mass
    
    @radius.setter
    def radius(self, value):
        self._radius = value

    @mass.setter
    def mass(self, value):
        self._mass = value

    @property
    def volume(self):
        return (4/3) * np.pi * np.power(self._radius, 3)
    
    @property
    def density(self):
        try:
            assert self.volume!=0, "Volume cannot be zero. Modify Radius"
        except AssertionError:
            raise ValueError
        return self._mass/self.volume
    
    @property
    def star_type(self):
        if self.density > 0 and self.density < 1000:
            return "Giant"
        elif self.density > 1000 and self.density < 1e8:
            return "Main Sequence"
        elif self.density > 1e8 and self.density < 1e16:
            return "White Dwarf"
        else:
            return "Neutron"
    
    @property
    def star_class(self):
        if self.star_type == "Main Sequence":
            if self.temperature > 2300 and self.temperature < 3900:
                return "M"
            elif self.temperature > 3900 and self.temperature < 5300:
                return "K"
            elif self.temperature > 5300 and self.temperature < 6000:
                return "G"
            elif self.temperature > 6000 and self.temperature < 7300:
                return "F"
            elif self.temperature > 7300 and self.temperature < 10000:
                return "A"
            elif self.temperature > 10000 and self.temperature < 33000:
                return "B"
            elif self.temperature > 33000:
                return "O"
            else:
                raise ValueError("Temperature of Star Too Low and uncharacteristic of Stars. Modification to Temperature needed")
        elif self.star_type == "Giant":
            if self.temperature > 3700 and self.temperature < 10000:
                return "Red"
            elif self.temperature > 10000:
                return "Blue"
            else:
                raise ValueError("Temperature of Star Too Low and uncharacteristic of Giants. Modification to Temperature needed")
        
    


class Galaxy():
    def __init__(self, mass:float, radius:float, black_hole_mass:float, star_number:int,
                 init_position:list, init_velocity:list):
        try:
            assert isinstance(mass, float), "Galaxy Property 'mass' needs to be of type float"
            assert isinstance(radius, float), "Galaxy Property 'radius' needs to be of type float"
            assert isinstance(black_hole_mass, float), "Galaxy Property 'black_hole_mass' needs to be of type float"
            assert isinstance(star_number, int), "Galaxy Property 'star_number' needs to be of type int"
            assert isinstance(init_position, list), "Galaxy Property 'init_position' must be a list"
            assert isinstance(init_velocity, list), "Galaxy Property 'init_velocity' must be a list"
        except AssertionError:
            print("Galaxy Initialization failed")
            sys.exit()
        try:
            assert mass>0, "Mass must be a positive value"
            assert radius>0, "Radius needs to be a positive value"
            assert (black_hole_mass>0 and black_hole_mass<mass), "Black Hole Mass needs to be a positive value less than the Galaxy mass"
            assert star_number>0, "Number of Stars need to be a positive integer"
        except AssertionError:
            print("Galaxy Initialization failed")
            sys.exit()
        self.mass = mass
        self.radius = radius
        self.black_hole_mass = black_hole_mass
        self.star_number = star_number
        self.init_position = np.array(init_position)
        self.init_velocity = np.array(init_velocity)
    
    def create_galaxy(self):
        left_over_mass = self.mass - self.black_hole_mass
        self.star_mass_list = (np.ones(self.star_number) + np.random.normal(0,0.2,self.star_number)) \
                                * left_over_mass / self.star_number
        self.star_positions = self.init_position + np.random.normal(0,0.45,(self.star_number,2)) * self.radius
        self.star_velocity = np.full((self.star_number,2), self.init_velocity)


class BlackHole():
    def __init__(self, mass:Union[float,int], init_position:list, init_velocity:list):
        try:
            assert isinstance(mass, (float, int)), "Black Hole property 'mass' must be of type float/int"
            assert(init_position, list), "Black Hole property 'init_position' must be of type list"
            assert(init_velocity, list), "Black Hole property 'init_velocity' must be of type list"
        except AssertionError:
            raise TypeError
        try:
            assert mass>0, "Black Hole property 'mass' needs to be a positive value"
        except AssertionError:
            raise ValueError
        self.mass = mass
        self.init_position = np.array(init_position)
        self.init_velocity = np.array(init_velocity)
    
    @property
    def mass(self):
        return self._mass
    
    @mass.setter
    def mass(self, value):
        self._mass = value

    @property
    def radius(self):
        radius = 2 * const.G * self._mass / (const.c * const.c)
        return radius

