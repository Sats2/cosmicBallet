import numpy as np
from typing import Union
import sys
import math


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
            print("Planet Initialization Failed")
            sys.exit()
        try:
            assert mass>0, "Planet Property 'mass' can not be zero or non-negative"
            assert radius>0, "Planet Property 'radius' can not be zero or non-negative"
        except AssertionError:
            print("Planet Initialization Failed")
            sys.exit()
        self.mass = mass
        self.radius = radius
        self.planet_type = planet_type
        self.planet_contour = planet_contour
        self.init_position = init_position
        self.init_velocity = init_velocity
    
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
            print("Star Initialization Failed")
            sys.exit()
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
            print("Star Initialization Failed")
            sys.exit()
        self.mass = mass
        self.temperature = temperature
        self.init_position = init_position
        self.init_velocity = init_velocity
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
        return self._mass/self.volume
    