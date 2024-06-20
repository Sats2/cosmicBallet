from Constants import material_property_list

class MaterialProperty():
    """A class that is used to generate a Planet's material properties based on the abundent elements present within.

    The class uses a weighted aggregate approach to determine a planet's Poisson Ratio, Young's Modulus and Yield Strength
    that is used for collision handling where fragmentation and shock waves require material properties. The most abundent
    materials present in planets (or moons) are available to choose from. The list of materials are given in the list:
        [silicates, iron_nickel, water_ice, ammonia_ice, methane_ice]
    
    Attributes:
        material_list (list): List of Materials in the planet. Must be available in the library materials.
        material_fraction (list): List of the fractions of the materials in the planet. Sum of all elements must be 1
    
    Methods:
        planet_material_property(): Calculates the aggregated material properties of the planets and returns it as a list
    """

    def __init__(self, material_list:list, material_fraction:list)->None:
        """Constructor of the MaterialProperty Class.

        Args:
            material_list (list): List of Materials in the planet.
            material_fraction (list): List of fractions of the materials in the planet.

        Raises:
            TypeError: Raised when the datatype of the input arguements are not of the supported types
            ValueError: Raised when there is value mismatch or iinconsistency in the number of materials.
        """
        try:
            assert isinstance(material_list, str), "Material List must be a list of material names"
            assert isinstance(material_fraction, str), "Material Fraction must be a list of material fractions"
            for item in material_list:
                assert isinstance(item, str), "Each item in material_list must be of type string"
            sum_fraction = 0
            for item in material_fraction:
                assert isinstance(item, (float,int)), "Each item in material_fraction must be of type float"
                sum_fraction += item
        except AssertionError:
            raise TypeError
        try:
            assert (len(material_list) == len(material_fraction)), "Lists material_list and material_fraction must be of same length"
            assert (sum_fraction == float(1)), "The sum of all items in material_fraction must be 1"
        except AssertionError:
            raise ValueError
        self.material_list = material_list
        self.material_fraction = material_fraction
    
    def planet_material_property(self)->list:
        """Method that determines the material properties of the planet using weighted aggregates technique.

        Raises:
            ValueError: When the material in the input material list is not available in the material library.

        Returns:
            list: A list containing the Poisson Ratio, Young's Modulus and Yield Strength of the planet in the same order.
        """
        material_name_list = ["silicates", "iron_nickel", "methane_ice", "water_ice", "ammonia_ice"]
        try:
            for item in self.material_list:
                assert (item in material_name_list), "Unidentified Material"
        except AssertionError:
            raise ValueError
        planet_poisson = 0
        planet_youngs_modulus = 0
        planet_yield_strength = 0
        for i in range(len(self.material_list)):
            material_index = material_name_list.index(self.material_list[i])
            material = material_property_list[material_index]
            planet_poisson += self.material_fraction[i] * material["poisson_ratio"]
            planet_yield_strength += self.material_fraction[i] * material["yield_strength"]
            planet_youngs_modulus += self.material_fraction[i] * material["youngs_modulus"]
        planet_property_list = [planet_poisson, planet_youngs_modulus, planet_yield_strength]
        return planet_property_list