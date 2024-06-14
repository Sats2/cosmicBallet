import Constants as const
from typing import Union


def mass_conversion(mass_list:list, unit:str)->list:
    """Function that converts the mass to SI Units

        Converts the mass values in the mass_list to SI Units. All elements within the mass_list need to be of a single unit. The
        units are provided to the function via shorthand notation.
            Accepted Units are: lbs (pounds), Me (Earth Mass) and Ms (Solar Mass)

    Args:
        mass_list (list): Mass values needed for conversion. Single values need to be input as a list
        unit (str): The current unit of the mass value

    Raises:
        TypeError: Raised if the mass_list or unit are not of specified data type.
        ValueError: Raised if the unit is unrecognized

    Returns:
        list: Mass list in SI Units.

    Example:
        test_mass = [1,2]
        test_mass_si = mass_conversion(mass_list=test_mass, unit="Ms")
        print(test_mass_si)
        >> [1.9884e30, 3.9768e30]
    """
    try:
        assert isinstance(mass_list, list), "Arguement mass_list must be of type list"
        assert isinstance(unit, str), "Arguement unit must be of type list"
    except AssertionError:
        raise TypeError
    if unit == "lbs":
        conversion_factor = const.lbs_to_kg
    elif unit == "Me":
        conversion_factor = const.Me
    elif unit == "Ms":
        conversion_factor = const.Ms
    else:
        raise ValueError("Unit not recognized. Select from 'lbs', 'Me' or 'Ms'.")
    for item in mass_list:
        item *= conversion_factor
    return mass_list


def velocity_conversion(velocity_list:list, unit:str)->list:
    """Function that converts the velocity to SI Units

        Converts the velocity values in the velocity_list to SI Units. All elements within the velocity_list need to be of a single 
        unit. The units are provided to the function via shorthand notation.
            Accepted Units are: ftps (feet/second), mph (miles/hour) and kmph (kilometer/hour)

    Args:
        velocity_list (list): Velocity values needed for conversion. Single values need to be input as a list
        unit (str): The current unit of the velocity value

    Raises:
        TypeError: Raised if the velocity_list or unit are not of specified data type.
        ValueError: Raised if the unit is unrecognized

    Returns:
        list: Velocity list in SI Units.

    Example:
        test_velocity = [100,200]
        test_velocity_si = velocity_conversion(velocity_list=test_velocity, unit="Ms")
        print(test_velocity_si)
        >> [44.704, 89.408]
    """
    try:
        assert isinstance(velocity_list, list), "Arguement speed_list must be of type list"
        assert isinstance(unit, str), "Arguement unit must be of type list"
    except AssertionError:
        raise TypeError
    if unit == "ftps":
        conversion_factor = const.ftps_to_mtps
    elif unit == "mph":
        conversion_factor = const.mph_to_mtps
    elif unit == "kmph":
        conversion_factor = const.kmph_to_mtps
    else:
        raise ValueError("Unit not recognized. Select from 'ftps', 'mph', or 'kmph'.")
    for item in velocity_list:
        item *= conversion_factor
    return velocity_list


def distance_conversion(distance_list:list, unit:str)->list:
    """Function that converts the distance to SI Units

        Converts the distance values in the distance_list to SI Units. All elements within the distance_list need to be of a 
        single unit. The units are provided to the function via shorthand notation.
            Accepted Units are: AU (astronomical unit), ly (light year), mi (mile), km (kilometer)

    Args:
        distance_list (list): Distance values needed for conversion. Single values need to be input as a list
        unit (str): The current unit of the distance value

    Raises:
        TypeError: Raised if the distance_list or unit are not of specified data type.
        ValueError: Raised if the unit is unrecognized

    Returns:
        list: Distance list in SI Units.

    Example:
        test_distance = [1,2]
        test_distance_si = distance_conversion(distance_list=test_distance, unit="ly")
        print(test_distance_si)
        >> [9.4607e15, 1.89214e16]
    """
    try:
        assert isinstance(distance_list, list), "Arguement distance_list must be of type list"
        assert isinstance(unit, str), "Arguement unit must be of type list"
    except AssertionError:
        raise TypeError
    if unit == "AU":
        conversion_factor = const.AU
    elif unit == "ly":
        conversion_factor = const.ly
    elif unit == "mi":
        conversion_factor = const.mile_to_mt
    elif unit == "km":
        conversion_factor = const.km_to_mt
    else:
        raise ValueError("Unit not recognized. Select from 'AU', 'ly', 'mi' or 'km'.")
    for item in distance_list:
        item *= conversion_factor
    return distance_list


def temperature_conversion(temperature:Union[float,int], unit:str)->float:
    """Function that converts the temperature to SI Units

        Converts the temperature to SI Units (kelvin). The accepted units for conversion and method of conversion are:
            C (celcius) -> temp_kelvin = temp_celcius + 273.15
            F (Fahrenheit) -> temp_kelvin = ((temp_fahrenheit - 32)/1.8) + 273.15
            R (Rankine) -> temp_kelvin = temp_rankine * 5 / 9

    Args:
        temperature (float/int): Temperature value for conversion to SI Units
        unit (str): Unit of the temperature value

    Raises:
        TypeError: Raised if the temperature or unit are of unspecified datatypes
        ValueError: Raised if the unit is unspecified

    Returns:
        float: Temperature value in SI Units

    Example:
        test_temp = 212
        test_temp_si = temperature_conversion(temperature=test_temp, unit="F")
        print(test_temp_si)
        >> 373.15
    """
    try:
        assert isinstance(temperature, (float,int)), "Arguement temperature must be of type float/int"
        assert isinstance(unit, str), "Arguement unit must of type string"
    except AssertionError:
        raise TypeError
    if unit == "C":
        temperature += 273.15
    elif unit == "F":
        temperature = ((temperature - 32) / 1.8) + 273.15
    elif unit == "R":
        temperature *= (5/9)
    else:
        raise ValueError("Unit not recognized. Select from 'C', 'F', or 'R'.")
    return temperature