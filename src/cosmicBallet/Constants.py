# All universal constants are in SI Units unless explicitly stated

G = 6.67430e-11             # Universal Gravitational Constant
MS = 1.9884e30              # Solar Mass (Mass of the Sun)
ME = 5.972e24               # Earth Mass (Mass of the Earth)
AU = 1.49598e11             # Astronomical Unit (Distance between the Earth and the Sun)
LY = 9.4607e15              # Light Year (Distance Travelled by Light in 1 year)
C = 2.99792458e8            # Speed of Light


# Conversion Factors

LBS_TO_KG = 0.4535924       # Conversion from pounds (lbs) to kilograms (kg)
FTPS_TO_MTPS = 0.3048       # Conversion from feet/second (ft/s) to meter/second (m/s)
MPH_TP_MTPS = 0.44704       # Conversion from miles/second (mile/s) to meter/second (m/s)
KMPH_TO_MTPS = 0.2777778    # Conversion from kilometers/hour (km/hr) to meter/second (m/s)
MILE_TO_MT = 1609.344       # Conversion from miles (mi) to meters (m)
KM_TO_MT = 1000.0           # Conversion from kilometer (km) to meters (m) 
HOUR_TO_SEC = 3600          # Conversion from hours to seconds (s)
DAY_TO_SEC = 24*HOUR_TO_SEC # Conversion from days to seconds (s)
MON_TO_SEC = 30*DAY_TO_SEC  # Conversion from months to seconds (s)
YEAR_TO_sEC = 12*MON_TO_SEC # Conversion from years to seconds (s)


# Materials Dictionaries with material properties, all values in pascals (Pa)

silicates = {
            "name":"silicates",
            "poisson_ratio":0.25,
            "youngs_modulus":70e9,
            "yield_strength":150e6
            }
iron_nickel = {
            "name":"iron_nickel",
            "poisson_ratio":0.3,
            "youngs_modulus":150e9,
            "yield_strength":300e6
            }
water_ice = {
            "name":"water_ice",
            "poisson_ratio":0.33,
            "youngs_modulus":9.3e9,
            "yield_strength":5e6
            }
ammonia_ice = {
            "name":"ammonia_ice",
            "poisson_ratio":0.33,
            "youngs_modulus":7e9,
            "yield_strength":5e6
            }
methane_ice = {
            "name":"methane_ice",
            "poisson_ratio":0.33,
            "youngs_modulus":4.5e9,
            "yield_strength":5e6
}
material_property_list = [silicates, iron_nickel, water_ice, ammonia_ice, methane_ice]