import numpy as np

class constants():
    """
    Contains physical constants, both in SI and cgs units.
    """
    def __init__(self, cgs):
        if cgs:
            self.R = 83144598.0             # erg / (K mol)
            self.m_p = 1.6726219e-24        # g
            self.m_e = 9.10938356e-28       # g
            self.k_B = 1.38064852e-16       # erg / K
            self.ec = 4.80320467299766e-10  # statcoul
            self.c = 2.99792458e+10         # cm / s
            self.Rsun = 6.957e+10           # cm
        else:
            self.R = 8.3144598              # J / (K mol)
            self.m_p = 1.672621898e-27      # kg
            self.m_e = 9.10938356e-31       # kg
            self.k_B = 1.38064852e-23       # J / K
            self.ec = 1.6021766208e-19      # C
            self.c = 299792458.0            # m / s
            self.Rsun = 695700000           # m


class units(constants):
    """
    Sets the unit normalisations for a given dataset. Initially at default values when loading in the dataset.
    Correct physical values can be set by calling set_units.
    """
    def __init__(self, header):
        self.cgs = header.get("cgs", True)
        self.unit_length = header.get("unit_length", 1.0)
        self.unit_numberdensity = header.get("unit_numberdensity", 1.0)
        self.unit_temperature = header.get("unit_temperature", 1.0)
        self.unit_velocity = header.get("unit_velocity", 0)
        self.gamma = header.get("gamma", 5.0/3.0)
        super(units, self).__init__(self.cgs)

        self.mu = 0
        self.unit_pressure = None
        self.unit_time = None
        self.unit_magneticfield = None
        self.unit_luminosity = None
        self.unit_dlambdadT = None
        self.unit_conduction = None

        self._calculate_units()

    def set_units(self, unit_length, unit_numberdensity, unit_temperature=1.0, unit_velocity=0):
        self.unit_length = unit_length
        self.unit_numberdensity = unit_numberdensity
        self.unit_temperature = unit_temperature
        self.unit_velocity = unit_velocity
        self._calculate_units()

    def _calculate_units(self):
        He_abundance = 0.1
        H_abundance = 1 - He_abundance
        self.mu = 1. / (2 * H_abundance + 0.75 * He_abundance)

        self.unit_density = (1.0 + 4.0 * He_abundance) * self.m_p * self.unit_numberdensity
        if self.unit_velocity == 0:  # <-- unit temperature is defined
            self.unit_pressure = (2.0 + 3.0 * He_abundance) * self.unit_numberdensity * self.k_B * self.unit_temperature
            self.unit_velocity = np.sqrt(self.unit_pressure / self.unit_density)
        else:  # <-- unit velocity is defined
            self.unit_pressure = self.unit_density * self.unit_velocity ** 2
            self.unit_temperature = self.unit_pressure / ((2.0+3.0*He_abundance) * self.unit_numberdensity * self.k_B)
        self.unit_magneticfield = np.sqrt(4 * np.pi * self.unit_pressure)
        self.unit_time = self.unit_length / self.unit_velocity

        # normalisation for the cooling function
        self.unit_luminosity = self.unit_pressure / (self.unit_numberdensity**2 * self.unit_time*(1.0+2.0*He_abundance))
        # when calculating dL/dT, unit_luminosity must be divided by unit_temperature
        self.unit_dldt = self.unit_luminosity / self.unit_temperature

        # normalisation for thermal conduction coefficients
        self.unit_conduction = self.unit_density * self.unit_length * self.unit_velocity ** 3 / self.unit_temperature

    def check_default_units(self):
        if (self.unit_length == 1.0 and self.unit_numberdensity == 1.0 and
                (self.unit_temperature == 1.0 or self.unit_velocity == 0)):
            print("[WARNING] Unit normalisations are still at their default settings. For a consistent calculation "
                  "these must be manually set to their respective values.\n"
                  "          This can be done through "
                  "ds.set_units(unit_lengh=..., unit_numberdensity=..., unit_temperature=...)")
