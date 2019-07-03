"""
Created on 6 Dec 2018

Module containing the unit normalizations from MPI-AMRVAC.

@author: Niels Claes
"""

import numpy as np
import settings

#  =========================================
try:
    import astropy.constants as ctes
    R    = ctes.R.cgs.value
    m_p  = ctes.m_p.cgs.value
    m_e  = ctes.m_e.cgs.value
    k_B  = ctes.k_B.cgs.value
    ec   = ctes.e.gauss.value
    c    = ctes.c.cgs.value
    Rsun = ctes.R_sun.cgs.value
except ImportError:
    # If astropy import fails (eg. not installed),
    # hardcode physical parameters
    R    = 83144598.0            # erg / K 
    m_p  = 1.6726219e-24         # g
    m_e  = 9.10938356e-28        # g
    k_B  = 1.38064852e-16        # erg / K
    ec   = 4.80320467299766e-10  # statcoul
    c    = 2.99792458e+10        # cm / s
    Rsun = 6.957e+10             # cm

He_abundance = 0.1
H_abundance = 1 - He_abundance
mu = 1. / (2 * H_abundance + 0.75 * He_abundance)

unit_length        = settings.unit_length
unit_temperature   = settings.unit_temperature
unit_numberdensity = settings.unit_numberdensity

# Other normalisations
unit_density = (1.0 + 4.0 * He_abundance) * m_p * settings.unit_numberdensity
unit_pressure = (2.0 + 3.0 * He_abundance) * unit_numberdensity * k_B * unit_temperature
unit_velocity = np.sqrt(unit_pressure / unit_density)
unit_magneticfield = np.sqrt(4 * np.pi * unit_pressure)
unit_time = unit_length / unit_velocity





