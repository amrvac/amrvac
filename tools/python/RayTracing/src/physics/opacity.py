"""
Created on 7 Dec 2018

Module to calculate the opacity for the H-alpha views. This opacity is then integrated along
the line of sight. Works in 1d, 2d and 3d, along any of the given coordinate axes.
@note: At the moment, only integration along one of the coordinate axes is supported.

@author: Niels Claes
"""

from maths import integration
from physics import ionization
from physics import unit_normalizations as units
import settings
import numpy as np


def gaussian(temp):
    """
    Calculates absorption line profile using Gaussian broadening.
    Equations used are (5) and (6) in paper.
    @param temp: Temperature, in Kelvin
              Type is np.ndarray of dimension ndim.
    @return: Gaussian function for H-alpha absorption line profile.
             Type is np.ndarray of dimension ndim, units = seconds
    """
    # Microturbulence
    ksi = 5 * 1e5  # cm/s

    # H-alpha wavelength is 6562.8 Angstrom
    nu_0 = units.c / (6562.8 * 1e-8)  # s-1

    delta_nu = 0

    delta_nuD = (nu_0 / units.c) * np.sqrt(2 * units.k_B * temp / units.m_p + ksi ** 2)

    phi_nu = 1.0 / (np.sqrt(np.pi) * delta_nuD) * np.exp(-delta_nu / delta_nuD) ** 2

    return phi_nu


def get_kappa(data):
    """
    Calculates the absorption coefficient for the H-alpha line using
    equation (4) from paper.
    @param data: Reduced data object.
    @return: Absorption coefficient for H-alpha line (in cgs)
             Type is np.ndarray of dimension ndim.
    """
    # Dimensionalize relevant variables
    temp = data.T * units.unit_temperature  # K
    pg   = data.p * units.unit_pressure  # dyn cm-2

    # Retrieve ionization degree and parameter f from data
    ion, f = ionization.get_i_f(data, settings.altitude)

    # Get number density of electrons using (10) from paper
    n_e = pg / ((1 + 1.1 / ion) * units.k_B * temp)  # cm-3

    # Calculate n2 from ne^2 = f*n2
    n2 = n_e ** 2 / f  # cm-3

    f23 = 0.6407  # Oscillator strength H-alpha (dimensionless)

    # Calculate absorption coefficient using (4) from paper
    kappa = (np.pi * units.ec ** 2 / (units.m_e * units.c)) * f23 * n2 * gaussian(temp)

    return kappa
    

def get_opacity(data, axis):
    """
    Returns the opacity along the line of sight.
    @param data: Reduced data object
    @param axis: (String) Line of sight, default = x
    @return: - Opacity interpolated from input if ndim < 3
             - Integrated opacity along line of sight if ndim == 3
             Type is np.ndarray of dimension 2
    """
    #Retrieve absorption coefficient
    kappa = get_kappa(data)
    
    #If 1D or 2D, simply return the grid itself
    if not data._ndim == 3:
        return kappa
    else:
        return integration.integrate_los(data, kappa, axis)
        
        
