"""
Created on 6 Dec 2018

Module to calculate the absorption coefficient for given physical parameters.

@author: Niels Claes
"""
from physics import ionization
from physics import unit_normalizations as units
import numpy as np
import settings

def gaussian(temp):
    """
    Calculates absorption line profile using Gaussian broadening.
    Equations used are (5) and (6) in paper.
    @param temp: Temperature, in Kelvin
              Type is np.ndarray of dimension ndim.
    @return: Gaussian function for H-alpha absorption line profile.
             Type is np.ndarray of dimension ndim, units = seconds
    """
    #Microturbulence
    ksi = 5*1e5 #cm/s
    
    nu_0 = units.c / (6562.8 * 1e-7)        #s-1
    
    delta_nu  = 0
    
    delta_nuD = (nu_0 / units.c) * np.sqrt(2 * units.k_B * temp / units.m_p  +  ksi**2)
    
    phi_nu = 1.0 / (np.sqrt(np.pi) * delta_nuD) * np.exp(-delta_nu / delta_nuD)**2
    
    return phi_nu
  

def get_kappa(data):
    """
    Calculates the absorption coefficient for the H-alpha line using
    equation (4) from paper.
    @param data: Reduced data object.
    @return: Absorption coefficient for H-alpha line (in cgs)
             Type is np.ndarray of dimension ndim.
    """
    #Dimensionalize relevant variables
    temp = data.T * units.unit_temperature  # K
    pg   = data.p * units.unit_pressure     # dyn cm-2
    
    #Retrieve ionization degree and parameter f from data
    ion, f = ionization.get_i_f(data, settings.altitude)
    
    #Get number density of electrons using (10) from paper
    n_e = pg / ( (1 + 1.1/ion) * units.k_B * temp )     #cm-3
    
    #Calculate n2 from ne^2 = f*n2
    n2  = n_e**2 / f    #cm-3
    
    f23 = 0.6407    #Oscillator strength H-alpha (dimensionless)
    
    #Calculate absorption coefficient using (4) from paper
    kappa = (np.pi * units.ec**2 / (units.m_e * units.c) ) * f23 * n2 * gaussian(temp)
    
    return kappa



