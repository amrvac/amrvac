"""
Created on 6 Dec 2018

Module containing the table values for the ionization degree i and the parameter f for
different temperatures and pressures at different altitudes.
Contains methods to interpolate data both in 2 and 3 dimensions.
Values taken from Table 1 in Heinzel et al. (2015), Astronomy and Astrophysics v.579, p. A16.
Table values are in cgs units.

@note: Altitude means height above the solar disk. The differences in values are due to the scattering of
       incident radiation from the solar disk, which is dependent on height.

@author: Niels Claes
"""

from scipy.interpolate import RectBivariateSpline
import numpy as np
from physics import unit_normalizations as units
import print_tools
import os
import settings

pg_table = np.asarray([0.01, 0.02, 0.05, 0.10, 0.20])           # dyn cm-2
T_table  = np.asarray([6000., 8000., 10000., 12000., 14000.])   # Kelvin
alt_table = np.asarray([10000, 20000, 30000])                   # km

# 10 000 km altitude table
ionization_10k = np.asarray([[0.74, 0.62, 0.44, 0.31, 0.20],
                             [0.83, 0.72, 0.55, 0.44, 0.35],
                             [0.87, 0.79, 0.70, 0.69, 0.73],
                             [0.91, 0.85, 0.82, 0.85, 0.89],
                             [0.93, 0.89, 0.89, 0.92, 0.94]])

f_10k = np.asarray([[5.0, 4.6, 4.2, 4.0, 4.0],
                    [6.7, 5.8, 5.0, 4.8, 4.7],
                    [8.1, 6.8, 5.3, 5.1, 5.3],
                    [9.1, 7.0, 5.0, 4.9, 5.2],
                    [9.8, 7.1, 4.9, 4.8, 5.0]])

# 20 000 km altitude table
ionization_20k = np.asarray([[0.73, 0.60, 0.41, 0.29, 0.18],
                             [0.81, 0.70, 0.52, 0.41, 0.33],
                             [0.86, 0.78, 0.68, 0.68, 0.72],
                             [0.90, 0.84, 0.81, 0.84, 0.88],
                             [0.92, 0.88, 0.89, 0.91, 0.94]])

f_20k = np.asarray([[4.7, 4.2, 3.8, 3.7, 3.6],
                    [6.3, 5.4, 4.6, 4.4, 4.3],
                    [7.6, 6.3, 4.8, 4.7, 4.9],
                    [8.6, 6.5, 4.6, 4.5, 4.8],
                    [9.2, 6.4, 4.5, 4.4, 4.6]])

# 30 000 km altitude
ionization_30k = np.asarray([[0.71, 0.58, 0.39, 0.27, 0.17],
                             [0.80, 0.68, 0.50, 0.39, 0.32],
                             [0.85, 0.76, 0.66, 0.67, 0.71],
                             [0.89, 0.83, 0.81, 0.84, 0.88],
                             [0.92, 0.88, 0.88, 0.91, 0.93]])

f_30k = np.asarray([[4.5, 4.0, 3.5, 3.4, 3.4],
                    [6.1, 5.1, 4.3, 4.0, 4.0],
                    [7.4, 6.0, 4.5, 4.3, 4.5],
                    [8.3, 6.1, 4.2, 4.2, 4.5],
                    [8.9, 6.0, 4.2, 4.1, 4.3]])

def _round_to_base(nb, base):
    """
    Method to round a given number to the nearest base
    :param nb: (int/dble) Number to round
    :param base: (int) Base to round to
    :return: Number rounded to given base as integer
    """
    return int(base * round(float(nb)/base))


def get_i_f(data, altitude=20000):
    """
    Returns the degree of ionization for each datapoint in the given input.
    :param data: Reduced data object.
    :param altitude: Altitude at which to evaluate (in km).
                     Default is 20000, otherwise rounded to nearest base 1e4 integer.
                     Type is double or integer
    :return: i | Degree of ionization at each data point of the input matrix.
                 Type is np.ndarray of dimension ndim.
             f | f-value at each data point of the input matrix,
                 multiplied by 1e16 (see table in paper).
                 Type is np.ndarray of dimension ndim.
    """
    #Perform altitude check
    if not altitude == 20000:
        altitude = _round_to_base(altitude, 10000)
        
    #Select ionization table
    if int(altitude) == 10000:
        i_table = ionization_10k
    elif int(altitude) == 20000:
        i_table = ionization_20k
    else:
        i_table = ionization_30k
        
    #Select f table
    if int(altitude) == 10000:
        f_table = f_10k
    elif int(altitude) == 20000:
        f_table = f_20k
    else:
        f_table = f_30k
        
    #Create bivariate spline approximation over rectangular mesh
    #Has to be done only once, then use it to evaluate (T, p) point
    spline_ion = RectBivariateSpline(T_table, pg_table, i_table)
    spline_f   = RectBivariateSpline(T_table, pg_table, f_table)
        
    #Re-dimensionalize temperature and pressure       
    temp = data.T * units.unit_temperature
    pg   = data.p * units.unit_pressure

    # Prevent double loading during the same run
    if data.ion is not None and data.f_param is not None:
        return data.ion, data.f_param
    else:
        #See if data is already stored to disk:
        filen = print_tools.trim_filename(settings.filename)
        if os.path.isfile("interpolated_files/ion_" + filen + ".npy") and os.path.isfile("interpolated_files/f_param_" + filen + ".npy"):
            print("Interpolated data already exists -- loading files.")
            print("    Loading Numpy files...")
            ion = np.load("interpolated_files/ion_" + filen + ".npy")
            f = np.load("interpolated_files/f_param_" + filen + ".npy")
            print("    Done.")
            data.ion = ion
            data.f_param = f*1e16
            return ion, f*1e16

    #Create matrix of same size of input
    ion = np.zeros_like(temp)
    f   = np.zeros_like(temp)
    
    #Fast iteration over array elements (calls the C array operator API)
    it = np.nditer(temp, flags=['multi_index'])
    
    print("Interpolating matrix for ionization and f.")
    
    if data._ndim == 2:
        tot_points = len(data.T[0, :])
    else:
        tot_points = len(data.T[0, 0, :])

    ctr = 0
    while not it.finished:
        #Get current index of iterator, Type = Tuple
        idx = it.multi_index          
        
        #Get temperature and pressure at current index
        t_i = temp[idx]
        p_i = pg[idx]
        
        #Interpolate ionization degree, save to current index
        ion[idx] = spline_ion.ev(t_i, p_i)
        #Interpolate f, save to current index
        f[idx] = spline_f.ev(t_i, p_i)
        
        #Advance iterator
        ctr += 1
        #Print out progress
        if ctr % 250 == 0:
            print_tools.progress(idx[-1], tot_points, '-- interpolating...')
        it.iternext()
    print_tools.progress(tot_points, tot_points, '-- completed.')
    print("\n")

    # save Numpy arrays for easy acces later on
    if settings.saveFiles:
        np.save("interpolated_files/ion_" + filen, ion)
        np.save("interpolated_files/f_param_" + filen, f)
        print("Interpolated arrays saved to")
        print("    interpolated_files/ion_" + filen + ".npy")
        print("    interpolated_files/f_param_" + filen + ".npy")

    data.ion = ion
    data.f_param = f * 1e16
        
    # Parameter f is tabulated in units of 10^16 cm-3
    return ion, f*1e16
