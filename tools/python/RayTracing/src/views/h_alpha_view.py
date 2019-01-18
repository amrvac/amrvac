"""
Module to plot H-alpha views.

Created on 7 Dec 2018

@author: Niels Claes
"""

from physics import unit_normalizations as units
from physics import opacity
from views import plotting
import print_tools
import numpy as np
import settings
import sys


def source_function(data):
    """
    Calculates the source function for given data using Eq. (7) in paper.
    @return: Value of the source function, dependent on altitude above solar surface.
             Type is double
    """
    #Altitude in cm
    H = settings.altitude * 1e5
    #Calculate dilution factor W
    W = 0.5 * ( 1 - np.sqrt(1 - (units.Rsun**2 / (units.Rsun + H)**2)) )
    
    return W * 0.17 * 4.077 * 1e-5

def get_intensity(data, axis):
    """
    Returns the integrated intensity over the given axis.
    @param data: Reduced data object.
    @param axis: (String) Line of sight (axis) to integrate over, default = x
    @return: Integrated intensity over the given axis.
             Type is np.ndarray of dimension 2
    """
    tau = opacity.get_opacity(data, axis)
    S   = source_function(data)
    
    if settings.type_simulation == "filament":
        Ibgr = 2.2 * S
        intensity = S * (1 - np.exp(-tau)) + Ibgr * np.exp(-tau)
    else:
        # Then type = prominence, default
        intensity = S * (1 - np.exp(-tau))
    
    # Calculation and integration happens correctly, but rotate
    # 2D array 90 degrees in order to be consistent with visualizations
    # through ParaView. Possibly due to Numpy - Fortran array ordening.
    # See @note in plotting.plot_surface()
    intensity = np.rot90(intensity)
        
    return intensity


def plot_h_alpha(data):
    """
    Plot the h-alpha line intensity over a given line of sight.
    @param data: Reduced data object
    """
    if data._ndim == 3:
        line_of_sight = settings.line_of_sight
        if line_of_sight == "all":
            intensity = get_intensity(data, "x")
            print("Plotting H-alpha intensity along x")
            plotting.plot_surface(data=data, matrix_2d=intensity, line_of_sight="x", title="H-alpha - x")
            intensity = get_intensity(data, "y")
            print("Plotting H-alpha intensity along y")
            plotting.plot_surface(data=data, matrix_2d=intensity, line_of_sight="y", title="H-alpha - y")
            intensity = get_intensity(data, "z")
            print("Plotting H-alpha intensity along z")
            plotting.plot_surface(data=data, matrix_2d=intensity, line_of_sight="z", title="H-alpha - z")
        elif line_of_sight == "x" or line_of_sight == "y" or line_of_sight == "z":
            intensity = get_intensity(data, line_of_sight)
            print("Plotting H-alpha intensity along %s \n" % line_of_sight)
            plotting.plot_surface(data=data, matrix_2d=intensity, line_of_sight=line_of_sight,
                                  title="H-alpha - %s" % line_of_sight)
        else:
            print("Line of sight parameter is not known. Should be 'x', 'y', 'z' or 'all'.")
            sys.exit(1)
    else:
        intensity = get_intensity(data, "x")
        print("Plotting H-alpha intensity for 2D dataset.")
        plotting.plot_surface(data=data, matrix_2d=intensity, line_of_sight="x", title="H-alpha")
    return    
    
    
    
    
        