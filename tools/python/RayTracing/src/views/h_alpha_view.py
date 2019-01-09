"""
Created on 7 Dec 2018

@author: Niels Claes
"""

import matplotlib.pyplot as plt
from physics import unit_normalizations as units
import numpy as np
import matplotlib as mpl
import sys
from physics import opacity
from matplotlib.ticker import LogFormatter 
import settings

def plot_surface(data, intensity, hori_ax, vert_ax, logscale, cmap, plotB):
    """
    Creates a surface plot of the integrated data.
    @param data: Reduced data object.
    @param intensity: Integrated intensity along the line of sight.
                      Type is np.ndarray of dimension 2
    @param hori_ax: Horizontal axis, dependent on line of sight.
                    Type is np.ndarray of dimension 1
    @param vert_ax: Vertical axis, dependent on line of sight.
                    Type is np.ndarray of dimension 1
    @param logscale: Plot data on log scale or not. Default = True
    @param cmap: Colormap for surface plot. Default = Reds
    @param plotB: Whether or not to impose magnetic field lines.
                  Angles with axes must be given in unit_normalizations. Default = False
    @note: - X integration: visual is from outside axis to origin, z is upwards, y is to the right
           - Y integration: visual is from origin to outside axis, z is upwards, x is to the right
           - Z integration: visual is from outside axis to origin, y is upwards, x is to the right
    """
    print("Plotting H-alpha intensity along %s \n" % line_of_sight)
    bounds = [np.min(hori_ax), np.max(hori_ax), np.min(vert_ax), np.max(vert_ax)]
    
    #Create figure
    fig, ax = plt.subplots(1)
    if logscale:
        im = plt.imshow(intensity, norm=mpl.colors.LogNorm(), cmap=cmap, extent=bounds)
    else:
        im = plt.imshow(intensity, cmap=cmap, extent=bounds)
    im.set_interpolation("bilinear")
    
    colorb = fig.colorbar(im)  
    colorb.set_label(r"$I(\nu)$ [erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$ Hz$^{-1}$]")
    ax.set_xlabel("[%.0e cm]" % units.unit_length)
    ax.set_ylabel("[%.0e cm]" % units.unit_length)
    title = "Integrated H-alpha line intensity along %s-axis \n %s" % (line_of_sight, settings.filename)
    ax.set_title(title)
    ax.set_xlim([bounds[0], bounds[1]])
    ax.set_ylim([bounds[2], bounds[3]])
    
    if plotB:
        impose_fieldlines(data, fig, ax)

    return

def impose_fieldlines(data, fig, ax):
    """
    Imposes magnetic field lines on surface plot, requires setting the angles
    theta and beta in unit_normalizations. Theta is the angle between the x-axis and vertical y-axis,
    while beta is the angle from the x-axis towards the y-axis.
    @param data: Reduced data object.
    @param fig: Matplotlib figure object
    @param ax: Matplotlib axis object
    """
    theta = settings.theta
    phi   = settings.phi
    
    if data._ndim == 2:
        hori = data._x
        vert = data._y
        slope = np.tan(theta)
    else:
        if line_of_sight == "x":
            hori = data._y
            vert = data._z
            # Projection of position vector on y-z plane, angle with horizontal y-axis
            angle = np.arctan(np.sin(theta) / (np.cos(theta) * np.sin(phi)))
            # NOTE: according to comparison with paraview, this should be pi/2 - angle??
            slope = np.tan(np.pi/2 - angle)
        elif line_of_sight == "y":
            hori = data._x
            vert = data._z
            # Projection of position vector on x-z plane, angle with horizontal x-axis
            angle = np.arctan(np.sin(theta) / (np.cos(theta) * np.cos(phi)))
            slope = np.tan(angle)
        else:
            hori = data._y
            vert = data._x
            # Projection of position vector on x-y plane, angle with horizontal x-axis
            slope = np.tan(phi)

    for c in np.arange(-2*np.max(vert), np.max(vert), 0.4):
        y = slope * hori + c
        ax.plot(hori, y, color="white", lw=2)
    return

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
    # through ParaView. See @note in plot_surface()
    intensity = np.rot90(intensity)
        
    return intensity

def get_arrays(data, axis):
    """
    Returns the arrays to plot over, i.e. the ones not equal to the line of sight.
    @param data: Reduced data object.
    @param axis: (String) Line of sight, default = x
    @return x1, x2: Horizontal and vertical axis, respectively.
                    Types are np.ndarrays of dimension 1
    """
    if data._ndim == 2:
        x1 = data._x
        x2 = data._y
    else:
        if axis == "x":
            x1 = data._y
            x2 = data._z
        elif axis == "y":
            x1 = data._x
            x2 = data._z
        else:
            x1 = data._x
            x2 = data._y
    return x1, x2

def plot_h_alpha(data, axis="x", cmap="Reds", logscale=True, plotB=False):
    """
    Plot the h-alpha line intensity over a given line of sight.
    @param data: Reduced data object
    @param axis: (String) Line of sight, default = x
    @param cmap: (String) Color map to plot, default = Reds
    @param logscale: Plot intensity on logscale or not. Default = True
    @param plotB: Draw magnetic field lines. Default = False
    """
    #Define global variable for access across methods in this module.
    global line_of_sight
    line_of_sight = axis
    
    x1, x2 = get_arrays(data, axis)
    
    intensity = get_intensity(data, axis)
    plot_surface(data, intensity, hori_ax=x1, vert_ax=x2, cmap=cmap, logscale=logscale, plotB=plotB)
    
    return    
    
    
    
    
        