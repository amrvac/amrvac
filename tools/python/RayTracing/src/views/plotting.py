"""
Module with the main plotting function.

Created on 09 Jan 2018

@author: Niels Claes
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import settings
import numpy as np
import print_tools
from physics import unit_normalizations as units

def plot_surface(data, matrix_2d, line_of_sight, title=None, xlabel=None, ylabel=None, cblabel=None):
    """
    Creates a surface plot of the integrated data.
    @param data: Reduced data object.
    @param matrix_2d: Integrated matrix along the line of sight.
                      Type is np.ndarray of dimension 2
    @param line_of_sight: (String) Line of sight corresponding to matrix_2d.
                          Can be 'x', 'y' or 'z', default = 'x'.
    @param title: (Optional, String) Title for the figure.
    @param xlabel: (Optional, String) Label for horizontal axis.
    @param ylabel: (Optional, String) Label for vertical axis.
    @param cblabel: (Optional, String) Label for the color bar.
    @note: - X integration: visual is from outside axis to origin, z is upwards, y is to the right
           - Y integration: visual is from origin to outside axis, z is upwards, x is to the right
           - Z integration: visual is from outside axis to origin, y is upwards, x is to the right
    """
    hori_ax, vert_ax = get_arrays(data, line_of_sight)
    bounds = [np.min(hori_ax), np.max(hori_ax), np.min(vert_ax), np.max(vert_ax)]

    # Create figure
    fig, ax = plt.subplots(1)
    if settings.logscale:
        im = plt.imshow(matrix_2d, norm=mpl.colors.LogNorm(), cmap=settings.cmap, extent=bounds)
    else:
        im = plt.imshow(matrix_2d, cmap=settings.cmap, extent=bounds)
    im.set_interpolation("bilinear")
    colorb = fig.colorbar(im)

    #Labels
    if cblabel is not None:
        colorb.set_label(cblabel)
    else:
        colorb.set_label("Intensity [cgs]")
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    else:
        ax.set_xlabel("[%.0e cm]" % units.unit_length)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel("[%.0e cm]" % units.unit_length)

    filename = print_tools.trim_filename(settings.filename)
    if title is not None:
        ax.set_title(title + "\n" + filename)
    else:
        ax.set_title(filename)

    #Bounds
    ax.set_xlim([bounds[0], bounds[1]])
    ax.set_ylim([bounds[2], bounds[3]])

    #Magnetic field lines
    if settings.plot_Blines:
        impose_fieldlines(data, line_of_sight, ax)

    return


def impose_fieldlines(data, line_of_sight, ax):
    """
    Imposes magnetic field lines on surface plot, requires setting the angles
    theta and beta in unit_normalizations. Theta is the angle between the x-axis and vertical y-axis,
    while beta is the angle from the x-axis towards the y-axis.
    @param data: Reduced data object.
    @param line_of_sight: Line of sight
    @param ax: Matplotlib axis object
    """
    theta = settings.theta
    phi = settings.phi

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
            # @note: according to comparison with paraview, this should be pi/2 - angle
            slope = np.tan(np.pi / 2 - angle)
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

    for c in np.arange(-2 * np.max(vert), np.max(vert), 0.4):
        y = slope * hori + c
        ax.plot(hori, y, color="white", lw=2)
    return

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