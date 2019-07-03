"""
Module with the main plotting function.

Created on 09 Jan 2018

@author: Niels Claes
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors
import settings
import numpy as np
import print_tools
from physics import unit_normalizations as units

# Outside of method definitions

mpl.rcParams['xtick.labelsize'] = settings.axis_fontsize
mpl.rcParams['ytick.labelsize'] = settings.axis_fontsize


def plot_surface(data, matrix_2d, line_of_sight, title=None, xlabel=None, ylabel=None, cblabel=None, cmap=None,
                 logscale=True, cbformat=None):
    """
    Creates a surface plot of the integrated data.
    :param data: Reduced data object.
    :param matrix_2d: Integrated matrix along the line of sight.
                      Type is np.ndarray of dimension 2
    :param line_of_sight: (String) Line of sight corresponding to matrix_2d.
                          Can be 'x', 'y' or 'z', default = 'x'.
    :param title: (Optional, String) Title for the figure.
    :param xlabel: (Optional, String) Label for horizontal axis.
    :param ylabel: (Optional, String) Label for vertical axis.
    :param cblabel: (Optional, String) Label for the color bar.
    :param cmap: (Optional, String) Color map to use in plotting. Default is 'Reds'.
    :param logscale: (Optional, boolean) Plot on logscale or not. Default is True.
    :param cbformat: (Optional, String) Format string for colorbar ticks.
    """
    hori_ax, vert_ax = get_arrays(data, line_of_sight)
    bounds = [np.min(hori_ax), np.max(hori_ax), np.min(vert_ax), np.max(vert_ax)]

    if cmap is not None:
        cmap = cmap
    else:
        cmap = "Reds"

    # Create figure
    fig, ax = plt.subplots(1)
    if logscale:
        im = plt.imshow(matrix_2d, norm=mpl.colors.LogNorm(), cmap=cmap, extent=bounds)
    else:
        im = plt.imshow(matrix_2d, cmap=cmap, extent=bounds)
    im.set_interpolation("bilinear")

    # Colorbar formatting
    if cbformat is not None and not logscale:
        colorbar = fig.colorbar(im, format=cbformat)
    else:
        colorbar = fig.colorbar(im)


    #Labels
    if cblabel is not None:
        colorbar.set_label(cblabel, fontsize=settings.cbar_fontsize)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=settings.axis_fontsize)
    else:
        ax.set_xlabel("[%.0e cm]" % units.unit_length, fontsize=settings.axis_fontsize)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=settings.axis_fontsize)
    else:
        ax.set_ylabel("[%.0e cm]" % units.unit_length, fontsize=settings.axis_fontsize)

    if title is not None:
        ax.set_title(title, fontsize=settings.title_fontsize)

    #Bounds
    ax.set_xlim([bounds[0], bounds[1]])
    ax.set_ylim([bounds[2], bounds[3]])

    #Magnetic field lines
    if settings.plot_Blines:
        impose_fieldlines(data, line_of_sight, ax)

    fig.tight_layout()

    return


def impose_fieldlines(data, line_of_sight, ax):
    """
    Imposes magnetic field lines on surface plot, requires setting the angles
    theta and beta in unit_normalizations. Theta is the angle between the x-axis and vertical y-axis,
    while beta is the angle from the x-axis towards the y-axis.
    :param data: Reduced data object.
    :param line_of_sight: Line of sight
    :param ax: Matplotlib axis object
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
            slope = np.tan(angle)
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
    :param data: Reduced data object.
    :param axis: (String) Line of sight, default = x
    :return x1, x2: Horizontal and vertical axis, respectively.
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