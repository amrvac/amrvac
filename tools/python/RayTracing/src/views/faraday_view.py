"""
Module to plot the Faraday rotation effect.

Created on 09 Jan 2019

@author: Niels Claes
"""

from physics import unit_normalizations as units
from physics import ionization
from maths import integration
from views import plotting
import print_tools
import settings
import numpy as np
import sys


def get_faraday_strength(data, axis):
    """
    Method to calculate the overall strength of the Faraday rotation effect, i.e.
    the rotation measure (RM) times the light wavelength squared.
    For 3D, the RM is integrated along the line of sight, for 2D it is simply interpolated.
    :param data: Reduced data object.
    :param axis: (String) Line of sight
    :return: - Strength of Faraday effect interpolated from input if ndim < 3
             - Strength of Faraday effect over the line of sight if ndim = 3
             Type is np.ndarray of dimension 2
    """
    prefactor = units.ec**3 / (2 * np.pi * units.m_e**2 * units.c**4)

    # Re-dimensionalize data
    temp = data.T * units.unit_temperature
    pg   = data.p * units.unit_pressure
    if axis == 'x':
        b_para = data.b1 * units.unit_magneticfield
    elif axis == 'y':
        b_para = data.b2 * units.unit_magneticfield
    else:
        b_para = data.b3 * units.unit_magneticfield

    ion, f = ionization.get_i_f(data, settings.altitude)
    # Get number density of electrons using (10) from paper
    n_e = pg / ((1 + 1.1 / ion) * units.k_B * temp)  # cm-3

    integral_argument = n_e * b_para

    if not data._ndim == 3:
        rm_measure = prefactor * integral_argument
    else:
        rm_measure = prefactor * integration.integrate_los(data, integral_argument, axis)

    # Same argument as in h_alpha_view, rotate 90 degrees
    faraday = np.rot90(rm_measure)

    # NOTE: rotation measure (rm_measure) is in rad/cm2, but faraday itself is in radians!
    if settings.include_lambda:
        faraday = faraday * settings.faraday_lambda**2
    else:
        # Convert from rad/cm2 to rad/m2
        faraday = faraday * 1e4

    return faraday


def plot_faraday(data):
    """
    Plot the faraday rotation strength along a given line of sight.
    :param data: Reduced data object.
    """
    if settings.include_lambda:
        cbtitle=r"Faraday Rotation $\beta$ [rad]" +"\n$\lambda$ = {} cm".format(settings.faraday_lambda)
    else:
        cbtitle=r"Faraday RM [rad m$^{-2}$]"
    if data._ndim == 3:
        # Check on line of sight
        if settings.line_of_sight == "all":
            line_of_sight = ["x", "y", "z"]
        elif settings.line_of_sight == "x" or settings.line_of_sight == "y" or settings.line_of_sight == "z":
            line_of_sight = [settings.line_of_sight]
        else:
            sys.exit("Line of sight parameter is not known. Should be 'x', 'y', 'z' or 'all'."
                     "\nCurrent value: '{}'".format(settings.line_of_sight))

        for los in line_of_sight:
            faraday = get_faraday_strength(data, los)
            print("Plotting Faraday effect along {}".format(los))
            title = "Faraday - {}\n{}".format(los, print_tools.trim_filename(settings.filename))
            plotting.plot_surface(data=data, matrix_2d=faraday, line_of_sight=los, title=title,
                                  cmap=settings.cmap_faraday, logscale=settings.logscale_faraday, cblabel=cbtitle,
                                  cbformat=settings.cbar_format_faraday)

    else:
        faraday = get_faraday_strength(data, "x")
        print("Plotting Faraday effect for 2D dataset.")
        plotting.plot_surface(data=data, matrix_2d=faraday, line_of_sight="x",
                              title="Faraday\n{}".format(print_tools.trim_filename(settings.filename)),
                              cmap=settings.cmap_faraday, logscale=settings.logscale_faraday, cblabel=cbtitle,
                              cbformat=settings.cbar_format_faraday)
    return

