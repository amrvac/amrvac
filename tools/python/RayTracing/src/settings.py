"""
Module containing the setup parameters for the program. The unit normalisations here must be changed according to
the desired values, to ensure a correct re-dimensionalization of the units.
The filename of the .dat file can be set here or passed as an argument from the command line.
In the latter case the value defined here is overruled.

@author: Niels Claes
"""

import numpy as np


filename = "dat_files/ot_2d0002.dat"

# User defined normalizations
unit_length             = 1e9    # cm
unit_temperature        = 1e6    # K
unit_numberdensity      = 1e9    # cm-3
altitude                = 20000  # km

#Can be filament or prominence
type_simulation         = "prominence"

# angle between x-y plane and z-axis
theta = np.pi / 5
# angle from x-axis towards y-axis
phi   = np.pi / 3

# Line of sight.
line_of_sight = "z"
# Plot on logscale or not.
logscale = True
# Color map to use when plotting.
cmap = "Reds"

#Calculate H-alpha view.
halpha = True

#Calculate Faraday effect.
faraday = True
faraday_lambda = 460  # nm (blue visible)

#Save files or not. Default = True
saveFiles = True

#Use multiprocessing to interpolate blocks.
multiple_procs = False
# Number of processors
nb_of_procs    = 4

#Whether or not to plot magnetic field lines.
plot_Blines = True