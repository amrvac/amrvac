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

theta = np.pi / 5   # angle between x-y plane and z-axis 
phi   = np.pi / 3   # angle from x-axis towards y-axis

#Save files or not. Default = True
saveFiles = True