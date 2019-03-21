"""
Module containing the setup parameters for the program. The unit normalisations here must be changed according to
the desired values, to ensure a correct re-dimensionalization of the units.
The filename of the .dat file can be set here or passed as an argument from the command line.
In the latter case the value defined here is overruled.

@author: Niels Claes
"""

import numpy as np


filename = "PATH_TO_DATFILE"


# ========== NORMALIZATION SETTINGS ==========
# User defined normalizations
unit_length             = 1e9    # cm
unit_temperature        = 1e6    # K
unit_numberdensity      = 1e9    # cm-3


# ========== USER PARAMETERS ==========
# Altitude above the solar disk
altitude                = 20000  # km
#Type of simulation, 'filament' or 'prominence'
type_simulation         = "prominence"
# These are only needed when plot_Blines = True
# angle between x-y plane and z-axis
theta = np.pi / 5
# angle from x-axis towards y-axis
phi   = np.pi / 3


# ========== ANIMATION SETTINGS ==========
# Create animated views
create_animation = True
# Folder containing .dat files
folderpath = "PATH_TO_FOLDER"
# Starting index
animate_start = 0
# Stopping index
animate_end   = 100
# Overwrite .png files or not
overwrite_png = True
# Duration of each frame in seconds
frame_duration = 0.2
# Remove .png files (confirmation will be asked)
remove_pngfiles = False
# Remove .npy files after finishing (confirmation will be asked)
remove_npyfiles = False


# ========== PLOTTING SETTINGS ==========
# Line of sight.
line_of_sight = "y"
# Plot on logscale or not.
logscale_halpha  = False
logscale_faraday = False
# Color map to use when plotting.
cmap_halpha  = "Reds_r"
cmap_faraday = "jet"
# Whether or not to plot magnetic field lines.
plot_Blines = False


# ========== FONT SIZES AND FORMATTING ==========
# Color bar formatting (not for log scale)
cbar_format_faraday  = "%0.2e"
cbar_format_halpha   = "%0.3e"
# Font size color bar
cbar_fontsize = 12
# Font size axis
axis_fontsize = 12
# Font size title
title_fontsize = 14


# ========== CALCULATION SETTINGS ==========
#Calculate H-alpha view.
halpha = True
#Calculate Faraday effect.
faraday = True
# If false, use RM. If true, use RM*lambda**2
include_lambda = False
faraday_lambda = 20.5  # cm


# ========== SAVE SETTINGS ==========
#Save files or not. Default = True
saveFiles = True


# ========== MULTIPROCESSING SETTINGS ==========
#Use multiprocessing to interpolate blocks.
multiple_procs = True
# Number of processors
nb_of_procs    = 4

