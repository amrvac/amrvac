# Setting up the python tools

[TOC]

# Introduction {#introduction}

This document describes the use of Python for processing MPI-AMRVAC simulation data.
Both the reading of native `.dat`- and `vtk`-filetypes (`.vtu`, `.pvtu`, `.vti`) are supported.
The goal of these Python tools is to provide additional data-analysis support with main focus
on 1D and 2D simulation data, although there is functionality to read in and process 3D datasets.

All Python tools are located in the directory `$AMRVAC_DIR/tools/python`, while the underlying
subdirectories `datfiles` and `vtkfiles` contain the tools for their respective filetypes.

These tools were originally initiated by _Oliver Porth_ (in particular the `vtk` readers),
and have recently been refactored and extended to include datfile support by _Niels Claes_.

# Installation {#installation}
## Prerequisites

The Python tools depend on a few standard packages:

* [SciPy](https://www.scipy.org/) which bundles essential numerical tools
* [matplotlib](https://matplotlib.org/) as a plotting backend
* [IPython](https://ipython.org/) which provides a command-line interface for interactive use
* [vtk](https://pypi.org/project/vtk/) for image processing and visualisation

Additionally, these tools are developed with Python 3 in mind, so it is possible that Python 2
users may experience limited functionality.

## Setting up
The tools are installed by adding the directory `$AMRVAC_DIR/tools/python` to your python path.
This can be done in a number of ways:

1. _Expanding the PYTHONPATH environment variable_

   * [Mac, Linux] Simply type the following in your bash shell:

         export PYTHONPATH=$AMRVAC_DIR/tools/python:$PYTHONPATH

   * [Windows] Windows users can modify the environment variables using the control panel:

         Start > Control panel > Edit the system environment variables >
         tab 'advanced' > environment variables

     Create or edit the variable `PYTHONPATH` in the _system variables_ list, where you add the
     full path to `.../amrvac/tools/python`. Multiple paths for a single environment variable are separated by a semicolon ';'.
2. _Adding the path directly in your script_

   Add these two lines to the beginning of every Python script where the tools are to be used:

       import sys
       sys.path.append('../amrvac/tools/python')

   where you replace '`..`' by the full path to the folder.
