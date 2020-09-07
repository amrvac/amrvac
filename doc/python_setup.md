# Setting up the python tools

[TOC]

# Introduction {#introductionpython}

This document describes the use of Python for processing MPI-AMRVAC simulation data.
Both the reading of native `.dat`- and `vtk`-filetypes (`.vtu`, `.pvtu`, `.vti`) are supported.
The goal of these Python tools is to provide additional data-analysis support with main focus
on 1D and 2D simulation data, although there is functionality to read in and process 3D datasets.

All Python tools are located in the directory `$AMRVAC_DIR/tools/python`. This folder contains a setup
script to install the package, a small readme with some instructions and finally the subdirectory `amrvac_pytools` itself,
containing the Python tools for both the .dat files and .vtk files.

These tools were originally initiated by _Oliver Porth_ (in particular the `vtk` readers),
and have recently been refactored and extended to include datfile support by _Niels Claes_ and _Clément Robert_.

# Installation {#installation}
## Prerequisites

The Python tools depend on a few standard packages:

* [SciPy](https://www.scipy.org/) which bundles essential numerical tools
* [matplotlib](https://matplotlib.org/) as a plotting backend
* [IPython](https://ipython.org/) which provides a command-line interface for interactive use
* [vtk](https://pypi.org/project/vtk/) for image processing and visualisation

## Setting up
The tools are packed into a handy Python package, which is continuously further developed and improved.
We recommend installing it with _develop_ arguments, so that it does not need to be reinstalled with each update.
This can be done in three distinct ways:

1. Using `conda` (**recommended**)

        cd $AMRVAC_DIR/tools/python
        conda install conda-build
        conda develop .
	
2. Using `pip`

        cd $AMRVAC_DIR/tools/python
        pip install -e .
	
3. Simply vanilla python

        cd $AMRVAC_DIR/tools/python
        python setup.py develop
	
## Alternative method
Ofcourse, adding the folder directly to your PYTHONPATH works as well but is not recommended.

1. _Expanding the PYTHONPATH environment variable_

   * [Mac, Linux] Simply type the following in your bash shell:

         export PYTHONPATH=$AMRVAC_DIR/tools/python:$PYTHONPATH

   * [Windows] Windows users can modify the environment variables using the control panel:

         Start > Control panel > Edit the system environment variables >
         tab 'advanced' > environment variables

     Create or edit the variable `PYTHONPATH` in the _system variables_ list, where you add the
     full path to `.../amrvac/tools/python`. Multiple paths for a single environment variable are separated by a semicolon ';'.
