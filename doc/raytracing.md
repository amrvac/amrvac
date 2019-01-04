# Ray-Tracing Algorithm

[TOC]

# Introduction {#raytracing_intro}

The source code itself is written in Python and can be found in the `tools/python/RayTracing/src` folder.
It is compatible with both Python 2.7 and Python 3.7.
The currently implemented tables and methods are based on results found in
[Heinzel et al. (2015)](https://www.aanda.org/articles/aa/pdf/2015/07/aa25716-15.pdf "Paper" ).

# Structure {#raytracing_code_structure}

The code is written with possible future extensions in mind. This means that most methods are
contained within their respective package, and can be called easily from other scripts.
When calling `main.py`, the code first checks the existence of a folder called `interpolated_files`, which is used to store interpolated data obtained during runtime for possible re-runs. If this folder does not exist, a new one is created.

The code then creates a `ProcessData` class object, which eventually contains all the header and data information from the given `.dat` file. This object can easily be passed on to other methods, eliminating the need for various global variables.

Before this object can be created however, the data must first be obtained from the `.dat` file,
which is done in `src/dataIO/dat_reader.py`. A check is performed to see if the entire mesh is uniform. If this is not the case, a _regridding_ is performed to the finest mesh resolution using linear interpolation on each grid block. Depending on the size of the mesh and the number of dimensions this can take a while. When the regridding is finished, the code stores the regridded data in a Numpy file in the folder `dat_files`. This file is loaded when doing re-runs, such that the regridding only has to be done once per data set.

The header information is also extracted from the `.dat` file, and used to initialize the data object. This ensures that the code works for a general MPI-AMRVAC `.dat` file, for any given number of dimensions in both HD and MHD.
All the information contained in the header, together with the raw data itself (both the conservative _and_ primitive variables), is available at all times inside the data object.

Once the data is obtained, the ionization degree at each grid point is interpolated using the local pressure and temperature.
As the interpolation is done for _each_ grid point, this can take a while for large (3D) data sets. The interpolated data is again saved for later use.

The ionization degree is then used to calculate the H-alpha absorption coefficient (kappa) using the results from the aforementioned paper. In three dimensions, the kappa-matrix is integrated along a given line of sight to produce a single 2-dimensional image. Finally, this result is used in the radiative transfer equation to obtain the H-alpha line intensity. Only line-of-sight views along each coordinate axis are supported at the moment, but this is usually sufficient.


# Usage {#raytracing_code_usage}
The user must first define the appropriate parameters in `settings.py`, such as the unit normalizations and magnetic field angle with the axes (if in mhd).
Running the code is possible from inside an IDE (provided everything is passed to the `settings.py` file) or from the command line via

    python main.py -d DAT_FILE -ns

where `-d` is followed by the `.dat` file. The optional argument `-ns` stands for `no_save`. If this is provided, the code will not save any data during runtime. If no arguments are given, the code defaults to the values in `settings.py`.

This is a standalone tool, as there is no dependence on the `amrvac` folder itself. The `RayTracing/src` folder can be copied to any directory, as long as the underlying structure in `src` is not modified.


