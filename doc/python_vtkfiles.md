# Reading the vtk files

[TOC]

# Introduction on reading vtk files {#introductionvtk}

This document describes the use of Python for plotting MPI-AMRVAC 1D and 2D simulation data using the `vtk` filetypes.
All required packages and instructions on how to set up the tools folder are described in @ref python_setup.md.
The vtk-file tools can be found in the folder `$AMRVAC_DIR/tools/python/vtkfiles`.

# Reading datasets {#reading_vtk}
Currently only the reading of **cell-centered** `vtk` file-types (`.vtu`, `.pvtu`, `.vti`) is supported,
together with some ascii `.csv` files. Reading the data is handled in the file `read.py`, present in the `vtkfiles` subdirectory.
This file contains various classes for different data.

Reading the `vtk` files is best done using an interactive environment. Hence, navigate to your directory containing your data and
fire up `IPython`:

    ipython --pylab

Next, import the reading and plotting classes:

    from amrvac_pytools.vtkfiles import read, amrplot

## vtu and pvtu filetypes
These filetypes use the class `read.load()`, loading the data can be done as follows:

    offset=15
    ds = read.load_vtkfile(offset, file='KH', type='vtu')

which will read the file `KH0015.vtu`. The argument `offset` is zero-padded to four digits
and is used to specify the file number. The argument `type` is either '`vtu`' (default) or `'pvtu'`,
which is used to append the file extension and select the correct reading methods.
The variables are contained in the instance `ds`, thus

    ds.rho

holds the Numpy array for the density, `ds.v1` holds the first component of the velocity
(assuming the data is in primitive quantities and your velocity variable is called '`v1`') and so on.
The available variables will be printed to the terminal when loading in the dataset.
Different datasets can be loaded by supplying different offsets.

## vti filetypes
This filetype uses the class `read.loadvti()`, loading the data can be done as follows:

    offset=15
    ds = read.loadvti(offset, file='KH')

which in turn reads the file `KH0015.vti`. In this case no filetype has to be specified.
Accessing the variables in these uniform `.vti` files is done in a similar way as described above,
the available variables are again printed to the terminal when loading the dataset.

## Particle csv files
Loading particle ensemble csv (comma-separated-value) files uses the class `read.ensemble()`.
To load the file `data_ensemble000015.csv`, one can do

    offset=15
    ds = read.ensemble(offset, file='data', npayload=1)

In this case, `offset` is zero padded to six digits. The argument `npayload` stands for the
number of payload arrays in the data (default is one).

## Slice csv files
To load a 1D slice in `.csv` format from a 2D simulation, you can use

    offset=15
    ds = read.loadcsv(offset, file='data', dir=1, coord=0.)

which loads in `data_d1_x+0.00D00_n0015.csv`. The argument `offset` refers to
the last 4 digits of the filename, the argument `coord` translates to `0.00D00` and
`dir=1` refers to the `d1` part of the filename (setting `dir=2` would mean `data_d2_x...`).

# Plotting data {#plotting_vtk}
## 1D plotting

Once the data is loaded into the instance `ds`, you can access it with
`ds.varname` where `varname` is a variable name chosen in your
AMRVAC parameter-file.  Pythons powerful introspection makes it really
easy to see what is available in the data: just type `ds.` and hit the *tab* key to see the
available variable names.

For uniform data (`.vti`), the arithmetic centre of the cells can be obtained via the command

    ds.getCenterPoints()

Hence, a simple 1D plot can be created using matplotlib by doing

    import matplotlib.pyplot as plt
    plt.plot(ds.getCenterPoints(), ds.rho)
    plt.show()

At this point, `ds.getCenterPoints()` and `ds.rho` are simple Numpy arrays, so you can use all of matplotlib for
your plotting experience.

For the structured `.vti` files, the attribute `ds.x` (or `.y` or `.z`) denotes the cell centre coordinate, such
that a plot can be obtained with

    import matplotlib.pyplot as plt
    plt.plot(ds.x, ds.rho)
    plt.show()

## 2D plotting
The plotting of two-dimensional AMR data is also supported. These plotting tools come in two
flavours:
1. The class `amrplot.polyplot` to show each cell in one patch
2. The class `amrplot.rgplot` to show data interpolated onto a uniform mesh

Both classes offer the same functionality. In order to display a 2D image of the dataset `ds`, invoke a new instance

    p1 = amrplot.polyplot(ds.rho, ds)

This will open a new window, containing a plot of the requested data. Drawing the window for the first time might take a few seconds,
as a list of the cell patches has to be generated. If all of this is done in interactive mode another quantity can be
plotted using the same window (for example the temperature) with the command

    p1.show(ds.p/ds.rho, ds)

The class `amrplot.polyplot` has several convenient members. One example is `amrplot.polyplot.save`, which saves the image to a file and takes the filename as parameter.
The axis ranges can be fixed with the attributes `p1.xrange` and `p1.yrange`, while the window zoom can be controlled by setting `p1.fixzoom` to a value other than `None`.
Additionally there is a nice interactive feature when using `amrplot.polyplot`: if the middle mouse button is clicked to mark a cell in the plotting window, a pink crosshair
will be displayed and all variables at the chosen location will be returned to the terminal. This may take a few seconds the first time this feature is executed,
as the centerpoints must be calculated first.

# Convenient auxiliary functions {#vtk_auxiliary}
The class `amrplot` has some additional auxiliary functions which can be quite useful:
* `amrplot.line`: selects variables over a straight path across a 2D slice
* `amrplot.velovect`: provides velocity vectors for a flow field
* `amrplot.contour`: overlays contours onto the plotting window
* `amrplot.streamlines`: shows streamlines of a vector field
* `amrplot.streamline`: plots a single streamline with given starting values

See the respective methods and classes in `vtkfiles/amrplot` for more information.
