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
Saving of these files can be disabled, if needed, by setting `saveFiles=False` in `settings.py`.

The header information is also extracted from the `.dat` file, and used to initialize the data object. This ensures that the code works for a general MPI-AMRVAC `.dat` file, for any given number of dimensions in both HD and MHD.
All the information contained in the header, together with the raw data itself (both the conservative _and_ primitive variables), is available at all times inside the data object.

Once the data is obtained, the ionization degree at each grid point is interpolated using the local pressure and temperature.
As the interpolation is done for _each_ grid point, this can take a while for large (3D) data sets. The interpolated data can again be saved for later use.

The ionization degree is then used to calculate eg. the H-alpha absorption coefficient (kappa) using the results from the aforementioned paper. In three dimensions, the kappa-matrix is integrated along a given line of sight to produce a single 2-dimensional image. Finally, this result is used in the radiative transfer equation to obtain the H-alpha line intensity. Only line-of-sight views along each coordinate axis are supported at the moment, but this is usually sufficient.
The Faraday effect can also be calculated, then the number density over the mesh is used together with the component of the magnetic field parallel to the line of sight, which is then integrated.

Letting the code save the `Numpy` files has multiple advantages. First, the `load` command from Python's `Numpy` package can load these data sets almost instantaneous, such that re-runs (for example when creating multiple figures for the same data set) are much, much faster.
Second, these `.npy` files can easily be loaded in with stand-alone scripts, and can be used for other kinds of data analysis. The size of these files is usually less than that of the `.dat` file itself.


# Usage {#raytracing_code_usage}
The user must first define the appropriate parameters in `settings.py`, such as the unit normalizations. In principle these latter ones can be trivial, and are solely used for quantifying the intensity.
Running the code is possible from inside an IDE (provided everything is passed to the `settings.py` file) or from the command line via

    python main.py -d DAT_FILE -ns -np x

where `-d` is followed by the `.dat` file. The optional argument `-ns` stands for `no_save`. If this is provided, the code will not save any data during runtime. If no arguments are given, the code defaults to the values in `settings.py`.
The code supports parallelization to perform the grid refinement if required. This option can be activated by adding the optional argument `-np x`, where `x` is the number of processors to use. Each block that has to be regridded will then be added
to a pool of workers (the processors), and each processor will refine a new block if it becomes available.

This is a standalone tool, as there is no dependence on the `amrvac` folder itself. The `RayTracing/src` folder can be copied to any directory, as long as the underlying structure in `src` is not modified.

# Variables in settings.py {#settings_information}
The `settings.py` file has to be modified before use. This can be done by changing the desired parameters, as not all of them can be called from the command line. After modifying, save the file and run the program (see Usage).

## Normalizations
variable name  |  type   |   description
--- | --- | ---
`unit_length`           |  float    |   Unit normalization for the length scale. Should be in cm.
`unit_temperature`      |  float    |   Unit normalization for the temperature. Should be in K.
`unit_numberdensity`    |  float    |   Unit normalization for the number density. Should be in cm-3.

## User parameters
variable name  |  type   |   description
--- | --- | ---
`altitude`              |  integer  |   Altitude above the solar disk. Should be in km.
`type_simulation`       |  string   |   The type of simulation. Possibilities are `prominence` or `filament`.
`theta`                 |  float    |   The angle between the x-y plane and z-axis. Only used when `plot_Blines = True`, draws straight lines under this angle.
`phi`                   |  float    |   Same as `theta`, this is the angle from the x-axis towards the y-axis. These angles are only useful with a straight, homogeneous magnetic field. Should be left off by default.

## Animation settings
variable name  |  type   |   description
--- | --- | ---
`create_animation`      |  boolean  |   If `True`, searches for .dat files specified in `folderpath` for a given range (see `start` and `end`). Loads the .dat files, creates .png files and finally merges them into a single .gif file.
`folderpath`            |  string   |   Path to the folder containing the .dat files needed to create the animation.
`animate_start`         |  integer  |   Index at which the animation should start (eg for 20, the file `xxx0020.dat` is the starting file).
`animate_end`           |  integer  |   Index at which the animation should stop.
`overwrite_png`         |  boolean  |   If `True`, overwrites all presently saved .png files in `animations/png_files` for the chosen simulation. If `false`, skip the files that are already present.
`frame_duration`        |  float    |   Duration of each frame in the .gif file, in seconds.
`remove_pngfiles`       |  boolean  |   If `True`, removes all the .png files that were created after finishing the .gif animation. Confirmation will be asked before removal starts.
`remove_npyfiles`       |  boolean  |   If `True`, removes all .npy files in interpolated_files after finishing the .gif animation. Confirmation will be asked before removal starts. 

## Plotting settings
variable name  |  type   |   description
--- | --- | ---
`line_of_sight`         |  string   |   The line of sight to integrate over. Possible values are `x`, `y`, `z` and `all`. At the moment the `all` parameter is not supported if `create_animation = True`.
`logscale_halpha`       |  boolean  |   If `True`, plots the h-alpha view on a logscale.
`logscale_faraday`      |  boolean  |   If `True`, plots the Faraday view on a logscale.
`cmap_halpha`           |  string   |   Colormap to use when plotting h-alpha views. `Reds_r` is suggested.
`cmap_faraday`          |  string   |   Colormap to use when plotting Faraday views, `jet` is suggested.
`plot_Blines`           |  boolean  |   If `True`, plots straight magnetic field lines over the grid according to the angles `theta` and `phi.`

## Font sizes and formatting
variable name  |  type   |   description
--- | --- | ---
`cbar_format_faraday`   |  string   |   Format parameter to use for the Faraday view colorbar (only on non-logscale). `%0.2e` is suggested.
`cbar_format_halpha`    |  string   |   Format parameter to use for the h-alpha view colorbar (only on non-logscale). `%0.2e` is suggested.
`cbar_fontsize`         |  integer  |   Fontsize to use for the color bar ticks.
`axis_fontsize`         |  integer  |   Fontsize to use for the axes labels.
`title_fontsize`        |  integer  |   Fontsize to use for the title.

## Calculation settings
variable name  |  type   |   description
--- | --- | ---
`halpha`                |  boolean  |   If `True`, calculates and plots the h-alpha view.
`faraday`               |  boolean  |   If `True`, calculates and plots the Faraday view.
`include_lambda`        |  boolean  |   If `False`, shows the Rotation Measure (RM) in rad/m2 for the Faraday rotation on the plot. If `True`, shows lamda^2 * RM on the plot with lambda user-defined.
`faraday_lambda`        |  float    |   Wavelength to use if `include_lambda = True`. Should be in cm.

## Others
variable name  |  type   |   description
--- | --- | ---
`filename`              |  filename |   The filename of the .dat file to read.
`saveFiles`             |  boolean  |   If `True`, save the interpolated and regridded results as Numpy (.npy) files in `interpolated_files` and `dat_files`, respectively. If `False`, interpolation and regridding has to be redone each run.
`multiple_procs`        |  boolean  |   If `True`, activates the `multiprocessing` module in Python to perform the regridding.
`nb_of_procs`           |  integer  |   Number of processors to use when the multiprocessing module is activated.



