# Reading .dat files

[TOC]

# Scripts {#dat_reader_scripts}

All scripts needed to read MPI-AMRVAC `.dat` files in Python are in the `tools/python/read_datfile` folder.
The reader is compatible with both Python 2.7 and Python 3.7.
The files contained in the folder are modified versions of the ones in the `RayTracing` package, with dependencies on other files removed.
There are four files in total: `conversions.py`, which contains methods to calculate the primitive variables from the conservative ones. The methods to extract the header, block tree, perform regridding, etc. are contained in `dat_reader.py`.
The file `reduce_data.py` is responsible for processing the `.dat` file, in essence a class object is created which contains header information and both the conservative and primitive variables.
The main file is `read_file.py`, which contains the main method. All extensions or data analysis go in here.

# Usage {#dat_reader_usage}
Simply import `reduce_data.py` in the main file and call

    data = reduce_data.ProcessData(filename, nbprocs)
    
The parameter `filename` (the path to the `.dat` file) _must_ be given. The second parameter, `nbprocs`, is optional, and corresponds to the number of processors to use when performing regridding of the data.
If the latter is not given, the code will default to the use of all available processors.
Once the data object is created, required variables can simply be called by doing

    temperature = data.T
    density     = data.rho
    bx          = data.b1

and others. Performing regridding when loading in AMR datasets is done by default. If this is not desired, one can also iterate over the AMR block tree by calling `get_block_data(dat)`,
where the argument `dat` is the `.dat` file opened in binary mode (which can be done through `dat_file = open(filename, "rb")`). This will return a dictionary containing every block in the AMR tree, on which data analysis can be performed.
The header is retrieved in `get_header(dat)`, and returns the `.dat` file header as a dictionary. For example, if in MHD, calling `hdr["physics_type"]` will return the string `"mhd"`.
Every method is documented with docstrings in the code itself if more information is needed.






