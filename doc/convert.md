# Data file conversion

[TOC]

# Introduction {#convert_intro}

The standard [MPI-AMRVAC dataformat](fileformat.md), i.e. the _*.dat_
files usable for restart, contain all the conservative variables in
all gridblocks, and hence suffice for visualization.
> Since late 2019, it can be directly read, visualised and
> analysed with the Python package [yt](https://yt-project.org),
> see [our documentation to get started with yt](yt_usage.md).

However, in many instances, one would like to use data formats that
are directly readable by some of the more widespread visualization
software packages. Therefore, we created the _convert.t_ module, which
ensures that this post-process data file conversion can be done with
the same executable (but possibly even on a different platform). The
many possibilities include conversion to _*.vtu_ (*VTK* *U*nstructured
data format) directly readable by [Paraview](http://www.paraview.org/)
(or [ViSiT](https://wci.llnl.gov/codes/visit/)), to _*.plt_ format for
the commercial package [Tecplot](http://www.tecplot.com/).  Also,
**this info will not explain you how to use the mentioned software for
visualization, but just explain how to do the conversion.**
Furthermore, this part of the code is subject to continuous change and
improvement, and we welcome extra contributions.

We now give some brief info on how to use the same executable _amrvac_ (which
you already compiled and used to obtain output _*.dat_ files with), to convert
a single or all _*.dat_ file(s) to one of these formats.

# Converting (on a single CPU) {#converting}

** Note that all the steps below assume you're running on a single CPU. The same steps are to be taken for obtaining any of the other precoded data formats. One important warning is due: when you run a simulation for some reason twice, and you did not erase the previously created _*.dat_ files, these files are overwritten (if the base_filename has not changed). Then, it may be that conversion fails, since the end of the file may contain some leftover data from the previous time, if the filelength has changed due to some other reason. The only remedy to this is that one should always remove old _*.dat_ files, or never forget to change the name for the files accordingly, by setting _base_filename_ in the _&amp;filelist;_.**

We will assume that you ran the standard 2D advection problem used for test
purposes, i.e. that you did the following steps beforehand:

    cd $AMRVAC_DIR/tests/rho/vac
    $AMRVAC_DIR/setup.pl -d=2
    make
    mpirun -np 1 amrvac

We also assume that in the parameter file amrvac.par, the namelist
_&amp;filelist;_ was stating (note that the end of the namelist is indicated
as usual by a backslash)

     &filelist
            base_filename='vaclogo'
     /

If all went well, you then have created as many _*.dat_ files as requested
through the settings you provided in the combined _&amp;savelist;_ and
_&amp;stoplist;_ namelists from the [par-file](par.md). For the example,
they normally default to asking a full data dump at time zero, as well as
every time the time has increased by 0.05, and this till _tmax=1.0d0_, such
that we actually have 21 snapshots in total. You should thus have files like
_vaclogo0000.dat_ up to _vaclogo0020.dat_. You can now
individually convert such _*.dat_ file to a _*.vtu_ file by doing the
following. Edit the amrvac.par file, to select a visualization data format like

     &filelist
            base_filename='vaclogo'
            convert_type='vtuBCC'
            convert=.true.
            restart_from_file='vaclogo0000.dat'
     /

you can then convert the single _vaclogo0000.dat_ file simply
running again

    mpirun -np 1 amrvac

or, which is actually equivalent (single CPU)

    amrvac

Note that this will create a new file, namely _vaclogo0000.vtu_, which
can be directly imported in Paraview. It will, under the settings above, just
contain the density on the grid hierarchy at time zero. The
_convert_type='vtuBCC'_ indicates that the data is exactly as the code
interprets and updates the values, namely as cell-centered quantities. 

Realizing that you typically want to convert multiple data files, you can do
this by repeating the above as many times as there are _*.dat_ files, by
raising/changing the _restart_from_file_ identifier. Since you typically want to
convert all data files between a minimum and maximum number of similarly named
files, the script **aiconvert** is added. If you have a line
`PATH="$AMRVAC_DIR:$AMRVAC_DIR/tools:./:$PATH"` in `~/.bash_profile` (or
`~/.bashrc`), typing _aiconvert_ will tell you its intended usage.
In the example case at hand, where we created 21 data files from running the
advection problem, **this _aiconvert_ script needs the intended _base_filename_ 
and the executable _amrvac_ to exist in the same directory.** It will complain when the parfile
does not exist, and obviously requires the existence of all files between the
start and stopindex (0 and 20 here). With paraview, you will then be able to
immediately import all 21 _*.vtu_ files with the same base filename, and
directly make movies or still images from them.

For example, to convert snapshots from number 10 to number 20:

    aiconvert 10 20

or to convert the snapshot number 12

    aiconvert 12

or just type

    aiconvert

to convert all the snapshots! You can also specify a parfile other than amrvac.par as:

    aiconvert newname.par 0 20

or to convert the snapshot number 12

    aiconvert newname.par 12

For details of aiconvert, please read the header of the
$AMRVAC_DIR/tools/aiconvert.

# Parallel conversion options {#parallel_convert}

For very large simulations (typically 3D, and/or runs achieving high effective
resolutions), even the data conversion may need to be done in parallel (and
ultimately, the visualization itself too). The _convert.t_ allows to perform
some of the _*.dat_ conversions in parallel, in particular, this is true for
the _*.vtu_ format, and for the _*.plt_ format. You should then select one of

    convert_type='vtumpi'
    convert_type='vtuCCmpi'
    convert_type='vtuBmpi'
    convert_type='vtuBCCmpi'
    convert_type='vtimpi'
    convert_type='vtiCCmpi'
    convert_type='pvtumpi'
    convert_type='pvtuCCmpi'
    convert_type='pvtuBmpi'
    convert_type='pvtuBCCmpi'
    convert_type='tecplotmpi'
    convert_type='tecplotCCmpi'

Here, the prefix _p_ stands for the parallel file format, where each process
is allowed to dump its data into a single (e.g. _*.vtu_) file and a master
file (e.g. _*.pvtu_) is stored by rank zero. This has the advantage that the
write operation on suitable file systems is sped up significantly. In a
visualization software, only the _*.pvtu_ files need to be imported and also
the reading process is sped up in case of parallel visualization.

You can again use aiconvert as explained above, and type in the number
of processors to use by answering a popup question:

    How many processors do you want to use? (default=1) 4

# Autoconvert {#autoconvert}

In addition to the conversion after the run, AMRVAC now offers also to
directly output files ready to use for visualization along with the
simulation. A parallel run will however only be capable to provide the file-
types ready for parallel conversion (see parallel conversion). To enable this
capability, simply set the switch _autoconvert=.true._. The example above
would then read

    &filelist;
            base_filename='vaclogo'
            autoconvert=.true.
            convert_type='vtuBCCmpi'
     /

and when the code is run via

    mpirun -np 2 amrvac

three new output files (_vaclogoXXXX0000.vtu,
vaclogoXXXX0001.vtu_) will appear simultaneous to the _vaclogoXXXX.dat_
files, stored at given intervals. All functionality of the usual convert is
conserved, e.g. derived quantities and primitive variables (using the
_saveprim=.true._ option) can be stored in the output files.

# Notes on conversion possibilities {#notes}

## Cell center versus cell corner values {#cell_vs_corner}

In all cases mentioned below, the difference between convert-types with or
without _CC_ relates to the difference between cell center (`CC') versus cell
corner values. For the cell center case, no interpolation of the computed data
is needed, since the (conservative) variables actually represent volume
averages stored at cell centers. For the other cases (without 'CC'), the
_convert.t_ tries to interpolate from cell center to the cell corners, and
then the data will be known at the grid nodes. This will be noticable on
reading in the data in _paraview_, which reports whether data is cell data
(cell centered) or point data (nodes or corners). In principle, the conversion
from cell centered (or cell data) to cell corner (or point data) types can
also be achieved in paraview itself, with the filter _Cell data to Point
data_. There may be subtle differences between the way MPI-AMRVAC does this
interpolation, and the way it happens internally to paraview, so we provide
both options as output _*.vtu_ files. Similar observations hold for the
Tecplot format.

## Conservative/primitive storage and adding derived quantities {#cons_vs_prim}

The **saveprim** logical allows you to select whether the conservative or
primitive variables need to be stored in the resulting output file. The names
for the conservative variables and primitive ones are hard coded depending on
the physics.

Another very useful option is to specify which variables actually need to be
converted: by default all conservative variables available in the _*.dat_ file
will be included, but then again filesizes may become restrictive. For that
purpose, the logical array _w_write_ allows to select which variable(s) to
store (and this in combination with saveprim, possibly). You can then create
different files for selected variables, knowing that the output filename will
start with _base_filename_.

We allow the possibility to compute derived variables from the _*.dat_ file in
the userfile, by setting how many you add beyond the _nw_ variables typcial
for the physics module at hand, in the integer _nwauxio_. Correspondingly that
many variables, you should then compute and store in the _w(*,nw+1)_ ...
_w(*,nw+nwauxio)_ entries, in the user-written subroutine _ specialvar_output_
(as defined in _mod_usr_methods.t_). The names for these variables then
need to be provided in the corresponding _specialvarnames_output_ subroutine.
This feature is very useful, for the same reason as above: you can let the 
code compute gradients of scalar fields, divergence of vector quantities, curls 
of vectors, etc, using the precoded subroutines for that purpose found in _geometry.t_.
You then do not rely on visualization software to do interpolations or
discretizations, which may not reflect those actually taken in MPI-AMRVAC.

Another useful feature is the possibility to select the output AMR level. You
can let the code compute from the _*.dat_ file an output residing on a
specified level _level_io_. This then uses the MPI-AMRVAC internal means to
perform restriction and prolongations, and you can then make sure to have a
single uniform grid output too.

### _convert_type='vtu'_ or _convert_type='vtuCC'_

Does the conversion to _*.vtu_ data files. This option works on 1 CPU. The
resulting _*.vtu_ files contain data in ASCII format.

### _convert_type='vtuB'_ or _convert_type='vtuBCC'_

Does the conversion to _*.vtu_ data files. This option works on 1 CPU. The
resulting _*.vtu_ files are in binary format.

### _convert_type='vtumpi'_ or _convert_type='vtuCCmpi'_

Does the conversion to _*.vtu_ data files. This option works on multiple CPUs.
The resulting _*.vtu_ files contain data in ASCII format.

### _convert_type='vtuBmpi'_ or _convert_type='vtuBCCmpi'_

Does the conversion to _*.vtu_ data files. This option works on multiple CPUs.
The resulting _*.vtu_ files contain data in binary format. (recommended)

### _convert_type='tecplot'_ or _convert_type='tecplotCC'_

This realizes conversion to _*.plt_ files, which can be read in directly by
Tecplot. Note that the output is in ASCII format, which may mean huge data
sets, but Tecplot has its own **preplot** command that will allow you to
operate on such a file, and thereby make a binary version for it. The above is
for single CPU execution, allows to add user-defined variables with _nwauxio_,
and renormalization using the _normt_ and _normvar_ array.

### _convert_type='tecplotmpi'_ or _convert_type='tecplotCCmpi'_

Same as above, but allows to perform the conversion in parallel. One can add
user -defined variables with _nwauxio_, and renormalize using the _normt_ and
_normvar_ array. The current implementation is such that tecplotmpi and
tecplotCCmpi will create different length output ASCII files when used on 1
versus multiple CPUs. In fact, on 1 CPU, there will be as many (tecplot) zones
as there are levels, while on on multiple CPUs, there will be a number of
zones up to the number of levels times the number of CPUs (or less, when some
level does not exist on a certain CPU).

### onegrid(mpi), oneblock(B), ...

Extra possibilities to allow creation of a single uniform grid level output.
Please inspect the workings in _convert.t_.
