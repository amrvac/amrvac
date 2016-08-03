# Data file conversion

[TOC]

# Introduction {#convert_intro}

The standard [MPI-AMRVAC dataformat](fileformat.html), i.e. the _*.dat_ files
usable for restart, contain all the conservative variables in all gridblocks,
and hence suffice for visualization, in principle. However, in many instances,
one would like to use data formats that are directly readable by some of the
more widespread visualization software packages. Therefore, we created the
_convert.t_ module, which ensures that this post-process data file conversion
can be done with the same executable (but possibly even on a different
platform). The many possibilities include conversion to _*.vtu_ (VTK
unformatted data format) directly readable by
[Paraview](http://www.paraview.org/) (or
[ViSiT](https://wci.llnl.gov/codes/visit/)), to _*.plt_ format for the
commercial package [Tecplot](http://www.tecplot.com/), or the _*.dx_ format
for the opensource package [openDX](http://www.opendx.org/). We also provide
possibilities to convert to a format _*.out_ suitable for Idl, for which some
automated macro's can be made in analogy with those used for the [Versatile
Advection Code ](http://grid.engin.umich.edu/~gtoth/VAC), but be forewarned
that 3D (and even 2D) AMR data for Idl will require you to program your own
macro's. Also, **this info will not explain you how to use the mentioned
software for visualization, but just explain how to do the conversion.**
Furthermore, this part of the code is subject to continuous change and
improvement, and we welcome extra contributions.

We now give some brief info on how to use the same executable _amrvac_ (which
you already compiled and used to obtain output _*.dat_ files with), to convert
a single or all _*.dat_ file(s) to one of these formats.

# Converting (on a single CPU) {#converting}

** Note that all the steps below assume you're running on a single CPU. The same steps are to be taken for obtaining any of the other precoded data formats. One important warning is due: when you run a simulation for some reason twice, and you did not erase the previously created _*.dat_ files, these files are overwritten (if the filenameout has not changed). Then, it may be that conversion fails, since the end of the file may contain some leftover data from the previous time, if the filelength has changed due to some other reason. The only remedy to this is that one should always remove old _*.dat_ files, or never forget to change the name for the files accordingly, by setting _filenameout_ in the _&amp;filelist;_.**

We will assume that you ran the standard 2D advection problem used for test
purposes, i.e. that you did the following steps beforehand:

    cd $AMRVAC_DIR/tests/rho/vac
    mkdir datamr
    $AMRVAC_DIR/setup.pl -d=22 -phi=0 -z=0 -p=rho -g=16,16 -s
    make
    ln -s testrho_vac22 amrvac.par
    mpirun -np 1 amrvac

We also assume that in the parameter file mentioned above, the namelist
_&amp;filelist;_ was stating (note that the end of the namelist is indicated
as usual by a backslash)

     &filelist;
            filenamelog='datamr/vaclogo'
            filenameout='datamr/vaclogo'
            primnames='rho'
     /

If all went well, you then have created as many _*.dat_ files as requested
through the settings you provided in the combined _&amp;savelist;_ and
_&amp;stoplist;_ namelists from the [par-file](par.html). For the example,
they normally default to asking a full data dump at time zero, as well as
every time the time has increased by 0.05, and this till _tmax=1.0d0_, such
that we actually have 21 snapshots in total. You should thus have files like
_datamr/vaclogo0000.dat_ up to _datamr/vaclogo0020.dat_. You can now
individually convert such _*.dat_ file to a _*.vtu_ file by doing the
following. Edit the par-file, to modify the _&amp;filelist;_ to something like

     &filelist;
            filenamelog='datamr/vaclogo'
            filenameout='datamr/vaclogo'
            primnames='rho'
            filenameini='datamr/vaclogo'
            convert=.true.
            convert_type='vtuCC'
            saveprim=.false.
            snapshotini=0
     /

Assuming that this par-file is still known through the symbolic link
_amrvac.par_ as above, you can then convert a single _*.dat_ file (here the
_datamr/testrho/vaclogo0000.dat_ file, as we select _snapshotini=0_) simply
running again

    mpirun -np 1 amrvac

or, which is actually equivalent (single CPU)

    amrvac

Note that this will create a new file, namely _datamr/vaclogo0000.vtu_, which
can be directly imported in Paraview. It will, under the settings above, just
contain the density on the grid hierarchy at time zero. The
_convert_type='vtuCC'_ indicates that the data is exactly as the code
interprets and updates the values, namely as cell-centered quantities. The
_saveprim=.false._ has for the example here no real meaning, as for advection
conservative and primitive variables coincide (just density _rho_ exists).

Realizing that you typically want to convert multiple data files, you can do
this by repeating the above as many times as here are _*.dat_ files, by
raising/changing the _snapshotini_ identifier. Since you typicallly want to
convert all data files between a minimum and maximum number of similarly named
files, the script **doconvert** is added. If you have a line
`PATH="$AMRVAC_DIR:$AMRVAC_DIR/tools:./:$PATH"` in `~/.bash_profile` (or
`~/.bashrc`), typing _doconvert_ will tell you its intended usage, namely

    doconvert testrho_vac22 0 20

in the example case at hand, where we created 21 data files from running the
advection problem. **This _doconvert_ script does assume that you actually
edited the par-file manually once as above (such that the needed lines for
conversion are in the _&amp;filelist;_ namelist), and that the executable
_amrvac_ exists in the same directory.** It will complain when the parfile
does not exist, and obviously requires the existence of all files between the
start and stopindex (0 and 20 here). With paraview, you will then be able to
immediately import all 21 _*.vtu_ files with the same base filename, and
directly make movies or still images from them. Furthermore, there is a more
convenient script **aiconvert**, which will automatically modify the par-file,
do converting data, and resume the par-file, as long as the amrvac.par has
been linked and it has the parameters for converting like

     &filelist;
            filenamelog='datamr/vaclogo'
            filenameout='datamr/vaclogo'
            primnames='rho'
     /
            filenameini='datamr/vaclogo'
            convert=.true.
            convert_type='vtuCC'
            saveprim=.false.
            snapshotini=0

     &savelist;

For example, to convert snapshots from number 10 to number 20:

    aiconvert 10 20

or to convert the snapshot number 12

    aiconvert 12

or just type

    aiconvert

to convert all the snapshots! Besides, if amrvac.par is not linked just do
similarly as doconcert:

    aiconvert testrho_vac22 0 20

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

Also, you can then use the same strategy as explained above for converting on
a single CPU: you will always need to edit the [par-file](par.html#Filelist)
once to specify how to do the conversion, and then you may run interactively
on e.g. 4 CPU like

    mpirun -np 4 amrvac

or do this in batch (use a batch job script for that), to do multiple data
file conversions. We also provide a small script, called **doconvertpar**,
which works similar to the **doconvert** explained above, but takes one extra
parameter: the number of CPUs. Its usage is described by

    doconvertpar parfilename startindex stopindex nprocessor

Besides, you can again use aiconvert as expained above, and type in the number
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
            filenamelog='datamr/testrho/vaclogo'
            filenameout='datamr/testrho/vaclogo'
            primnames='rho'
            saveprim=.false.
            autoconvert=.true.
            convert_type='pvtuCCmpi'
     /

and when the code is run via

    mpirun -np 2 amrvac

three new output files (_vaclogoXXXX.pvtu, vaclogoXXXXp0000.vtu,
vaclogoXXXXp0001.vtu_) will appear simultaneous to the _vaclogoXXXX.dat_
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
nterpolation, and the way it happens internally to paraview, so we provide
both options as output _*.vtu_ files. Similar observations hold for the
Tecplot format.

## Conservative/primitive storage and adding derived quantities {#cons_vs_prim}

The **saveprim** logical allows you to select whether the conservative or
primitive variables need to be stored in the resulting output file. The names
for the conservative variables are taken from the _wnames_ string, and those
for the primitive need to be set in _primnames_.

Another very useful option is to specify which variables actually need to be
converted: by default all conservative variables available in the _*.dat_ file
will be included, but then again filesizes may become restrictive. For that
purpose, the logical array _writew_ allows to select which variable(s) to
store (and this in combination with saveprim, possibly). You can then create
different files for selected variables, knowing that the output filename will
start with _filenameout_, while the actual data file converted is known from
the combination _filenameini_ and _snapshotini_.

We allow the possibility to compute derived variables from the _*.dat_ file in
the userfile, by setting how many you add beyond the _nw_ variables typcial
for the physics module at hand, in the integer _nwauxio_. Correspondingly that
many variables, you should then compute and store in the _w(*,nw+1)_ ...
_w(*,nw+nwauxio)_ entries, in the user-written subroutine _ specialvar_output_
(as defined in _amrvacnul.speciallog.t_). The names for these variables then
need to be provided in the corresponding _specialvarnames_output_ subroutine,
which simply then extends the strings _wnames_ and _primnames_. This feature
is very useful, for the same reason as above: you can let the code compute
gradients of scalar fields, divergence of vector quantities, curls of vectors,
etc, using the precoded subroutines for that purpose found in _geometry.t_.
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

### _convert_type='idl'_ or _convert_type='idlCC'_

This will do the conversion to _*.out_ files, which are a generalization of
the [Versatile Advection Code ](http://grid.engin.umich.edu/~gtoth/VAC) output
files. For those VAC-style files, extensive macro's are provided with the VAC
code itself, allowing for fairly interactive visualization of quantities, or
computation of derived quantities etc. The Idl macro's that allow to read
_*.out_ files converted with MPI-AMRVAC, and that use similar commands as for
VAC files, are downloadable [here as a single gzipped tar
file](allidlmacros.tar.gz). It contains the hidden, to be adjusted file
_.idlrc_, and the directories _Idl_ and _Idl_amr_, with macro's inside.
However, they only allow a very limited visualization, with some (possibly
incomplete and inaccurate) [description here](idl.html), fine for 1D and
little for 2D runs, but with no support for 3D data analysis. This would mean
you need to write your own Idl macro's, after you fully understand what
dataformat is actually stored in a _*.out_ file. For that, just study the
source code in _convert.t_.

The Idl conversion does not work in parallel, it can handle the addition of
extra IO variables (_nwauxio_), and allows to renormalize the data using the
_normt_ and _normvar_ array, in case you directly want to have dimensional
quantities available. An additional script provided is the **doidlcat**
script, which basically concatenates all requested _*.out_ files in a single
file, that can be used with the (few) Idl macro's above. Its intended use for
the 2D advection example would be

    doconvert par/testrho/testrho_vac22 0 20
    doidlcat datamr/testrho/vaclogo 0 20 1

The first line creates the 21 files _datamr/testrho/vaclogo0000.out_ till
_datamr/testrho/vaclogo0020.out_ (assuming you edited the par-file and
indicated the proper _convert_type_ for idl), while the next line then gathers
them all in a single _datamr/testrho/vaclogoall.out_ file, ready for Idl
visualization with VAC-like macro's, like _.r getpict_, _.r plotfunc_ or _.r
animate_. The 3 integer parameters to **doidlcat** indicate the first and last
snapshot number, and a skiprate. If the latter is different from 1, you
include every so many files in the concatenation.

### _convert_type='dx'__

This (heritage) output format is suited for openDX conversion.

### onegrid(mpi), oneblock(B), ...

Extra possibilities to allow creation of a single uniform grid level output.
Please inspect the workings in _convert.t_.

## Endianness {#endianness}

Our _*.dat_ files contain the data in binary format. This is in principle
exchangeable between different platforms, with the possible exception for big
endian or little endian type machines (these differ in the ordering of bytes
for integers and reals). E.g. IBM SP platforms have different endianness than
the VIC3 cluster at K.U.Leuven. In that sense, you may need to do conversion
on the same (type) platform than you used for running the code.

For the binary _.vtu_ files, the endianness can be provided in the xml header
and vtu readers like paraview will then interprete the data correctly. The
default endianness of _.vtu_ files is little endian. To switch to big endian,
simply add the line

    #define BIGENDIAN

to the _definitions.h_ file.

There can be solutions on the machine at hand, using the _assign_ command
(whose syntax you will need to get info on). We would also like to hear if
anyone knows about a way to specify the endianness of the output in
MPI/Fortran itself, independent of the platform. Using Fortran compilers like
gfortran and intel fortran, it is now possible to output .dat data files in
the other endian, e. g. big endian, with a parameter endian_swap=.true. in the
filelist section of the parameter file.
