# Setting parameters

[TOC]

# Introduction {#par_intro}

This document describes the usage of a `.par` parameter file for MPI-AMRVAC.
Note that the default parameter values are set in `amrvacio/amrio.t`, look at
subroutine `readparameters` for details.

# Namelists {#par_namelists}

The parameter file consists of a sequence of namelists, which look like this:

     &listname
      var_a = value

      ! A comment
      var_b = value
      array = value1, value2, ...
      array(1,1) = value

      ! Repeat value 5 times in array
      array = 5*value
      ...
     /

If you do not define a variable the default value is used. The default values
for the (many) parameters are defined in the file `amrio.t`, specifically in the
subroutine `readparameters`.

The Fortran 90 standard for logical variable values is either `T` and `F` or
`.true.` and `.false.`, but some compilers accept only one of them. Text between
namelists is ignored, but it may result in a run time error on some
machines/compilers.

The following namelist examples contain all the possible variables to set,
choices are indicated by `|`. The first choice is the default value. Only the
parameters different from default need to be set. Names that should be replaced
are in capital letters. The `...` indicates optional extra elements for arrays,
or extra words in strings.

## Filelist {#par_filelist}

     &filelist
       <variable definitions, see below>
     /


name | type | default | description
---|---|---|---
filenameini | string | 'unavailable' | If not 'unavailable', resume from snapshot with this base file name
snapshotini | integer | - | Resume from the snapshot with this index. For example, if `filenameini` is "test" and `snapshotini` is 4, the simulation will resume from `test0004.dat`.
filenameout | string | 'data' | Base file name for simulation output, which will be followed by a 4-digit number
filenamelog | string | 'amrvac' | Log file name (without the `.log` extension)
typefilelog | string | 'default' | Use the value 'special' to enable user-defined log output
snapshotnext | integer | 0 | Start index for writing snapshots
slicenext | integer | 0 | Start index for writing slices
firstprocess | logical | F | If true, call `initonegrid_usr` upon restarting
changeglobals | logical | F | If true, call `initglobal` to specify global parameters upon restarting
resetgrid | logical | F | If true, rebuild the AMR grid upon restarting
typeparIO | integer | 1 | 1: Parallel MPI output, 0: master-slave parallel IO, -1: master-slave IO without MPI (no `MPI_FILE_WRITE`, `MPI_FILE_OPEN` etc)
addmpibarrier | logical | F | Enable additional `MPI_BARRIER` calls, useful when debugging on new platforms
convert | logical | F | If true and filenameini and snapshotini are given, convert snapshots to other file formats
convert_type | string | vtuBCCmpi | Which format to use when converting, options are: idl, tecplot, tecplotCC, vtu, vtuCC, vtuB, vtuBCC, dx,  tecplotmpi, tecplotCCmpi, vtumpi,  vtuCCmpi, pvtumpi, pvtuCCmpi, tecline, teclinempi, onegrid
autoconvert | logical | F | If true, already convert to output format during the run
sliceascii | logical | F | If true, enable ASCII output of @ref slices.md
saveprim | logical | F | If true, convert from conservative to primitive variables in output
primnames | string(1:nw) | ' ' | Names of the primitive variables
nwauxio | integer | 0 | !> Number of auxiliary variables that are only included in the output
normvar | real(0:nw) | 1.0 | Conversion factors for length (`normvar(0)`) and the primitive variables (`normvar(1:nw)`)
normt | real | 1.0 | Conversion factor for time unit
`level_io` | integer | - | When doing a convert, generate a uniform grid at this level
`level_io_min` | integer | 1 | Minimum grid level when doing a convert
`level_io_max` | integer | `nlevelshi` | Maximum grid level when doing a convert
nocartesian | logical | F | If true, do not convert the output to a Cartesian coordinate system
uselimiter | logical | F | If true and doing a 1D run, use a limiter to determine corner values
dxfiletype | string | 'lsb' | Data Explorer file endianness ('msb' or 'lsb' for last or most significant bit order)
writew | logical(1:nw) | all true | VTK: Only write variables for which `writew(iw)` is true
writelevel | logical(1:nlevelshi) | all true | VTK: only write these levels
writespshift | real(1:ndim,1:2) | all zero | clip this amount of the domain at the lower and upper side in each dimension

### The log file

By default, the logfile contains one line with a string
corresponding to the `fileheadout' given in methodlist (see below), and one line
with a string that is meant to identify the coordinate names, the conserved
variables (wnames) and other entries, and then follows a sequence of lines
containing numbers: i.e. a single line per requested output time, containing the
integer timestep counter _it_, the time _t_, the time step to be used in the
next time advance _dt_, the domain integrated value of each conserved variable
(nw real numbers, which allows to check perfect conservation across the grid
tree when the boundary conditions imply it), the percentage of the domain
covered by each allowed grid level (_mxnest_ real numbers between 0.0 and 1.0,
with 1.0 indicating 100% coverage: when all _mxnest_ numbers are summed, we get
1.0), and the number of grids per allowed grid level (hence, _mxnest_ integers).
The logfile is by default saved as an ASCII file. When a residual is calculated
(steady- state computations), the value of the residual is stored in the
logfile, just after the time step _dt_ en before all domain integrated values.

The order of saving snapshots, and regridding actions through the subroutine
_resetgridtree_ is fixed: regrid happens after the advance by one timestep,
then regrid, then save the data. This has consequences for the optional
variables beyond _nwflux_.

### Further info on output and postprocessing

The code can be used to postprocess the MPI-AMRVAC .dat files (which are the
only ones to be used for restarts) to some convenient data files for later
visualisation purposes. Such conversion of a single .dat file at a time is to be
done with the same executable (or at least one on a possibly different machine,
but with the same setamrvac settings), on a single processor (i.e. using _mpirun
-np 1 amrvac_). Only selected output types can be converted in parallel, namely
those whose name contains _mpi_ as part of the _convert_type_ string. Currently,
this includes the ASCII versions of _vtumpi_ and _vtuCCmpi_ (corner versus cell
center values), and similarly for tecplot (_tecplotmpi_ or _tecplotCCmpi_). In
addition, _pvtumpi_ and _pvtuCCmpi_ are possible which will result in a _*.vtu_
file for each processor.

In this conversion mode, the idea is to set the filenameini and the
snapshotini entries together with `convert=T`. You can ask the code during
conversion to change from conservative to primitive variable output by setting
`saveprim=T`, and then you should give the corresponding names for the
primitive variables in `primnames`. Just like `wnames` (the latter is
defined in the _methodlist_), this is a single string with space-seperated
labels for the variables. The `primnames` just should list names for the
_nw_ primitive variables, like _wnames_ does for the conservative ones. It is
also possible to perform the conversion step during run of the simulation with
the switch `autoconvert=T`. Naturally, this leads to more computational
overhead and IO, but using the _pvtu(CC)mpi_ filetype, this can be reduced to
a minimum.

For simulations on non-cartesian grids (cylindrical or spherical), there is
the option to output the non-cartesian cell coordinates and vector components,
which then forces you typically to do the conversion to cell center cartesian
grid coordinates and vector variables in your visualization session. By
default (i.e. `nocartesian=F`), the convert module does the conversion from
the orthogonal to locally cartesian coordinates and vector components for you.
You can overrule this default behavior by setting `nocartesian=T`. (note:
for tecplot format, the coordinate labels are then corrected in the converted
file as well).

The only variable that then further matters is `convert_type`. Selecting
_'idl'_ will generate a corresponding 'datamr/FILEINIBASExxxx.out' file, which
is stored in binary and can be handled with the Idl macros. For type
'tecplot', a corresponding `datamr/FILEINIBASExxxx.plt` file will be
generated, which is an ASCII file that stores the cell corner locations and
corner values for the conserved variables, to be handled with Tecplot. The
'onegrid' conversion type is just useful in 1D AMR runs, to generate a
single block file (extension '.blk'). Also particular to 1D data, and for
TecPlot purposes alone, is the 'tecline' option. This can also be done in
parallel mode, where it is called 'teclinempi'.

The type 'dx', will generate a Data Explorer file 'datamr/FILEINIBASExxxx.dx',
which can be used with the free [DX visualization software](www.opendx.org).
When `convert_type='dx'`, there is the additional `dxfiletype` variable. The dx
filetype stores in binary format, and stores cell center coordinates and values.
The binary format for dx files can differ from machine to machine, but will be
one of 'lsb' or 'msb' (for last or most significant bit order, respectively).

For visualization using [Paraview](www.paraview.org), the option to convert to
`convert_type='vtu'` can be used. For both _vtu_ and _tecplot_ formats, there
are also corresponding _vtuCC_ and _tecplotCC_ options, which store the data
with the actually computed cell-centered values. For the _vtu_ and _tecplot_
formats on the other hand, the code tries to already do some interpolation from
cell center to cell corner variables for you, but this may introduce some
artificial effects in non-cartesian geometries. The _vtuB_ and _vtuBCC_ do the
same as _vtu(CC)_ but save the data in binary format.

It is even possible to temporarily add additionally computed auxiliary
variables that are instantaneously computable from the data in the
`datamr/FILEOUTBASExxxx.dat` file to the converted snapshot. You should then
provide the number of such extra variables in `nwauxio` (see also
[this page](mpiamrvac_nw.md)), and a corresponding
definition for how to compute them from the available _nw_ variables in the
subroutine _specialvar_output_ whose default interface is provided in the
_amrvacnul.speciallog.t_ module. You can there compute variables that are not
in your simulation or data file, and store them in the extra slots
_nw+1:nw+nwauxio_ of the _w_ variable. For consistency, you should also then
add a meaningfull name to the strings that we use for identifying variables,
namely _primnames, wnames_. This has to be done in the subroutine
_specialvarnames_output_.

The output values are normally given in code units, i.e. in the dimensionless
values used throughout the computation (in the initial condition, we always
adhere to the good practice of choosing an appropriate unit of length, time,
mass and expressing everything in dimensionless fashion). One can, in the
convert stage only, ask to multiply the values by their respective dimensional
unit value. `normt` should then be the unit for time, while the array
`normvar` combines the unit of length (to be stored in _normvar(0)_) with
corresponding units for all primitive variables in _normvar(1:nw)_. The
corresponding values for conservative entries are computed in _convert.t_ when
_saveprim=F_. Note that these latter are not all independent and must be set
correctly by the user. In any case, they are not in the restart files (the
ones with _.dat_), and just used at conversion stage. See for details of their
use the _convert.t_ module.

There is a `uselimiter` logical variable, which is false by default, but can
be used when having a 1D Cartesian grid computation, where it then influences
the way in which the convert step computes corner values from cell center
values. Normally, it would just take the arithmetic average as in
_0.5(w_L+w_R)_ where _w_L_ is the cell center value at left, and _w_R_ is at
right. Activating uselimiter will compute _0.5(w_LC+w_RC)_ instead, where the
left and right edge values _w_LC, w_RC_ are computed by limited reconstruction
first.

Note that different formats for postprocess data conversion can be added in
the `convert.t` subroutine. See [convert](convert.md) for
details.

The _VTK_-based formats allow for saving only a part of the _nw_ variables, by
setting the logical array `writew`. The same is true for selecting levels by
using the array `writelevel`. Finally, you can clip away part of the domain,
for output in a selected region. This is done by filling the array
`writespshift`. That array should use (consistent) double precision values
(between 0 and 1) that specify the percentage of the total domain to be
clipped away from the domain boundary at the minimum side, and from the
maximum side, respectively, and this for all dimensions. The array has thus a
dimensionality _(1:ndim,1:2)_, with the first entry specifying the dimension,
and the second whether you clip for minimum (1) and maximum (2) sides,
respectively. Note that in the end, only the grids that are fully contained
within the clipped region are then converted and stored in the output.

The switches `level_io`, `level_io_min` and `level_io_max` are there to
restrict the AMR levels of the output file at the convert stage. These
switches do not work with autoconvert. E.g. setting _level_io=3_ will
coarsen/refine the data to level 3 everywhere, resulting in a uniform grid.
_level_io_min_ limits the minimum level for output by refining levels below
level_io_min until level_io_min is reached. Correspondingly, _level_io_max_
limits the maximum level of the output file. This can be useful to visualize
large datasets.

## Savelist {#par_savelist}

Example:

    &savelist
        itsave(1,1)=0
        itsave(1,2)=0
        ditsave(1)=1
        dtsave(2)=0.05d0
        dtsave(4)=0.1
    /


name | type | default | description
---|---|---|---
`ditsave(FILEINDEX)` | integer | `biginteger` | Repeatedly save output when `ditsave(FILEINDEX)` time steps have passed
`dtsave(FILEINDEX)` | double | `bigdouble` | Repeatedly save output when `dtsave(FILEINDEX)` simulation time has passed
`itsave(SAVEINDEX,FILEINDEX)` | integer | 1 | Save on these time steps
`tsave(SAVEINDEX,FILEINDEX)` | double | `bigdouble` | Save on these times
nslices | integer | 0 | Number of slices
`slicedir(INTEGER)` | integer | - | Slice direction, see @ref slices.md
`slicecoord(INTEGER)` | double | - | Slice coordinate, see @ref slices.md
`collapse(INTEGER)` | logical | F | See @ref collapsed.md
collapseLevel | integer | 1 | See @ref collapsed.md

Here FILEINDEX has the following meaning:

index | meaning
---|---
1 | Log output
2 | Normal output
3 | Slice output, see @ref slices.md
4 | Collapsed output, see @ref collapsed.md
5 | Call analysis subroutine, see @ref analysis.md

The times can be given in timesteps or physical time. Typical examples:
`ditsave=1,10` saves results into the log file at every timestep, and will
generate a .dat output file every 10-th step (however, the number at the end
of the `datamr/FILEOUTBASExxxx.dat` file will increase by one from
firstsnapshot, as explained above). `dtsave=0.1,12.5` saves into the log
file at times 0.1,0.2,... and will generate a .dat output file at time
12.5,25,37.5,... , assuming that we start at t=0. `ditsave(1)=10
tsave(1,2)=5.2 tsave(2,2)=7.` will save info into the log file every 10-th
timestep and snapshots at t=5.2 and 7. Actually, the first timestep when
physical time is greater than 5.2 (or 7.0) is saved. Mixing itsave and tsave
is possible, mixing dtsave (or ditsave) with tsave (or itsave) for the same
file should be done with care, since dtsave (and ditsave) will be offset by
any intervening tsave (and itsave). However, one may want to save snapshots
more frequently at the beginning of the simulation. E.g. `tsave(1,2)=0.1
tsave(2,2)=0.25 tsave(3,2)=0.5 dtsave(2)=0.5` could be used to save snapshots
at times 0.1, 0.25, 0.5, 1, 1.5, ... etc.

If no save condition is given for a file you get a warning, but _the final
output is always saved_ after the stop condition has been fulfilled. If
`itsave(1,2)=0` is set, the initial state is saved before advancing.

## Stoplist {#par_stoplist}

         &stoplist
    	itmax =INTEGER
    	tmax =DOUBLE
    	tmaxexact =F | T
    	dtmin =DOUBLE
    	it =INTEGER
    	t =DOUBLE
    	treset =F | T
    	itreset =F | T
    	residmin =DOUBLE
    	residmax =DOUBLE
     typeresid = 'relative' | 'absolute'
    /

You may use an upper limit `itmax` for the number of timesteps and/or the
physical time, `tmax`. If `tmaxexact=T` is set, the last time step will be
reduced so that the final time 't' is exactly 'tmax'.

Numerical or physical instabilities may produce huge changes or very small
time steps depending on the way `dt` is determined. These breakdowns can be
controlled by either setting a lower limit `dtmin` for the physical time
step, which is useful when `dt` is determined from the `courantpar`
parameter. If AMRVAC stops due to `dt &lt; dtmin`, a warning message is
printed.

You have to specify at least one of `tmax, itmax`. AMRVAC stops execution
when any of the limits are exceeded. The initial time value `t` and integer
time step counter `it` values can be specified here. However, when a restart
is performed from a previous .dat file, the values in that file will be used
unless you enforce their reset to the values specified here by activating the
logicals `treset=T`, or `itreset=T`.

In case the code is used for computing towards a steady-state, it is useful to
quantify the evolution of the residual as a function of time. The residual
will be computed taking account of all _nwflux_ variables (see also
[this page](mpiamrvac_nw.md)), over all grids. You can
tell the code to stop computing when a preset minimal value for this residual
is reached, by specifying `residmin=1.0d-6` (for example, a rather stringent
value when all variables are taken into account). Similarly, setting
`residmax=1.0d8` will force the code to stop when the residual change
becomes excessively large (indicating that you will not reach a steady state
for the problem at hand then, without explaining you why). The residual is
computed as specified in the subroutine _getresidual_ in the _setdt.t_ module,
and you may want to inspect what it actually computes there (and perhaps
modify it for your purposes), and see the distinction between
_typeresid='relative'_ or _typeresid='absolute'_. When either residmin or
residmax is specified here, the value of the residual will be added to the
logfile.

## Methodlist {#par_methodlist}

        &methodlist

    wnames=	'STRING'
    fileheadout= 'STRING'

    typeadvance='twostep' | 'onestep' | 'threestep' | 'rk4' | 'fourstep' | 'ssprk43' | 'ssprk54'
    typefull1=nlevelshi strings from: 'tvdlf','hll','hllc','hllcd','tvdmu','tvd','cd','fd','hll1','hllc1','hllcd1','tvd1','tvdlf1','tvdmu1','source','nul'
    typepred1=nlevelshi strings from: 'default','hancock','tvdlf','hll','hllc','tvdmu','cd','fd','nul'
    typelow1= nlevelshi strings from: 'default' | 'tvdlf1' | 'hll1' | 'hllc1' | 'hllcd1' | 'tvdmu1' | 'tvd1' | 'cd' | 'fd'

    typelimiter1= nlevelshi strings from: 'minmod' | 'woodward' | 'superbee' | 'vanleer' | 'albada' | 'ppm' | 'mcbeta' | 'koren' | 'cada' | 'cada3' | 'mp5'
    typegradlimiter1= nlevelshi strings from: 'minmod' | 'woodward' | 'superbee' | 'vanleer' | 'albada' | 'ppm' | 'mcbeta' | 'koren' | 'cada' | 'cada3'
    typelimited='original' | 'previous' | 'predictor'
    useprimitive= T | F
    loglimit= nw logicals, all false by default
    flatsh = F | T
    flatcd = F | T
    flatppm= T | F
    mcbeta= DOUBLE

    typeentropy= 'nul','powell','harten','ratio','yee'
    entropycoef= DOUBLE, DOUBLE, DOUBLE, ....

    typetvd= 'roe' | 'yee' | 'harten' | 'sweby'
    typeaverage='default' | 'roe' | 'arithmetic'

    typetvdlf= 'cmaxmean' | 'other'
    tvdlfeps = DOUBLE
    BnormLF = T | F

    flathllc= F | T
    nxdiffusehllc = INTEGER

    typeaxial= 'slab' | 'cylindrical' | 'spherical'
    typespherical= 1 | 0
    usecovariant= F | T

    ssplitdivb= F | T
    ssplitdust= F | T
    ssplitresis= F | T
    ssplituser= F | T
    typesourcesplit= 'sfs' | 'sf' | 'ssfss' | 'ssf'
    dimsplit= F | T
    typedimsplit= 'default' | 'xyyx'| 'xy'

    smallrho= DOUBLE
    smallp= DOUBLE
    fixsmall= T | F
    strictsmall= T | F
    strictgetaux= F | T
    nflatgetaux=1
    fixprocess= F | T

    typedivbfix= 'powel' | 'janhunen' | 'linde' | 'glm1' | 'glm2' | 'glm3'
    divbwave= T | F
    divbdiff= DOUBLE
    typedivbdiff= 'all' | 'ind'
    B0field= F | T
    Bdip= DOUBLE
    Bquad= DOUBLE
    Boct= DOUBLE
    Busr= DOUBLE
    compactres= F | T
    typegrad = 'central' | 'limited'
    typediv = 'central' | 'limited'

    useprimitiveRel= T | F
    maxitnr= INTEGER
    tolernr= DOUBLE
    absaccnr= DOUBLE
    dmaxvel= DOUBLE
    typepoly= 'meliani' | 'bergmans' | 'original' | 'gammie'
    strictnr= T | F
    strictzero= T | F

    ncool= INTEGER
    cmulti = INTEGER
    coolmethod= ' '
    coolcurve= ' '
    Tfix= F | T

    ptmass= DOUBLE
    x1ptms= DOUBLE
    x2ptms= DOUBLE
    x3ptms= DOUBLE

    dustmethod= 'Kwok','sticking','linear','none'
    dustzero = T|F
    dustspecies = 'graphite','silicate'
    dusttemp = 'constant','ism','stellar'
    smallrhod = DOUBLE'

    /

### wnames, fileheadout {#par_wnames}

`wnames` is a string of (conserved) variable names, which is only stored for use in the input and output data files. The `wnames` string of variable names is used in the default log file header, and should be e.g. 'rho m1 m2 e b1 b2' for MHD with 2 vector components. These labels are the ones passed on when doing the conversion step for tecplot, vtu formats.

`fileheadout` is the header line in the output file, and is an identifyer used in the (converted) data files. Both wnames and fileheadout are thus only relevant within amrio.t and convert.t, and only have effect for the later data visualization (e.g. when fileheadout='test_mhd22' is specified, the Idl macros will deduce from the name that it is for a 2D MHD run).

### typeadvance, typefull1, typepred1, typelow1 {#par_typeadvance}

The `typeadvance` variable determines the time integration procedure. The
default procedure is a second order predictor-corrector type 'twostep' scheme
(suitable for TVDLF, TVD-MUSCL schemes), and a simple 'onestep' algorithm for
the temporally second order TVD method, or the first order TVDLF1, TVDMU1,
TVD1 schemes. `It is not possible to mix different step size methods across
the AMR grid levels.` The temporally first order but spatially second order
TVD1 algorithm is best suited for steady state calculations as a 'onestep'
scheme. The TVDLF and TVD-MUSCL schemes can be forced to be first order, and
linear in the time step, which is good for getting a steady state, by setting
`typeadvance='onestep'`.

There is also a fourth order Runge-Kutta type method, when
`typeadvance='fourstep'`. It can be used with _dimsplit=.true._ and
_typelimited='original'_. These higher order time integration methods can be
most useful in conjunction with higher order spatial discretizations like a
fourth order central difference scheme (currently not implemented). See also
[discretization](discretization.md).

The array `typefull1` defines a spatial discretization
[method](methods.md) used for the time integration per activated grid level
(and on each level, all variables use the same discretization). In total,
_nlevelshi_ methods must be specified, by default _nlevelshi=8_ and these are
then all set by _typefull1=8*'tvdlf'_. Different discretizations can be mixed
across the _mxnest_ activated grid levels (but the same stepping scheme must
apply for all of the schemes).

Setting for a certain level the typefull1 to 'nul' implies doing no advance at
all, and 'source' merely adds sources. These latter two values must be used
with care, obviously, and are only useful for testing source terms or to save
computations when fluxes are known to be zero.

The `typepred1` array is only used when `typeadvance='twostep'` and
specifies the predictor step discretization, again per level (so _nlevelshi_
strings must be set). By default, it contains _typepred1=8*'default'_ (default
value _nlevelshi=8_), and it then deduces e.g. that 'cd' is predictor for
'cd', 'hancock' is predictor for both 'tvdlf' and 'tvdmu'. Check its default
behavior in the _amrio.t_ module. Thus `typepred1` need not be defined in
most cases, however `typefull1` should always be defined if methods other
than 'tvdlf' are to be used.

The `typelow1` is only used when the AMR strategy is based on Richardson
extrapolation. It then specifies the first order scheme used in this process,
to be done per level. By default, `typelow1=nlevelshi*'default'` will imply
the use of `'tvdlf1'` for all methods, except for 'cd' where it is 'cd'. It
is the same low order scheme that is then used in the Richardson process
across the entire grid hierarchy, and is decided on the basis of the
typefull1(l) for level l in amrio.t.

### Limiter type {#par_typelimiter}

For the TVDLF and TVD-MUSCL methods different limiter functions can be defined
for the limited linear reconstructions from cell-center to cell-edge
variables, and for the TVD method, for the characteristic variables. See the
`src/amrvacpar.t` file for the order of the characteristic variables. The
default limiter is the most diffusive `typelimiter1=nlevelshi*'minmod'`
limiter (minmod for all levels), but one can also use
`typelimiter1=nlevelshi*'woodward'`, or use different limiters per level.

The `typegradlimiter1` is the selection of a limiter to be used in computing
gradients (or divergence of vector) when the typegrad=limited (or
typediv=limited) is selected. It is thus only used in the gradientS
(divvectorS) subroutines in geometry.t (and has effect for the MHD modules).

The `typelimited` variable tells the TVD type methods what should be used as
a basis for the limiting. By default, the `original` value is used in 1D and
for dimensional splitting, while for dimensionally unsplit multidimensional
case (dimsplit=F), TVDLF and TVD-MUSCL uses the `previous` value from
`wold` for limiting.

The `useprimitive` variable decides whether TVDLF and TVDMUSCL schemes
should limit the slopes of the _primitive_ or _conservative_ variables. The
default behaviour is limiting the primitive variables. For the onestep TVD
scheme the 'useprimitive' parameter determines how the pressure jump is
calculated in the approximate Riemann solver. For the false value the jump is
calculated from the jumps in the conservative variables following strictly the
TVD algorithm, while for the default true value the jump is approximated by
the pressure difference, which is simpler and slightly faster. This subtle
change does not seem to influence the results obtained by the TVD scheme
significantly.

The use of `useprimitive=T` can be combined with the selection of
logarithmic transformation on (positive) variables. I.e., when e.g. having a
gravitational stratification, one might benefit from performing linear
reconstruction on the primitive variables log10(rho) and/or log10(p). This can
be done by setting the corresponding _loglimit(iw)=T_ with _iw_ the label of
the corresponding component in the _w_ array (for density, this is thus
_iw=1_).

When using PPM as a limiter, minor differences can be obtained using the
switches flatppm, flatcd, flatsh. The last two are meant to minimize potential
ripples around contact discontuinities (flatcd) or shocks (flatsh), but one
should first try without these flattenings (default behavior). PPM is actually
only used in a quadratic reconstruction from center to edge, requires the use
of a larger stencil (dixB=4), and can be used either in the methods (by
setting typelimiter1) or in the gradientS/divvectorS routines (when typegrad
or typediv is limited, and typegradlimiter1 is ppm). The latter is encoded in
geometry.t.

### Typeentropy, entropycoef {#par_typeentropy}

For Riemann solver based methods, such as TVD and TVD-MUSCL (but not TVDLF),
an entropyfix may be applied to avoid unphysical solutions. The fix is applied
to the characteristic variables, their order is defined in
`src/amrvacpar.t`. The default entropy fix is `'nul'`, i.e. no entropy
fix. When an expansion shock is formed, the entropy fix should be applied to
the non-degenerate characteristic waves, i.e. waves that can form shocks
(sound waves, fast and slow magnetosonic waves). The most powerful entropy fix
is called 'powell'. In practice, one may apply an entropy fix to all
characteristic waves, usually the slight extra diffusion makes the schemes
more robust. For Yee's entropyfix the minimum characteristic speed (normalized
by dt/dx) can be set for each characteristic wave using the `entropycoef`
array.

### Different TVD varians {#par_tvdvariants}

Both `tvd` and `tvdlf` have a few variants, these can be set in the
strings `typetvd` and `typetvdlf`, with defaults 'roe' and 'cmaxmean',
respectively. The default `typetvd='roe'` is the fastest of the four upwind
types. For the TVDLF, the 'cmaxmean' merely means whether the maximal physical
propagation speed is determined as the maximum speed for the mean state based
on averaging left centered and right centered states, or by taking the maximum
physical speed among both deduced from left and right centered states (the
latter is the Local Lax-Friedrichs variant). In the TVDLF flux, the diffuse
flux part has a coefficient `tvdlfeps` which is 1 by default. For steady-
state computations, one may gain in shock sharpness by reducing this factor to
a positive value smaller than 1.

_Only for the adiabatic hydro module_, the option to select an arithmetic, or
a roe average is available for use in the roe solver. This is set by the
`typeaverage`.

Just like in the TVDLF method, the slightly more involved HLL method has a
diffuse flux part with coefficient `tvdlfeps` which is 1 by default. For
steady-state computations, one may gain in shock sharpness by reducing this
factor to a positive value smaller than 1.

For calculations involving magnetic fields (variable b0_ is positive),
`BnormLF=T` actually uses the Lax-Friedrichs flux expression for the normal
magnetic field component in HLL and HLCC methods. This improves robustness.

When using the HLLC scheme variants, for HD, MHD, SRHD and SRMHD there is an
optional additional flattening in case the characteristic speed at the contact
is near zero. This is activated by setting `flathllc=T` (its default is
false). One can also solve some potential noise problems in the HLLC by
swithcing to the HLLCD variant, a kind of mix between HLLC and TVDLF. The
TVDLF is then used in a user-controlled region around a point where there is a
sign change in flux, whose width is set by `nxdiffusehllc` (an integer which
is 0 by default).

### Typeaxial, typespherical, usecovariant {#par_typeaxial}

`typeaxial` defines the type of curvilinear grid. For cylindrical coordinate systems, the _-phi=_ and _-z=_ flags have a meaning and should be used at _$AMRVAC_DIR/setup.pl_, to denote the order of the second and third coordinate. Together, they control the addition of the geometrical source terms as implemented in the various physics modules in the subroutine _addgeometry_. By default, `typeaxial='slab'` and Cartesian coordinates are used (with translational symmetry for ignorable directions).

In 1D where _-d=11_ we have 1 coordinate with one vector component, and for
cylindrical or spherical cases, the coordinate is the radial cylindrical or
spherical distance, respectively, with corresponding radial vector component.
For _-d=12_ and `typeaxial='cylindrical'`, the second vector component can
be either the phi or the z component, depending on the -phi and -z flags. For
_-d=12_ and `typeaxial='spherical'`, the second vector component is always
the theta component. Similarly for _-d=13_, only the cylindrical case has a
choice for the order of the vector components. Not all combinations make
sense, though.

In 2D 'slab' means translational symmetry when performing 2.5D simulations
(i.e. _$AMRVAC_DIR/setup.pl -d=23_).

For 2D and 'cylindrical' (which can be 2D or 2.5D) the grid and the symmetry
depend on the settings for the -phi and -z flags. When -d=22 -z=2, a cartesian
grid is used in a poloidal plane, with axial symmetry for the r- and z- vector
components. The same is true for -d=23 -z=2 -phi=3, when all three vector
components are then depending on (r,z) coordinates only. The vector components
then appear in the r,z,phi order. One can use 2.5D on a circular grid also
with translational symmetry in the third, i.e. axial, direction by the use of
_$AMRVAC_DIR/setup.pl -d=23 -phi=2 -z=3_. The vector components then appear in
the r,phi,z order.

For 2D and 'spherical', the coordinates denote radial and polar angle in a
poloidal plane, and similarly for 2.5D in combination with spherical where the
third vector component is then the phi- component, depending on (r,theta)
alone. This means that 2D and 'spherical' implies the use of a polar grid in a
poloidal cross-section (i.e. one containing the symmetry axis) and axial
symmetry for 2.5D runs.

In 3D the choice of curvilinear grid is Cartesian for 'slab', and the usual
Cylindrical and Spherical coordinate systems when setting one of those latter
values. Note that vector components are to be interpreted in the corresponding
coordinate system!

Please read [axial](axial.md) before you try to do
simulations in non-slab symmetry.

In case you select `typeaxial='spherical'`, the geometrical info is filled
in a slightly different way depending on the integer _typespherical_. See the
details in _geometry.t_ module.

The _usecovariant_ option is as yet inactive, but is meant to prepare for
general relativistic modules, and/or an alternative means to handle non-
cartesian geometries.

### Dimensional splitting {#par_dimsplit}

The sources, if any, can be added in a split or unsplit way according to the
logical variables `ssplitdivb`, `ssplitdust`, `ssplitresis`, and
`ssplituser` which correspond to divb source to maintain divergence free of
magnetic field, dust effect, resistivity, and other sources added by user,
respectively. Their default values are false meaning these sources are added
in a unplit way by default. The split sources are added according to
`typesourcesplit`. The meaning of the different options for
`typesourcesplit` is described in
[discretization](@ref disc-splitting). Under default
settings, we use unsplit sources only, and if one reverts to split sources,
`typesourcesplit='sfs'`.

In multidimensional calculations dimensional splitting can be used by setting
`dimsplit=T`, with an alternating order of the sweeps
`typedimsplit='xyyx'` by default. For AMRVAC simulations, it is best to use
`dimsplit=F`, the default value, but the TVD method needs a dimensionally
split strategy. The limitations on using dimensionally unsplit methods are
described in [methods](methods.md).

### Positivity fixes {#par_positivityfix}

Negative pressure or density caused by the numerical approximations can make the
code crash. For HD, MHD and all SRHD and SRMHD variants, this can be monitored
or even cured by the `smallvalues` subroutines in
`src/amrvacphys.correctaux[mod].t` modules, where `[mod]` denotes the physics
(hd, mhd, srhd, srmhd). This monitoring is active whenever `fixsmall=T`, its
default setting (hence, you can avoid all checks, but also all cures, by setting
`fixsmall=F`). The control parameters `smallrho, smallp` play a role here:
they can be set to small positive values, while their default is
`smallrho=-1.0d0` and `smallp=-1.0d0`, i.e. no replacements at all. They
take effect for HD, MHD, SRHD and SRMHD equations when set to a small value,
e.g. `smallp=1.0d-12`. Actually, they in turn determine _minp, minrho, smalle_
for HD and MHD modules (as set in the `initglobaldata` subroutine in the
physics module) and _minp, minrho, smalltau, smallxi_ for the relativistic
variants. These latter quantities appear in the
`src/amrvacphys.correctaux`.t` modules.

The actual treatment involves the _strictsmall_ parameter: Its default value
(T) causes a full stop when the smallvalues subroutine in the physics modules
would normally correct small densities or energies by some artificial vacuum
prescription. This corrective prescription is thus turned off by default. In
this way, you can use it for debugging purposes, to spot from where the actual
negative and unphysical value gets introduced. If it is somehow unavoidable in
your simulations, then you may rerun with a recovery process truned on as
follows. When setting _strictsmall=F_, two kinds of recovery procedures can be
selected, controlled by the logical _strictgetaux_. When _strictgetaux=T_, the
parameters smallp, smallrho (and derived values minrho, minp, smalle,
smalltau, smallxi) are used in an ad-hoc prescription for dealing with
`vacuum', as encoded in `src/amrvacphys.correctaux`.t`. In case you select
_strictgetaux=F_, the subroutine _correctaux_ is used instead, which uses some
kind of averaging from a user-controlled environment about the faulty cells.
The width of this environment is set by the integer _nflatgetaux_.

Setting `fixprocess=T` results in a call to the process subroutine before
the writing of a snapshot, and following the determination of the timestep
constraint by means of CFL and other restrictions. It then interfaces to the
_process_grid_usr_ subroutine whose default interface is found in
_amrvacnul.speciallog.t_. You can use this routine for doing computations of
non-local auxiliary variables (like the divergence of some vector fields,
etc), then using these in turn to do particle acceleration treatments (the
latter to be implemented in the process subroutine, using MPI!), etc.

### Magnetic field divergence fixes {#par_divbfix}

Depending on `typedivbfix`, sources proportionate to the numerical monopole
errors are added to momemtum, energy, and induction equation (the 'powel'
type), or to the induction equation alone (the 'janhunen' type). The latter
type can also be used for SRMHD cases. The `divbwave` switch is effective
for the Riemann type solvers for multi-D MHD only. The default true value
corresponds to Powell's _divergence wave_ which stabilizes the Riemann solver.
Naturally, if Powell's source terms are to be added and in practice it is best
to keep the sources out of the Richardson type error estimator (when used), so
that we advocate `ssplitdivb=T`.

Another source term strategy for monopole error control is to do parabolic
cleaning, i.e. add source terms which diffuse the local error at the maximal
rate still compliant with the CFL limit on the time step. This is activated
when `divbdiff` is set to a positive number, which should be less than 2,
and again comes in various flavors depending on which equations receive source
terms. The choice where only the induction equation gets modified, i.e.
`typedivbdiff='ind'` can again be used for both MHD and SRMHD, and is
advocated.

For MHD, we implemented the possibility to use a splitting strategy following
Tanaka, where a time-invariant, potential background magnetic field is handled
exactly, so that one solves for perturbed magnetic field components instead.
This field is taken into account when `B0field=T`, and the magnitude of this
field is controlled using the variables `Bdip, Bquad, Boct, Busr`. The first
three are pre-implemented formulae for a dipole, quadrupole and octupole field
in spherical coordinates only (the parameters set the strength of the dipole,
quadrupole and octupole field). This is coded up in the module _set_B0.t_.
This same module calls in addition the _specialset_B0_ subroutine when
`Busr` is non-zero, where it then should be used to quantify an additional
potential, time-independent field. This latter can be used for cartesian or
cylindrical coordinates as well. This splitting strategy can be extended to
linear force-free background field with non-zero current by adding a source
term, namely, perturbed magnetic field cross velocity then dot background
current, in the right hand side of the conservative energy equation as a
user's unsplit special source .

Resistive source terms for MHD can use a compact non-conservative formulation
of resistive source terms, by setting compactres=T. The default
`compactres=F` setting is normally preferred.

The `typegrad` can be selected to switch from simple centered differencing
on the cell center values, to limited reconstruction followed by differencing
when computing gradients. They call either of _gradient_ ('central') or
_gradientS_ ('limited') subroutines that are themselves found in the
_geometry.t_ module. Similarly, a switch for the divergene of a vector is the
`typediv` switch. These are as yet only used in the MHD modules (classical
and relativistic). When the 'limited' variant is used, one must set the
corresponding typegradlimiter1 array to select a limiter (per level).

### Using primite variables for reconstruction {#par_useprim}

For the SRHD and SRMHD modules only. The `useprimitiveRel=T` will ensure
that in combination with useprimitive, limited linear reconstruction is done
on the spatial four-velocity instead of the velocity. It is the default value.
The special relativistic modules involve a Newton-Raphson procedure for
switching from conservative to primitive, and the maximum number of NR
iterates is set to a default value `maxitnr=100`. This newton raphson has a
tolerance parameter and absolute accuracy parameter, by default set to
_tolernr=1.0d-13_ and _absaccnr=1.0d-13_. These can be changed if needed. The
logical _strictnr_ (T by default) will cause the code to stop when there is an
error related to this NR procedure. See the detailed implementation in the
_amrvacphys.t.sr(m)hd(eos)_ modules.

For the relativistic MHD modules, an additional parameter is the
_strictzero=.true._ parameter: it sets the limit for zero flow conditions as
measured by S.S (or even B.B) in the NR procedure to zero (T) or
smalldouble*smalldouble (F). Further, the `dmaxvel=1.0d-8` default value is
used in the process to define a maximal velocity allowed, namely 1-dmaxvel
(where velocities are always below 1, which is the speed of light value). For
SRMHD, the `typepoly` determines which formulation is used for the quartic
polynomial whose zeros determine the forward and backward slow and fast signal
speeds.

### ncool, cmulti, coolmethod, coolcurve, Tfix

These are only used in combination with a cooling module for HD and MHD, as
developed by AJ van Marle, described in
[mpiamrvac_radcool.md](mpiamrvac_radcool.md).

### ptmass, x1ptms,x2ptms,x3ptms

These are only used in combination with an optional point gravity module for
HD and MHD, as developed by AJ van Marle, described in
[mpiamrvac_pointgrav.md](mpiamrvac_pointgrav.md).

### dustmethod, dustzero,dustspecies,dusttemp,smallrhod

These are only used when one or more dustspecies is used for HD.

## Boundlist {#par_boundlist}

         &boundlist
     dixB= INTEGER
    	typeB= 'cont','symm','asymm','periodic','special','noinflow','limitinflow'
     ratebdflux = DOUBLE
     internalboundary = F | T
    	typeghostfill= 'linear' | 'copy' | 'unlimit'
    	typegridfill= 'linear' | 'other'
    /

The boundary types have to be defined for each conserved variable at each
physical edge of the grid, i.e. for 2D hydrodynamics they are in the order:
rho,m1,m2,e at the left boundary; rho,m1,m2,e at the right; rho,m1,m2,e at the
bottom; finally rho,m1,m2,e at the top boundary. In general, the order is
xmin, xmax, ymin, ymax, zmin, and zmax.

The default number of ghost cell layers used to surround the grid (and in fact
each grid at each level and location) is set by default to `dixB=2`. If
needed, this value can be increased.

The default boundary type is `cont` for all variables and edges, it means
that the gradient (of the conservative variables) is kept zero by copying the
variable values from the edge of the mesh into the ghost cells.

Other predefined types are the `symm` and `asymm` types, which are mostly
used for reflective boundaries, or at symmetry axes of the domain (the polar
or equatorial axis, e.g.). One then typically makes the momentum orthogonal to
the given boundary antisymmetric (`asymm`), the rest of the variables
`symm`. These boundary types can also be used to represent a perfectly
conducting wall (the orthogonal component of the magnetic field should be
antisymmetric, the transverse component symmetric) or the physical symmetry of
the physical problem.

The case of periodic boundaries can be handled with setting 'periodic' for all
variables at both boundaries that make up a periodic pair. Hence triple
periodic in 3D MHD where 8 variables are at play means setting
_typeB=8*'periodic',8*'periodic',8*'periodic',8*'periodic',8*'periodic',8*'periodic_.
For 3D cylindrical and spherical grid computations, the singular polar axis is
trivially handled using a so-called pi-periodic boundary treatment, where
periodicity across the pole comes from the grid cell diagonally across the
pole, i.e. displaced over pi instead of 2 pi. These are automatically
recognized from the typeaxial setting, and the corresponding range in angle
phi must span 2 pi for cylindrical, and theta must then start at zero (to
include the north pole) and/or end at pi (for the south pole) for spherical
grids. The user just needs to set the typeB as if the singular axis is a
symmetry boundary (using symm and asymm combinations).

The possibility exists to put a boundary condition mimicking zero or reduced
inflow across the computational boundary, by selecting _typeB='noinflow'_ or
_typeB='limitinflow'_ for the momentum vector components of your particular
application. This is in principle only relevant for the momentum component
locally perpendicular to the boundary (for others a continuous extrapolation
is done). The _noinflow, limitinflow_ extrapolates values that are outwardly
moving continuously, while clipping all values that are inwardly advecting
momentum to zero (noinflow) or to a user-controlled fraction of the inward
momentum (limitinflow). The latter fraction is set by _ratebdflux_ which is 1
by default, and should be set to a value between zero and 1 accordingly.

The `special` type is to be used for setting fixed values, or any time
dependent or other more complicated boundary conditions, and results in a call
to the `specialbound_usr` subroutine which has to be provided by the user in
the [AMRVACUSR module](@ref #amrvacusr_specialbound). The variables with
`special` boundary type are updated last within a given boundary region,
thus the subroutine may use the updated values of the other variables. The
order of the variables is fixed by the equation module chosen, i.e. _rho m1 m2
m3 e b1 b2 b3_ for 3D MHD, but by setting all typeB entries for a certain
boundary region to special, one is of course entirely free to fill the
boundary info in a user-defined manner.

Internal boundaries can be used to overwrite the domain variables with
specified values. This is activated with the switch `internalboundary=T`.
Internally, these are assigned before the ghost-cells and external boundaries
are applied (in subroutine routine get_bc). The user can provide conditions on
the conserved variables depending on location or time in the subroutine bc_int
which is defaulted in amrvacnul.specialbound.t.

The `typeghostfill='linear'` implies the use of limited linear
reconstructions in the filling of ghost cells for internal boundaries that
exist due to the AMR hierarchy. A first order 'copy' can be used as well, or
an unlimited linear reconstruction by setting it to 'unlimit'. To retain
second order accuracy, at least the default 'linear' type is needed.

The `typegridfill='linear'` implies the use of limited linear
reconstructions when filling newly triggered, finer grids from previous
coarser grid values. Setting it different from this default string will imply
mere first order copying for finer level grids (and is thus not advised when
second order is desired).

## Amrlist {#par_amrlist}

         &amrlist
    	mxnest= INTEGER
    	nxlone1= INTEGER
    	nxlone2= INTEGER
    	nxlone3= INTEGER
    	dxlone1= DOUBLE
    	dxlone2= DOUBLE
    	dxlone3= DOUBLE
     xprobmin1= DOUBLE
     xprobmax1= DOUBLE
     xprobmin2= DOUBLE
     xprobmax2= DOUBLE
     xprobmin3= DOUBLE
     xprobmax3= DOUBLE
     errorestimate= INTEGER
    	nbufferx1= INTEGER
    	nbufferx2= INTEGER
    	nbufferx3= INTEGER
     skipfinestep= F | T
     amr_wavefilter= nlevelshi DOUBLE values
     tol= nlevelshi DOUBLE values
     tolratio= nlevelshi DOUBLE values
     flags= INTEGER array, with at most nw+1 entries
     wflags= DOUBLE array, with at most nw values that must sum up to 1.0d0
     logflag= nw logical values, all F by default
     iprob= INTEGER
     prolongprimitive= F | T
     coarsenprimitive= F | T
     restrictprimitive= F | T
     amrentropy= F | T
     typeprolonglimit= 'default' | 'minmod' | 'woodward' | 'mcbeta' | 'koren'
     tfixgrid= DOUBLE
     itfixgrid= INTEGER
     ditregrid= INTEGER
    /

### mxnest, nxlone^D, dxlone^D, xprob^L {#par_mxnest}

`mxnest` indicates the maximum number of grid levels that can be used during the simulation, including the base grid level. It is an integer value which is maximally equal to the parameter _nlevelshi_ and minimally equal to 1 (which is the default value implying no refinement at all, but possibly a domain decomposition when the domain resolution is a multiple of the maximal grid resolution controlled by the -g= flag of $AMRVAC_DIR/setup.pl). The parameter _nlevelshi=8_ by default, a value set in _mod_indices.t_, so that if more than 8 levels are to be used, one must change this value and recompile. Note that when _mxnest&gt;1_, it is possible that during runtime, the highest grid level is temporarily lower than mxnest, and/or that the coarsest grid is at a higher level than the base level.

The computational domain is set by specifying the minimal and maximal
coordinate value per direction in the _xprob^L_ settings. When cylindrical or
spherical coordinates are selected with typexial, the angle ranges (for phi in
the cylindrical case, and for both theta and phi in the spherical case) are to
be given in 2 pi units.

The base grid resolution (i.e. on the coarsest level 1) is best specified by
providing the number of grid cells per dimension to cover the full
computational domain set by the _xprob^L_ ranges. This is done by specifying
these numbers in _nxlone^D_, where there are as many integers to be set as the
dimension of the problem. `Note that it is necessary to have a consistent
combination of base grid resolution and the _-g=_ setting for
_$AMRVAC_DIR/setup.pl_: the latter specifies the maximal individual grid
resolution which includes ghost cells at each side, while the _nxlone^D_
resolution must thus be a multiple of the individual grid resolution without
the ghost cells included`. An alternative way to specify the domain
resolution is to give the base grid cell size directly using _dxlone^D_, but
then again the same restrictions apply, and one must be sure that the step
size properly divides the domain size (from the _xprob^L_ pairs). It is
advocated to always use the integer setting through _nxlone^D_, from which the
code automatically computes the corresponding _dxlone^D_.

### errorestimate, nbufferx^D, skipfinestep, amr_wavefilter {#par_errest}

The code offers various choices for the error estimator used in automatically
detecting regions that need refinement.

When `errorestimate=0`, all refinement will only be based on the user-
defined criteria to be coded up in subroutine _specialrefine_grid_.

When `errorestimate=1`, refinement at time t_n will be based on a Richardson
procedure where two low order (using the typelow1 discretization) predictions
of the future time level t_(n+1) solution are computed and compared. The two
solutions are obtained by reversing the order of integrating and coarsening.
Note that only unsplit sources are taken into account in this process, and
that always a dimensionally unsplit scheme of first order is used (in
practice, always tvdlf1). If this comparison should rather be done on two such
solutions at time t_n itself, one must set the `skipfinestep=F`.

When `errorestimate=2`, we simply compare the previous time level t_(n-1)
solution with the present t_n values, and trigger refinement on relative
differences.

When `errorestimate=3`, the default value, errors are estimated using
current t_n values and their gradients following Lohner's prescription. In
this scheme, the `amr_wavefilter` coefficient can be adjusted from its
default value 0.01d0. You can set different values for the wavefilter
coefficient per grid level. This error estimator is computationally efficient,
and has shown similar accuracy to the Richardson approach on a variety of test
problems.

In all three latter cases, a call to the user defined subroutine
_specialrefine_grid_ follows the error estimator, making it possible to use
this routine for augmented user-controlled refinement, or even derefinement.

Depending on the error estimator used, it is needed or advisable to
additionally provide a buffer zone of a certain number of grid cells in width,
which will surround cells that are flagged for refinement by any other means.
It will thus trigger more finer grids. `nbufferx^D=2` is usually sufficient.
It can never be greater than the grid size specified with the -g setting of
_$AMRVAC_DIR/setup.pl_. For Lohner's scheme, the buffer can actually be turned
off completely by setting `nbufferx^D=0` which is default value.

### flags, wflags, logflag, tol, tolratio {#par_flags}

In all errorestimators mentioned above (except the errorestimate=0 case), the
comparison or evaluation is done only with a user-selected (sub)set of the
conserved variables. As there are _nw_ variables (which may include auxiliary
variables such as predefined in the special relativistic modules), the number
of variables to be used in the estimator is set by specifying it in
_flags(nw+1)_. Correspondingly, the first _flags(nw+1)_ entries for the flags
array select their index number, and corresponding weights are to be given in
wflags. E.g., 3D MHD has _nw=8_, and one can select refining based on density,
energy, and first component of the magnetic field by setting flags(9)=3,
flags(1)=1, flags(2)=5, flags(3)=6, since we predefined the order _rho m1 m2
m3 e b1 b2 b3_. The weights wflags(1), wflags(2), and wflags(3), must be
positive values between 0.0d0 and 1.0d0, and must add to unity. By default,
only the first conserved variable (typically density) is used in the
comparison.

The Lohner error estimation (errorestimate=3) may also decide to use
differences in the log10(rho), and this is then done by setting the
_logflag(1)=T_. This can be done per selected variables involved in the
estimation, but obviously only works for those that remain positive
throughout.

In the comparison involving the above selected variables, when the total error
exceeds the value set by `tol`, the grid is triggered for refining.
Reversely, if the error drops below `tolratio * tol`, the grid is coarsened.
`The user must always set a (problem dependent) value for `tol` (below 1),
while the default value for tolratio=1.0d0/8.0d0 has shown to be a rather
generally useful value. You can set tolerance values that differ per
refinement level. `

### iprob {#par_iprob}

As a possible integer switch for selecting multiple problem setups in the same
executable code, the integer switch `iprob` is provided. It is meant to be
used only in the user-written subroutines, for switching between e.g. multiple
initial conditions for the same executable.

### prolongprimitive, coarsenprimitive, restrictprimitive, amrentropy {#par_prolongprim}

It is possible to enforce the code to use primitive variables when coarsening
grid information (restriction), or filling new finer level grids
(prolongation). They are then used instead of the conservative variables,
which may not be a wise choice, but is perhaps better behaved with respect to
positivity of pressure etc. This is activated seperately for prolongation by
`prolongprimitive=T`, and for restriction by `restrictprimitive=T`. Also
in the Richardson error estimation process (when used), a coarsened grid is
created and this can be filled using the primitive variables when
`coarsenprimitive=T`. Again, this is to be used with care.

A much better strategy for handling positivity issues is offered by means of
`amrentropy=T`. If one selects `amrentropy=T`, the conserved energy
variable is replaced by the conserved entropy density (for HD and MHD only
now) (and both restrictprimitive and coarsenprimitive are overruled). The
switch from energy to entropy density is a general strategy, and is then used
in all restriction and prolongation phases. By default, all these flags are
set to false. For HD and MHD where the switch from energy to entropy is
already precoded in their _amrvacphys.t.EQUATION_ module, we advocate the use
of `amrentropy=T`.

The parameters `tfixgrid, itfixgrid` are available to fix the AMR
hierarchical grid from a certain time onwards (tfixgrid) or iteration (the it-
counter for the timesteps) onwards (itfixgrid). This may be handy for steady-
state computations, or for those cases where you know that the initial
conditions and physical setup is such that the AMR structure at t=0 will be
optimal for all times.The parameter `ditregrid` is introduced to reconstruct
the whole AMR grids once every ditregrid iteration(s) instead of regridding
once in every iteration by default.

## Paramlist {#par_paramlist}

         &paramlist
    	dtpar= DOUBLE
    	courantpar= DOUBLE
     typecourant= 'maxsum' | 'summax' | 'minimum'
    	dtdiffpar= DOUBLE
    	dtTCpar= DOUBLE
    	slowsteps= INTEGER
     time_accurate= T | F
     cfrac= DOUBLE

    /

### dtpar, courantpar, typecourant, dtdiffpar, dtTCpar, slowsteps {#par_dt}

If `dtpar` is positive, it sets the timestep `dt`, otherwise
`courantpar` is used to limit the time step based on the Courant condition.
The default is `dtpar=-1.` and `courantpar=0.8`.

For resistive MHD, the time step is also limited by the diffusion time: `dt
&lt; dtdiffpar*dx^2/eta`. The default is `dtdiffpar=0.5`. Further
restrictions on the time step can be put in the `getdt_special` subroutine in
the [AMRVACUSR module](@ref amrvacusr_specialsource). The library routines for
viscosity and div B diffusive cleaning, all use the coefficient `dtdiffpar` in
their stability conditions. `dtTCpar` with default value of 0.5 limits the
time step of thermal conduction. The `typecourant='maxsum'` means that the
time step limit for the CFL conditions takes the maximum over a summed
contribution to the maximum physical propagation speed for all dimensions. The
detailed formulae are found in setdt.t.

If the `slowsteps` parameter is set to a positive integer value greater than
1, then in the first `slowsteps-1` time steps `dt` is further reduced
according to the

                 2
        dt'= dt * [ 1 - (1-step/slowsteps) ]

formula, where `step=1..slowsteps-1`. This reduction can help to avoid
problems resulting from numerically unfavourable initial conditions, e.g. very
sharp discontinuities. It is normally inactive with a default value -1.

### time_accurate {#par_timeacc}

For steady state calculations, a grid-dependent value for the time step can be
used if the temporal evolution is not interesting. This local time stepping
strategy is presently under construction, and hence we always have
`time_accurate=T`.

### cfrac {#par_cfrac}

This is again specific for a cooling module for HD and MHD, as developed by AJ
van Marle, described in [mpiamrvac_radcool.md](mpiamrvac_radcool.md).
