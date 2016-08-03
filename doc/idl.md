# IDL visualization macros

# Introduction

This document describes the use of IDL for plotting results from VAC/MPI-
AMRVAC. IDL is not a free software, it is distributed by Research Systems Inc.
**You should check the VAC manual when really handling uniform grid data, it
is more extensive. The functionality for amr data is limited to 1D and some 2D
stuff. Some of the descriptions below may be disfunctional for true AMR
data.**

The IDL macros provided are all in the **Idl_amr** and **Idl** directories
except the **.idlrc** startup file, which is supposed to be in the main
directory.

# Startup file

It is useful to let IDL know about the existence of the startup file
**.idlrc** in the main MPI-AMRVAC/ directory. Under UNIX put

    setenv IDL_STARTUP .idlrc

into your ~/.login or ~/.cshrc file if you use the csh or tcsh shell. For
other UNIX shells, use

    IDL_STARTUP=.idlrc
    export IDL_STARTUP

in the **.profile** file. Read the IDL manual if your operating system is not
UNIX. If you start IDL in a directory where there is no .idlrc, IDL will give
a warning message.

# Running IDL

Assuming that the IDL_STARTUP variable is properly set, simply start IDL as

    idl

in the main MPI-AMRVAC directory. If IDL_STARTUP is not set, type

    @.idlrc

at the IDL&gt; prompt or start IDL like this

    idl .idlrc

In either case, the commands in the **.idlrc** file are executed: some global
variables like !path are set, and some procedures in
**Idl(_amr)/procedures.pro, Idl(_amr)/animfunc.pro** and **Idl/vel.pro** are
compiled. The script **Idl_amr/defaults** is also executed to set some global
variables to their default values. You can customize the startup of IDL by
editing **.idlrc** and **Idl_amr/defaults**, e.g. you can compile your own IDL
subroutines.

If you get trapped by an error inside some IDL routine,

    retall

will return to the main level. To exit IDL type

    exit

# Reading a Snapshot with getpict

To read a single frame from a file type at the "IDL&gt;" prompt of IDL

    IDL> .r getpict

The procedure will prompt you for the **filename**, and it determines the
**filetype** and **npictinfiles** (the number of snapshots in the file)
automatically in case it is a VAC file, or detect that it is an AMR file (for
which it can not determine the number of frames from the (variable)
filelength). Then it asks for the frame-number **npict** (1, 2,...
npictinfiles) of the snapshot to be read from the file. When
**npictinfiles=1**, the frame number is set to 1 automatically:

    filename(s)   ? data/exampleA22.ini
    filetype(s)   = ascii
    npictinfile(s)=      1
    npict=       1

The header of the file is read and echoed on the screen, and the type of
equation **physics** is asked, unless it is read from the headline of the file
given as e.g. **ExampleA_hdadiab22**:

    headline  =ExampleA
    ndim      = 2, neqpar= 3, nw= 3
    nx        =    50   50
    eqpar     =       2.0000000       1.0000000       0.0000000
    variables = x y h m1 m2 gamma ghalf coriolis
    physics (eg. mhd12)=hdadiab22
    it      =           0, time=       0.0000000
    Read x and w
    GRID            INT       = Array(50, 50)

At the end, the **x** and **w** variables (containing the coordinates and the
conservative variables respectively) are read from the file.

If your file contained data in Cartesian coordinates, you get the "IDL&gt;"
prompt back, and you can do whatever you want with **x**, **w**, and all the
other variables **headline, it, time, gencoord, ndim, neqpar, nw, nx, eqpar,
variables** defined by the header.

You can use the **plotfunc** script to get some sophisticated plots or you can
use any of the IDL procedures directly, e.g.

    print,time,it
    print,nx
    print,variables
    plot,w(*,10,0)
    surface,w_1(*,*,2)-w_0(*,*,2)
    contour,w(*,*,1)
    vel,w(*,*,1),w(*,*,2)

The functions **diff2, diff3**, and **diff4** are provided for taking spatial
derivatives of variables represented on 2D Cartesian mesh. To calculate the
curl of the velocity field you could type

    vx=w(*,*,1)/w(*,*,0)
    vy=w(*,*,2)/w(*,*,0)
    xx=x(*,*,0)
    yy=x(*,*,1)
    curl_v=diff2(1,vy,xx)-diff2(2,vx,yy)

The difference between diff2 and diff3 is minor: diff2 uses centered
differences, while diff3 relies on IDL built-in procedures. **diff4** is a
fourth order accurate centered difference formula. You can use these functions
in the **Idl/animfunc.pro** file to define functions plotted by **plotfunc**
and **animate**, for example the function 'divb4' uses the diff4 function.

For non-Cartesian 2D grids the functions **grad, div, curl** and **grad_rz,
div_rz, curl_rz** are provided, which use contour integrals to estimate the
cell averaged gradient, divergence, and curl. The subscripts **_rz** refer to
cylindrical symmetry in the ignored PHI direction. All these functions are
described in detail in **Idl/procedures.pro**, and some examples for their use
can be found in **Idl/animfunc.pro**. For example, 'curlv' is defined with the
**curl** function. **Note that the above may not work properly for AMR data,
it does work with single grid VAC files**

# Plotting the Data

Once the data is read by getpict or animate you can plot functions of **w**
with

    .r plotfunc

You will see the value of **physics** for the last file read and some
parameters with standard default values for different plotting routines:

    physics (e.g. mhd12)      =hdadiab22
    ======= CURRENT PLOTTING PARAMETERS ================
    ax,az=  30, 30, contourlevel= 30, velvector= 200, velspeed (0..5)= 5
    multiplot= 0 (default), axistype (coord/cells)=coord
    bottomline=2, headerline=2

The viewing angle is given by **ax** and **az** for the plotting modes
**plotmode=** 'surface', 'shade', and 'shadeirr'. The **contourlevel**
parameter determines the number of contourlevels for the plotting modes
'contour', 'contlabel', and 'contfill'. The **velvector** and **velspeed**
parameters are used by the **vel** procedure to set the number of vectors
drawn and their speed during animation in **animate**. The number of subplots
is usually determined by the number of files and the number of functions, but
you can override this by setting e.g. **multiplot=[2,3,0]**, which gives 2 by
3 subplots filled in row-wise. If the third element is 1, the subplots are
filled in column-wise. Setting **multiplot=3** is identical with
**multiplot=[3,1,1]**. **multiplot=0** gives the default behaviour. The plots
are normally shown in physical coordinates, i.e. **axistype='coord'**, but the
axes can also run in cell indices if **axistype='cells'** is set. The
variables **bottomline** and **headerline** control the number of values shown
at the bottom from **time, it, nx** and at the top from **headline, nx**. You
can change these values explicitly (e.g. **az=50**), or change their default
values in **Idl/defaults**. Now, you will be prompted for the name of
function(s) and the corresponding plotting mode(s):

    ======= PLOTTING PARAMETERS =========================
    wnames                     =  h m1 m2
    func(s) (e.g. rho p m1;m2) ? h m1
    2D plotmodes: contour/contlabel/contfill/shade/shadeirr/surface/tv/
                  stream/vel/velovect/ovelovect'
    plotmode(s)                ? surface
    plottitle(s) (e.g. B [G];J)=default
    autorange(s) (y/n)         =y
    GRID            INT       = Array(50, 50)

The function(s) to be plotted are determined by the **func** string parameter,
which is a list of function names separated by spaces. The number of functions
**nfunc** is thus determined by the number of function names listed in
**func**. These function names can be any of the variable names listed in the
string array **wnames**, which is read from the header of the file, or any of
the function name strings listed in the **case** statements in the animfunc
function in **Idl/animfunc.pro**, or any expressions using the standard
variable, coordinate and equation parameter names valid for the particular
geometry and physics. These variable names are a subset of

    xx yy zz r phi z   rho m1 m2 m3 mm e b1 b2 b3 bb  gamma eta adiab csound2

where "mm" and "bb" are the momentum and magnetic field squared, respectively.
Note that "xx ... bb" are arrays, while "gamma ... csound2" are scalars. For
example **r*m2/rho** is a valid expression for the angular momentum if the
second momentum points in the phi direction, or the maximum Alfven speed could
be given as **func='sqrt(bb/rho)'**, but this is already defined in
animfunc.pro as 'calfven'.

You may combine two function names with the **;** character representing two
components of a vector, e.g. **v1;v2**, which can either be plotted as a
vectorfield by the **velovect** and **vel** procedures, or as streamlines,
using **stream** plotmode, or otherwise the absolute value
**sqrt(v1**2+v2**2)** is plotted by **surface, contour** etc. The
**ovelovect** choice is identical with **velovect** except for the axis
ranges, which are better suited for overplotting with a **contour** or
**contfill** plot, for example. You can also put a minus sign in front of any
function name, which will simply multiply the value of the rest of the string
by -1. For example '-T' plots (-1)*temperature. _Beware! func='-m1+m2' will
actually plot -(m1+m2), so use '(-m1+m2)' or 'm2-m1' instead_.

For each function you may set the **plotmode** (in 1D there is no choice,
since only the **plot** routine is useful). If you give fewer plotmode(s) than
**nfunc**, the rest of the functions will use the last plotmode given, in the
above example **surface**. This padding rule is used for all the arrays
described by strings. The **plottitle** parameter is usually set to
**default** which means that the function name is used for the title, but you
can set it explicitly, e.g. **plottitle='Density;Momentum'**. Here the
separator character is a semicolon, thus the titles may contain spaces. No
titles are produced if **plottitle=' '** is set.

For each function you may set the plotting range by hand or let IDL to
calculate the minimum and maximum by itself. This is defined by the
**autorange** string parameter, which is a list of 'y' and 'n' characters,
each referring to the respective function. If you set 'n' for any of the
variables, the **fmin** and **fmax** arrays have to be set, e.g.

    fmin=[1. ,-1.]
    fmax=[1.1, 1.]

Make sure that the values are floating point and not integer. IDL remembers
the previous setting and uses it, unless the number of functions are changed.
You can always set fmin=0, fmax=0, and let IDL prompt you for the values.

It is possible to plot a part of the simulation domain. You should use
domainset for that, in case of AMR data.

The number of functions and the number of subplots can be any combination you
would like. In 1D plots, the line style is varied for the different functions,
so the curves can be distinguished. Note, that **plottitle** is set to avoid
the default titles 'h' and 'v1;v2' overlap on top of the plot. The
**multiplot=N** setting is equivalent with **multiplot=[N,1,1]**. In the
second plot 'ovelovect' is used (instead of 'velovect') for the velocity to
get good alignment with the 'contfill' plot.

# Plotting another Snapshot

If you type

    .r getpict
    .r plotfunc

again, the data will be read and plotted again without any questions asked,
since IDL remembers the previous settings.

If you want to read another frame, say the second, from the same file, type

    npict=2
    .r getpict

You can change the **func** and **plotmode** variables the same way:

    func='rho p'
    plotmode='contour surface'
    .r plotfunc

Note that we did not need to reread the data. Other variables, all listed in
**Idl/defaults**, can be set similarly. If you set

    doask=1

the macros will ask for all the parameters to be confirmed by a simple RETURN,
or to be changed by typing a new value. Set **doask=0** to get the default
behaviour, which is no confirmation asked. To overplot previous plots without
erasing the screen, set

    noerase=1

You can return to the default settings for all parameters by typing

    .r defaults

# Plotting and Animation with animate

This general procedure can plot, save into file(s), or animate (using IDL's
Xinteranimate) different functions of data read from one or more files. If a
single snapshot is read, the plot is drawn without animation. In essence
**animate** combines **getpict** and **plotfunc** for any number of files and
any number of snapshots.

    .r animate

will first prompt you for **filename(s)** unless already given. Animating more
than one input files in parallel is most useful for comparing simulations with
the same or very similar physics using different methods or grid resolution.
It is a good idea to save snapshots at the same _physical_ time into the data
files, so use **dtsave** and **tsave** or fixed time steps **dtpar** in the
parameter file. The headlines and the grid sizes will be shown in the
headerline for each file separately above the corresponding plots.

The function(s) to be animated and the plotmode(s) for the functions are
determined by the same **func, plotmode**, and **plottitle** strings as for
plotfunc. If any part of the **autorange(s)** string is set to **'y'**, the
data file(s) will be read twice: first for setting the common range(s) for all
the snapshots and the second time for plotting. The first reading of the
file(s) also determines **npict**, which is the number of snapshots to be
animated. It is limited by the end of file(s) and/or by the **npictmax**
parameter. With a formula

    npict=min( npictmax, (min(npictinfiles)-firstpict)/dpict + 1 )

If **autorange='n'** for all the functions the file is only read once.

The animation runs from **firstpict**, every **dpict**-th picture is plotted
and the total number of animated frames is at most **npictmax**. You can save
the frames into a series of GIF files **Movie/1.gif,Movie/2.gif,...** by
setting **savemovie='y'**, or into PostScript files
**Movie/1.ps,Movie/2.ps,...** by setting **printmovie='y'**. The number of
plots saved into files is limited by **nplotmax**. The GIF files can be put
together into a movie by some program, which you check with your local
sysadmin.

The **multiplot** array can be used to get some really interesting effects in
**animate**. Besides overplotting different functions, as explained for
**plotfunc**, the data of different files can also be overplotted for
comparison purposes. Probably it is a good idea to compare 1D slices rather
than full 2D plots, e.g.

    filename='data/exampleA22.ini data/exampleB22.ini'
    func='h v1'
    cut=grid(*,25)
    multiplot=2
    .r animate

will overplot height and velocity read from the two files. The lines belonging
to the two data files are distinguished by the different line styles.
Overplotting two data sets is especially useful when the two results are
supposed to be identical.

Timeseries can also be produced easily with **multiplot**.

    filename='data/example22.out'
    func='h'
    plotmode='surface'
    npictmax=6
    multiplot=[3,2,0]
    headerline=2
    bottomline=1
    .r animate

will show the first 6 snapshots of height in a single plot. Now the time is
shown for each plot individually, and setting **bottomline=1** limits the time
stamp to the most essential information, time. If **npict*nfile*nfunc** is
greater than the number of subplots defined by **multiplot**, an animation is
done. This can be combined with **printmovie='y'** to produce PostScript
figures containing 6 subplots each, which is convenient for printing an
animation. Type

    multiplot=0

to return to default behavior, which is one snapshot per plot.

Even after exiting from Xinteranimate, the animation can be repeated again
without rereading the data file(s) by typing

    xinteranimate,/keep_pixmaps

# Function Definitions in animfunc

Any function of the conservative variables **w**, the coordinates **x**, and
the equation parameters **eqpar** can be defined in the **Idl/animfunc.pro**
file. Extra information is provided by the equation type **physics** and the
variable names **wnames**. The function is identified by the string **f**. The
user can easily define new functions for a specific application by adding a
new **case** statement, for example the accretion rate in cylindrical symmetry
can be defined as

          f eq 'rmr':  result=w(*,*,1)*r

Make sure that the function is in the right branch of the **ndim** switch.
Note the use of the **r** array, which is extracted from **x(*,*,0)** for
convenience. The modified file should be recompiled and can be used as

    .r animfunc
    func='rmr'
    .r plotfunc

A few frequently used functions are defined in a rather general way, they work
for all the appropriate equations and number of dimensions. ** Note: the list
below is for VAC, you should check yourself whether all below still works
properly in AMR runs!**

    Function name	Physics	NDIM	Meaning
    ------------------------------------------------------------
    v		any 1D		velocity in 1st direction
    vx, v1		any		velocity in 1st direction
    vy, v2		any		velocity in 2nd direction
    vz, v3		any		velocity in 3rd direction
    p               mhd, mhdiso	total pressure
    p		hd, hdadiab	thermal pressure
    pth		any		thermal pressure
    pbeta           mhd, mhdiso     plasma beta: 2*p_thermal/B^2
    T		any		temperature
    s               any             entropy: p_thermal/rho^gamma
    csound		any		sound speed: sqrt(gamma*p_thermal/rho)
    cslow1		mhd, mhdiso     slow magnetosonic speed along 1st dimension
    cslow2		mhd, mhdiso     slow magnetosonic speed along 2nd dimension
    cslow3		mhd, mhdiso     slow magnetosonic speed along 3rd dimension
    calfven1	mhd, mhdiso     Alfven speed along 1st dimension: b1/sqrt(rho)
    calfven2	mhd, mhdiso     Alfven speed along 1st dimension: b2/sqrt(rho)
    calfven3	mhd, mhdiso     Alfven speed along 1st dimension: b3/sqrt(rho)
    calfven		mhd, mhdiso     maximum of Alfven speed: |B|/sqrt(rho)
    cfast1		mhd, mhdiso     fast magnetosonic speed along 1st dimension
    cfast2		mhd, mhdiso     fast magnetosonic speed along 2nd dimension
    cfast3		mhd, mhdiso     fast magnetosonic speed along 3rd dimension
    cfast		mhd, mhdiso     maximum of fast speed: sqrt(csound^2+calfven^2)
    mach 		any		Mach number: |v|/csound
    mach1		any		Mach number:  v1/csound
    mach2		any		Mach number:  v2/csound
    mach3		any		Mach number:  v3/csound
    Mslow1		mhd, mhdiso     slow Mach number along 1st dimension: v1/cslow1
    Mslow2		mhd, mhdiso     slow Mach number along 2st dimension: v2/cslow2
    Mslow3		mhd, mhdiso     slow Mach number along 3st dimension: v3/cslow3
    Malfven1	mhd, mhdiso     Alfven Mach numbe along 1st dim: v1/calfven1
    Malfven2	mhd, mhdiso     Alfven Mach numbe along 2nd dim: v2/calfven2
    Malfven3	mhd, mhdiso     Alfven Mach numbe along 3rd dim: v3/calfven3
    Malfven		mhd, mhdiso     maximum Alfven Mach number: |v|/calfven
    Mfast1		mhd, mhdiso     fast Mach number along 1st dimension: v1/cfast1
    Mfast2		mhd, mhdiso     fast Mach number along 2nd dimension: v2/cfast2
    Mfast3		mhd, mhdiso     fast Mach number along 3rd dimension: v3/cfast2
    Mfast		mhd, mhdiso     maximum fast Mach number: |v|/cfast
    machA		mhd, mhdiso	same as Mfast
    curlv		any 2D		curl of velocity
    j		mhd, mhdiso 2D	current for slab symmetry
    j_rz		mhd, mhdiso 2D  current for cylindrical symmetry
    j_rp		mhd, mhdiso 2D	current for polar coordinates
    divb		mhd, mhdiso 2D	div B for slab symmetry
    divb4		mhd, mhdiso 2D	div B fourth order for slab symmetric
                                    uniform Cartesian grid
    divb_CT         mhd, mhdiso 2D  div B for CT schemes on Cartesian grid
    divb_CD         mhd, mhdiso 2D  div B for CD schemes on generalized grid
    divb_rz		mhd, mhdiso 2D	div B for for cylindrical symmetry
    divb_rp		mhd, mhdiso 2D	div B for for polar coordinates
    A, AA, B	mhd, mhdiso 2D	vector potential for slab symmetry
    A_r, AA_r	mhd, mhdiso 2D	vector potential*r for cylindrical symmetry

Here _any_ physics means any of _hdadiab, hd, mhdiso, mhd_.

The vector potential A (or A*r in cylindrical coordinates) is most useful with
**plotmode='contour'** to plot the magnetic field lines. (In case of A*r the
contours are parallel to the field lines but their density is not proportional
to the field strength). It is calculated by integrating the components of
**B** along rows and columns of the possibly irregular grid. The vector
potential is well defined only if divergence B is zero, or negligible
everywhere. Since this is not always true, there are a few options to
calculate the vector potential. 'A' and 'A_r' are integrated upwards starting
from the bottom, while 'AA', 'B', and 'AA_r' are integrated from left to
right. This is important if the computational domain is irregular, and the
interpolated grid contains some unphysical regions, which should be at the end
of the integration paths. Finally 'B' differs from 'AA' in that the resulting
vector potential is smoothed a bit after the integration, which makes the
field lines look nicer.

# Reading the Logfile with getlog

One or more (at most three) logfiles can be read by

    .r getlog

which reads data from the file(s) determined by the **logfilename** parameter.
An initial guess for the name is made if the **filename** parameter has
already been given. The data in the logfile(s) is put into the
**step,t,dt,wlog** arrays in case of a single file, and into
**step0,t0,dt0,wlog0, step1,t1,dt1,wlog1** in case of two files. If the
[residual](par.html#Stoplist) is calculated and stored in the last column of
the logfile, the **resid** is also generated. The **wlog(nt,nwlog)** array
contains the rest of the columns in the logfile. A simple example is

    .r getlog
    logfilename(s) ? data/example22.log
    headline       =ExampleA_hdadiab22
    Reading arrays step,t,dt,wlog:
      wlog(*, 0)= h
      wlog(*, 1)= m1
      wlog(*, 2)= m2
    Number of recorded timesteps: nt=      20
    plot,t,wlog(*,0),xtitle='t',ytitle='h_mean'

which checks the global mass conservation for exampleA.

# Saving Plots into Postscript Files

In IDL printing a plot is possible through Postscript files. After the plot
looks fine on the screen, use for example

    set_plot,'PS'
    device,filename='myfile.ps',xsize=24,ysize=18,/landscape,/color,bits=8
    loadct,3
    .r plotfunc
    device,/close
    set_plot,'X'

For a non-color plot omit the **/color,bits=8** parameters and the loading of
the color table by the **loadct** command. For a _portrait_ picture use
**xsize=18,ysize=24** and omit the **/landscape** keyword. If the printout is
off the page, set **yoffset=3** too. You can use .r animate instead of .r
plotfunc (e.g. for multiple files or for time series), but make sure that only
one plot is produced by setting **npictmax=1**, and use **firstpict** to
select the snapshot. To save each frame of an animation into a different
Postscript file, set the **printmovie='y'** parameter as described above. More
snapshots can also be saved in a single plot using **multiplot**.

All these commands can be collected into a single file, like
**Idl/myfig.pro**, which can be run from IDL by

    @myfig

This is a convenient way to store the commands for producing complicated
figures.
