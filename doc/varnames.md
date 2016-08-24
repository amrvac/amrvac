# Variable names

This document describes how variable names are constructed in the source files
of MPI-AMRVAC. This was done rigorously in VAC, and in MPI-AMRVAC we have kept
most of the naming conventions, but introduced new structure to handle the AMR
grid hierarchy. The patterns starting with **^** are expanded by the VACPP
preprocessor. The meaning of some tokens is given below.

## Differential Equations and Physical Quantities

    pw(igrid)%w  	:	w	- conservative flow variables
    pwold(igrid)%w	:	w       - possibly at previous time step
    pwres(igrid)%w	:	residual- relative change in w
    px(igrid)%x	:	x	- spatial variables
    
    fC	- flux, defined at cell edges
    
    t	- time
    it	- integer time counter
    r_	- R   direction for polar or cylindrical symmetry
    phi_    - PHI direction for polar or cylindrical symmetry
    z_      - Z   direction for polar or cylindrical symmetry
    
    rho	- density
    p	- pressure
    pthermal- thermal pressure
    ptotal  - total pressure including magnetic pressure
    v       - velocity
    csound  - sound speed
    b	- magnetic field
    divb	- divergence of B
    current - current density
    res     - resistivity
    eta     - resistivity coefficient
    
    a	- dimensionless propagation speed or eigen value
    il	- index for characteristic variables or eigen values
    *W_	- index name for characteristic waves
    jump	- jump of characteristic variables
    roe	- Roe's or just arithmetic average in Roe solver
    phi	- the limiter function for characteristic variable
    
    zero	- 0.D0
    one	- 1.D0
    two	- 2.D0
    half	- 0.5D0
    quarter - 0.25D0
    dpi     - the value of pi (double precision)
    
    smalldouble -a small number above machine precision
    bigdouble   -a very large number

## Dimensions, Directions, and Limits

    dim	- dimension or direction
    ^D	- expands to 1,..,ndim
    ^ND	- expands to the value of ndim
    dir	- direction (of vector variables)
    ^C	- expands to 1,..,ndir (Component)
    ^NC	- expands to the value of ndir
    *min	- minimum of * (expected, or actual)
    *max	- maximum of * (expected, or actual)
    ^LIM	- expands to min,max
    ^L	- expands to min1,min2,min3,max1,max2,max3
    ^S	- expands to min1:max1,min2:max2,min3:max3
    *lo	- lowest  value of * (enforced by array boundaries)
    *hi	- highest value of * (enforced by array boundaries)
    ^LLIM	- expands to lo,hi
    ^LL	- expands to lo1,lo2,lo3,hi1,hi2,hi3
    ^T	- expands to lo1:hi1,lo2:hi2,lo3:hi3

## General Tokens

    i*	- index for *
    h*	- i*-1  Previous index
    j*	- i*+1  Next index
    n*	- number of *. Usually refers to an index i*=1..n*.
    d*	- difference of *, derivative of *
    grad*   - spatial derivative of *
    l*	- limited value
    *C	- add 1/2 to spatial index/indices
    *CT	- add 1/2 to temporal index
    *R	- on the right
    *L	- on the left
    *I	- input  (eg. ixI^L are the limits of useful part of input)
    *O	- output (eg. ixO^L are the limits of useful part of output)
    *2	- on the power two, squared
    sqr*	- square root of
    kr	- Kronecker delta 3*3 tensor
    lvc	- Levi-Civita 3*3*3 tensor
    ^LEN*	- length of string for *
    q*	- local variables to avoid name conflict with common ones
    special*- user defined subroutines and parameters
    *_	- indexname, integer parameters defined for legible code

## Grid Related Tokens

    G	- grid (boundaries included)
    grid	- grid (boundaries included)
    M	- mesh (without boundary layers)
    mesh	- mesh (without boundary layers)
    B	- boundary
    dixB	- boundary width
    volume	- volume of mesh
    dvolume	- volume of cell
    surface	- surface of cell edge
    normal	- normalized normal vector of cell surface 

## I/O Related Tokens

    file	- input/output files
    read    - reading from parameter or data files
    save	- writing into output files
    last	- last occasion (e.g. for saving into file)
    *ini	- file containing initial conditions
    *out	- file containing full results
    head	- header string in files
    name	- names (e.g. filename, variables, w)
    wnames	- names of all w variables (e.g. 'rho m1 m2')
    unit	- units for open, close, read and write statements
    unitterm- standard output

## Method Related Tokens

    type*	- types (e.g. boundary, flow variables)
    par	- parameters
    eq	- equation
    fix	- fixing certain problems (e.g. divbfix)
    error	- error quantification
    err	- same as error
    courant	- Courant condition
    step	- step in multistep methods
    pred	- predictor (or half) step
    full	- full step
    advance - advancing from time t to t+dt
    advect  - advecting in directions idimmin..idimmax
    advect1 - a single step of a multistep time integration in advect
    method	- spatial discretization method

## Subroutine Name Tokens

    get*	- subroutine names for getting values (e.g. getflux)
    add*	- subroutine names for adding values  (e.g. addsource)

## Examples

For example the variable name

    ditsave = d + i + t + save

is the **d**ifference in the **i**ndices of the **t**imes between **save**s
into the output files, i.e. the frequency of saving snap shots measured in
time steps. Another example:

    ixGhi^D = i + x + G + hi + ^D

is the highest possible indices for the coordinates for the grid for each
dimensions.

\todo Document all the variables described here **inside the code**
