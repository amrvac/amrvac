# Variable names

## Dimensions, Directions, and Limits

name | description
---|---
ndim	| number of spatial dimensions (1, 2, 3)
^D	| expands to 1,..,ndim
^ND	| expands to the value of ndim
ndir	| number of vector components (1, 2, 3)
^LIM	| expands to min,max
^L	| expands to min1,min2,min3,max1,max2,max3
^S	| expands to min1:max1,min2:max2,min3:max3

## Differential Equations and Physical Quantities

name | description
---|---
`pw(igrid)%%w`    |	w - conservative variables
`pw(igrid)%%wold` |	w at beginning of current time step
`pw(igrid)%%x`    |	x - spatial variables
global_time	| time
it | integer counter of temporal iterations

## Indices of variables in w array:

name | description
---|---
rho_	| density
e_ | energy density
p_ | pressure (primitive form of e_, same index)
mom(1:ndir) | momentum density
mag(1:ndir)	| magnetic field

## Coordinate system

name | description
---|---
r_	| R   direction for polar or cylindrical symmetry
phi_    | PHI direction for polar or cylindrical symmetry
z_      | Z   direction for polar or cylindrical symmetry

## Constants

name | description
---|---
dpi         | the value of pi (double precision)
smalldouble | a small number above machine precision
bigdouble   | a very large number

## General Tokens

name | description
---|---
i*	| index for *
h*	| i*-1  Previous index
j*	| i*+1  Next index
n*	| number of *. Usually refers to an index i*=1..n*.
d*	| difference of *, derivative of *
grad*   | spatial derivative of *
l*	| limited value
*C	| add 1/2 to spatial index/indices
*CT	| add 1/2 to temporal index
*R	| on the right
*L	| on the left
*I	| input  (eg. ixI^L are the limits of useful part of input)
*O	| output (eg. ixO^L are the limits of useful part of output)
*2	| on the power two, squared
sqr*	| square root of
kr	| Kronecker delta 3*3 tensor
lvc	| Levi-Civita 3*3*3 tensor
q*	| local variables to avoid name conflict with common ones
usr*    | user defined subroutines and parameters
*_	| indexname, integer parameters defined for legible code

## Grid Related Tokens

name | description
---|---
G	| grid (boundary ghost cells included)
M	| mesh (without boundary layers)
B	| boundary
nghostcells | block boundary cells width in number of cells
volume	| volume of mesh
dvolume	| volume of cell
surface	| surface of cell edge
normal	| normalized normal vector of cell surface 

## I/O Related Tokens

name | description
---|---
file	| input/output files
read    | reading from parameter or data files
save	| writing into output files
last	| last occasion (e.g. for saving into file)
*ini	| file containing initial conditions
*out	| file containing full results
head	| header string in files
unit	| units for open, close, read and write statements
unitterm| standard output

## Method Related Tokens

name | description
---|---
type*	| types (e.g. boundary, flow variables)
par	| parameters
eq	| equation
fix	| fixing certain problems (e.g. divbfix)
error	| error quantification
err	| same as error
courant	| Courant condition
step	| step in multistep methods
advance | advancing from time t to t+dt
advect  | advecting in directions idimmin..idimmax
advect1 | a single step of a multistep time integration in advect

## Subroutine Name Tokens

name | description
---|---
get*	| subroutine names for getting values (e.g. get_flux)
add*	| subroutine names for adding values  (e.g. add_source)

## Examples

For example the variable name

    ditsave = d + i + t + save

is the **d**ifference in the **i**ndices of the **t**imes between **save**s
into the output files, i.e. the frequency of saving snap shots measured in
time steps. Another example:

    ixGhi^D = i + x + G + hi + ^D

is the highest possible indices for the coordinates for the grid for each
dimensions.
