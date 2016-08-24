# The VAC preprocessor

This document briefly describes the use of the VAC PreProcessor (VACPP) which
converts dimension independent notation into Fortran 90. VACPP is a
specialized implementation of the general [LASY Preprocessor](http://www-
personal.umich.edu/~gtoth/Papers/lasy.html).

VACPP is implemented in Perl in the **src/vacpp.pl** script. The main
variables of VACPP, such as the number of dimensions **$ndim** and vector
components **$ndir**, are normally modified by another Perl script
**src/setamrvac**. Based on these variables and the LASY patterns in the
`.t` source files VACPP generates the Fortran 90 source code of the output
`.f` files.

The maximum length of the lines is determined by the _final_ translation step
and it can be set by the `-maxlen=...` flag. By default the output lines are
at most 78 character long. Most compilers have flags to accept such line
length, but if this is not the case, you may set '-maxlen=72' for the
preprocessor step in `src/makefile`.

The preprocessor is mainly used via the makefile, but one can translate a single
file directly, or even use the preprocessor interactively, when the dimension
independent notation is typed in line by line from the keyboard and the expanded
code appears on the screen. The interactive use is a very efficient way of
checking the syntax of complex dimension independent notation. For interactive
usage the number of dimensions and vector components can be overdefined
temporarily with the `-d` flag. For example type (**note the last dash!**)

    vacpp.pl -d=23 -
    ^D
    1,2
    ^C
    1,2,3
    x(ixI^S,1)
    x(ixImin1:ixImax1,ixImin2:ixImax2,1)
    ...
    Ctrl-D

You may use **vacpp.pl** for translation directly

    vacpp.pl FILENAME.t > FILENAME.f

or indirectly via the makefile with

    make FILENAME.f

You may change the maximum line length to e.g. 72 directly on the command line

    vacpp.pl -maxlen=72 FILENAME.t > FILENAME.f

or edit the **src/makefile** accordingly.

