# The VAC preprocessor

This document briefly describes the use of the VAC PreProcessor (VACPP) which
converts dimension independent notation into Fortran 90. VACPP is a
specialized implementation of the general [LASY Preprocessor](http://www-
personal.umich.edu/~gtoth/Papers/lasy.html).

VACPP is implemented in Perl in the **src/vacpp.pl** script. The main
variables of VACPP, the number of dimensions **$ndim**,
are normally modified by another Perl script
**src/setup.pl**. Based on these variables and the LASY patterns in the
`.t` source files VACPP generates the Fortran 90 source code of the output
`.f` files.

The preprocessor is mainly used via the makefile, but one can translate a single
file directly, or even use the preprocessor interactively, when the dimension
independent notation is typed in line by line from the keyboard and the expanded
code appears on the screen. The interactive use is a very efficient way of
checking the syntax of complex dimension independent notation. 

You may use **vacpp.pl** for translation directly

    vacpp.pl -d=3 FILENAME.t > FILENAME.f

You may change the maximum line length to e.g. 72 directly on the command line

    vacpp.pl -d=2 -maxlen=72 FILENAME.t > FILENAME.f

You can call vacpp.pl interactively with

    vacpp.pl -d=2 -

and then type a line of code to see how it is translated.
