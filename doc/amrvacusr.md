# Setting up an new problem

[TOC]

# Introduction {#user_intro}

This document describes how users can create new problems. It explains how to
set initial conditions, and how to perform other customizations.

# Creating a New Setup {#user_setup}

User code for a new problem should go into a file called either:

* `mod_usr.t`, if it contains dimension-independent code (Fortran+LASY)
* `mod_usr.f`, if it contains regular Fortran code

A template for this file is available in `src/mod_usr_template.t`. To create a
new 1D/2D/3D setup, run one of these commands:

    $AMRVAC_DIR/setup.pl -d=1
    $AMRVAC_DIR/setup.pl -d=2
    $AMRVAC_DIR/setup.pl -d=3

This will copy an AMRVAC makefile and specify the problem dimension in it. If no
user code is found, it will ask whether you want to copy the default template.
Alternatively, you can look for an existing problem (look in `tests/`) and
customize its `mod_usr.t` or `mod_usr.f` file to get started. To create a normal
Fortran file from `mod_usr.t`, you can type `make mod_usr.f`.

# Structure of mod_usr.t {#user_structure}

\include mod_usr_template.t


# TODO

@todo Finish documentation for user code, par files

## Specialbound part {#user_specialbound}

When the predefined boundary types provided by MPI-AMRVAC are not sufficient
the **specialbound** subroutine can solve the problem. It is called for each
boundary region and variable for which the boundary type
**typeB='special'...** is selected. It is likely that you can use the
predefined boundary types for at least some of the variables and/or regions.
The **specialbound** subroutine is called _after_ the variables with
predefined boundary types have been updated in the ghost cells. The calling
interface for the specialbound subroutine is found in the
_amrvacnul.specialbound.t_ file.

An example for the use of this _specialbound_usr_ subroutine is found in the
example user file **usr/amrvacusr.t.wchd22**, which realizes the standard 2D
hydro Woodward and Collela shock reflection problem.

## Specialsource part {#user_specialsource}

There are lots of possible physical source terms for the same basic equation.
Rather than writing a new physics module for each, it is simpler to define a
problem dependent source term in the AMRVACUSR module. There are already a
number of implemented source terms in the  src/amrvacusr.LIBRARY.t files,
which can be studied for example cases.

The **specialsource** subroutine is called at least once in every time step by
MPI-AMRVAC. The number of calls depends on the time integration scheme defined
by **typeadvance**, the parameter **sourcesplit** and also on the number of
dimensions if the parameter **dimsplit=T**.

In any case the subroutine should integrate **w** from time **qt** to
**qt+qdt** for the variables listed in **iw^LIM** in the region **ixO^S**. The
source terms should be evaluated for the **wCT** array, which corresponds to
the physical time **qtC**. In case of explicit time dependence, **qtC** should
be used as time. Only elements within **ixI^S** can be used from **wCT**.

The **getdt_special** subroutine can limit **dt** for numerical stability if
the source term requires that. This subroutine is called after the CFL
condition **getdt_courant** and the **getdt** subroutine of the AMRVACPHYS
module have been executed.

The **specialeta** subroutine is used for MHD to set the resistivity array
**eta** when it is not constant in time and/or space. The **current** array
must then be computed, when anomalous resistivity is described as a function
of **J**.

The **specialrefine_grid** subroutine allows to add user controlled
(de)refinement, by setting the integers _refine,coarsen_. You have all info
available to do this refining (grid level, physical values in _w_, coordinates
in _x_, time in _qt_). Similarly, the **specialvarforerrest** subroutine
allows to compute a (local) new variable, to be stored in _var_, on which you
can then base refinement as well. This is true for the lohner error estimator
only.

## Speciallog part {#user_speciallog}

The _amrvacnul/speciallog.t_ file contains additional subroutines more related
to special I/O requests. The default log-file may be altered, for which you
need to code up the _printlog_special_ subroutine. For parallel execution,
this invariably means the use of MPI constructs, so you should copy in the
default version from _amrio.t_ and then study it, and modify accordingly.

The _process_grid_usr_ is a subroutine which allows to compute auxiliary
variables which happen to be non-local (like div v), and are in no way used
for flux computations. As auxiliaries, they are also not advanced. This
functionality was added to allow for separate particle treatments using
stochastic differential equations, where the particle dynamics is only relying
on local compression values etc.

The _specialvar_output_ is extremely handy to compute variables from the
actually computed conserved variables, that can then be visualized directly.
It is only used in combination with the conversion subroutines. E.g., one may
here compute current density components using the actual code discretizations
for computing a curl, and then visualize those with any of the visualization
tools applicable. You then also need to specify a label for this variable, in
_specialvarnames_output_.
