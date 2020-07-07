# Setting up a new problem

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

    setup.pl -d=1
    setup.pl -d=2
    setup.pl -d=3

This will copy an AMRVAC makefile and specify the problem dimension in it. If no
user code is found, it will ask whether you want to copy the default template.
Alternatively, you can look for an existing problem (look in `tests/`) and
customize its `mod_usr.t` or `mod_usr.f` file to get started. To create a normal
Fortran file from `mod_usr.t`, you can type `make mod_usr.f`.
Specify other user routines, for a list see mod_usr_methods.t

# Structure of mod_usr.t {#user_structure}

\include mod_usr_template.t

## Special boundary {#user_specialbound}

When the predefined boundary types provided by MPI-AMRVAC are not sufficient
, a pointed **usr_special_bc** subroutine can solve the problem. It is called 
for each boundary region and all variable in the boundary region
**typeboundary_min1=8*'special'...** is selected. 

## Special source {#user_specialsource}

There are lots of possible physical source terms for the same basic equation.
Rather than writing a new physics module for each, it is simpler to define a
problem dependent source term in the user module. 

A pointed **usr_source** subroutine will be called at least once in every time 
step. The number of calls depends on the time integration scheme defined
by **time_integrator**, the parameter **source_split_usr** and also on the 
number of dimensions if the parameter **dimsplit=T**.

In any case the subroutine should integrate **w** from time **qt** to
**qt+qdt** for the variables listed in **iw^LIM** in the region **ixO^S**. The
source terms should be evaluated for the **wCT** array, which corresponds to
the physical time **qtC**. In case of explicit time dependence, **qtC** should
be used as time. Only elements within **ixI^S** can be used from **wCT**.

A pointed **usr_get_dt** subroutine can limit **dt** for numerical stability if
the source term requires that. 

A pointed **usr_special_resistivity** subroutine is used for MHD to set the 
resistivity array **eta** when it is not constant in time and/or space.

A pointed **usr_refine_grid** subroutine allows to add user controlled
(de)refinement, by setting the integers _refine,coarsen_. You have all info
available to do this refining (grid level, physical values in _w_, coordinates
in _x_, time in _qt_). Similarly, a **usr_var_for_errest** subroutine
allows to compute a (local) new variable, to be stored in _var_, on which you
can then base refinement as well. 

## Special output and analysis {#user_speciallog}

The default log-file may be altered, for which you need to code up  
a **usr_print_log** subroutine. For parallel execution,
this invariably means the use of MPI constructs, so you should copy in the
default version **printlog_default** from _src/amrvacio/mod_input_output.t_ 
and then study it, and modify accordingly.

A pointed **usr_aux_output** is extremely handy to compute variables from the
actually computed conserved variables, that can then be visualized directly.
It is only used in combination with the conversion subroutines. e.g., one may
here compute current density components using the actual code discretizations
for computing a curl, and then visualize those with any of the visualization
tools applicable. You then also need to specify a name for each variable, in
a pointed **usr_add_aux_names**.

A pointed **usr_process_grid** and a pointed **usr_process_global** are 
subroutines called in each iteration, which allow to compute auxiliary
variables which happen to be non-local (like div v), and are in no way used
for flux computations. As auxiliaries, they are also not advanced. This
functionality was added to allow e.g. for separate particle treatments using
stochastic differential equations, where the particle dynamics is only relying
on local compression values.
