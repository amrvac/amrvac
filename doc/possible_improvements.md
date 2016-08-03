# Suggestions and improvements

The suggestions and/or improvements below have been divided in major and minor
ones. This division does not reflect their importance, but only how much work is
involved. Please change/update this list according to your own ideas. It is
probably good to put your name next to your suggestions (as I have done below).

# Major improvements or suggestions

## Documentation

**Jannis**: The source code currently has very little documentation. I've already started to
use Doxygen to document the source code extensively (as was done for the python
tools). The work comes down to:

* Bring the manual (the html files have been converted to markdown) up to date
* Document what a subroutine/function does (if it is non-trivial)
* Document non-trivial parts of code / algorithms
* Document the important variables and parameters
* Document the contents and purpose of modules
* Document important examples

Doxygen can currently not parse those files that were only intended to be
included (because they are not valid Fortran by themselves). This conversion is
discussed elsewhere on this page.

## Automated tests

**Jannis**: We should have some basic tests in place to reduce the number of regressions
(bugs) that we introduce when changing the code. After a change, these tests
should compile and run, and their output should match previously checked
*correct* values.

As of now, I'm not sure what kinds of tests we should include, and whether we
already have them. Some ideas:
  * Perhaps we can run some programs through valgrind to detect memory problems,
    although I'm not sure how well this works with MPI
  * We could include a number of basic tests; for example to individually test
    the mesh refinement procedure, the time integration procedures, the
    convergence of the hyperbolic solvers etc.
  * Run all tests in *debuging* mode, i.e., compiled with all relevant
    error-checking and warnings enabled. For example, with `gfortran` you can
    use finit-real=snan to initialize all reals with a signalling `NaN`.

## HPC aspects

**Jannis**: A number of things that caught my eye thus far are listed below

### Fixed number of grids and levels

There seems to be a fixed upper bound on the number of grids (`ngridshi`) and
levels (`nlevelshi`). There are a number of problems with such fixed bounds:

  * You have to manually set them to a sensible number for your application
  * One cannot easily change the grid/domain size (which affects the mesh
    adaptivity) without also changing them
  * There is currently a lot of code like this:

        do igrid = 1, ngridshi

        allocate(recvstatus(MPI_STATUS_SIZE,ngridshi))

    which run slower when you increase `ngridshi`, regardless of the number of
    grids in use.

### Global grid structure

The grid topology/connectivity is currently known on all processors. Would it
make sense to use a distributed approach, where each processor only knows about
its own grid and its connection to other processors? Then you have to implement
some sort of parallel sort algorithm for the load balancing, I believe reading
about this in a paper on the Dendro code.

### Hybrid OpenMP/MPI

Scattered throughout the code there are OpenMP statements. Is the potential
performance gain worth it? And more importantly, do we have the debugging and
maintenance manpower to handle hybrid MPI/OpenMP?

I've tried to look up some literature on hybrid MPI/OpenMP vs pure MPI
performance, and found the following:

[link](http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1302919)
*hybrid models are not efficiently supported by existing MPI implementations,
resulting to an imbalanced message passing that is performed solely by the
master thread.*

[link](http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4912964)
*Mismatch problems, i.e. the unsuitability of current hybrid hardware for
running highly parallel workloads, are often hard to solve, let alone in a
portable way. However, the potential gains in scalability and absolute
performance may be worth the significant coding effort.*

[link](https://goparallel.sourceforge.net/openmp-vs-mpi-better/)
*The paper indicated that there can be significant improvements with a hybrid
approach but found that the older memory model used in OpenMP reduces
performance and they have devised extensions to OpenMP to add “locality”.*

[link](https://computation.llnl.gov/casc/people/chow/pubs/hpaper.pdf)
*My summary: hybrid not advantageous now, may change in future*

[link](https://www.epcc.ed.ac.uk/sites/default/files/PDF/mixedMode3.pdf)
*Our results indicate that performance portability of mixed OpenMP/MPI programs
is likely to be a significant issue on current systems: different
implementations have widely differing performance characteristics*

Therefore I would suggest we focus on a pure MPI approach for the time being.

## Division of the code into modules

**Jannis**: I suggest we convert the physics modules to actual Fortran modules.
They could each define the variables and parameters that they need, and come
with a custom initialization routine.

The use of modules would have several advantages:
  * You can more easily use multiple modules
  * The modules can use/extend each other
  * Constants, variables and parameters can be localized and documented in
 modules, reducing the number of global variables
  * Modules are standard Fortran, so Doxygen and other tools can parse them

A change to modules could significantly affect how users generate their code
(although we could keep it backwards-compatible). We should therefore carefully
think about the (future) structure of MPI-AMRVAC. Some of the points that need
to be addressed are:

  * Do we want to give users the flexibility to combine physics modules? For
    example, when I add a multigrid solver, that could be used for various
    purposes. Perhaps we can define a *hierarchy of modules*, so that higher
    ones 'know' how to include lower ones.
  * How much code do we want to automatically generate with the setup script?
    Now it does quite a lot for the user.
  * Do we want to give different names to the routines in the modules and use
    e.g., function pointers?

When we *modernize* the code to use modules, we could also think about
(re)defining a couple of basic data structures. For example, I think a single
mesh block would seem a natural data structure in MPI-AMRVAC.

## Code standards

**Jannis**: My goal is that new modelers can learn best practices from the
MPI-AMRVAC code. Perhaps we can bundle a style guide with MPI-AMRVAC and try to
adhere to it when we add or change code. Here are already some links to already
existing Fortran style guides and best practices:

* http://www.fortran90.org/src/best-practices.html
* http://www.fortran.com/Fortran_Style.pdf
* http://research.metoffice.gov.uk/research/nwp/numerical/fortran90/f90_standards.html
* https://github.com/Fortran-FOSS-Programmers/Best_Practices

We could also think about more general coding practices (not Fortan-specific), for example:

  * Clear purpose. For example, it should be clear what the physics modules can
    be used for, and what their limitations are.
  * Robustness: the code should not fail on valid input, and should return
    meaningful error messages. Crashes without error messages should be
    prevented.
  * Maintainability: code should be as easy to maintain as possible, which means
    code should be simple (or as simple as possible), testable and documented.

## Limit the usage of the preprocessor (Jannis)

The advantages of the **VACPP** preprocessor are outlined in
[this paper](http://www-personal.umich.edu/~gtoth/Papers/lasy.html):

* Compact code
* No separate 2D and 3D versions (which can lead to bugs)
* Efficiency (no complicated abstractions)

However, in my view, extensive usage of a preprocessor also has a major
drawback: **the code is more complicated to understand**. The project's entry
barrier is also increased, because knowing the programming language (in this
case Fortran) is not enough -- one also has to learn the preprocessor language.

Take for example the following lines and their expansion in 3D:

```
    pflux(iside,^D,igrid)%flux(1^D%:^DD&,1:nwflux) = &
                      fC(ixMhi^D^D%ixM^T,1:nwflux,^D)

    pflux(iside,1,2,3,igrid)%flux(1,:,:,:,1,:,:,:,1,1:nwflux) = fC(ixMhi1,&
       ixMlo2:ixMhi2,ixMlo3:ixMhi3,ixMlo1:ixMhi1,ixMhi2,ixMlo3:ixMhi3,&
       ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMhi3,1:nwflux,1,2,3)
```

```
    itag=4**^ND*(ineighbor-1)+{inc^DD*4**(^DD-1)+}

    itag=4**3*(ineighbor-1)+inc1*4**(1-1)+inc1*4**(1-1)+inc1*4**(1&
       -1),itag=4**3*(ineighbor-1)+inc2*4**(2-1)+inc2*4**(2-1)+inc2*4**(2&
       -1),itag=4**3*(ineighbor-1)+inc3*4**(3-1)+inc3*4**(3-1)+inc3*4**(3-1)
```

For the next example the surrounding braces have been omitted because they broke
the markdown formatting:

```
    do ic^DDB=1+int((1-i^DDB)/2),2-int((1+i^DDB)/2)
    inc^DDB=2*i^DDB+ic^DDB

    do ic3=1+int((1-i3)/2),ic2=1+int((1-i2)/2),ic1=1&
       +int((1-i1)/2),2-int((1+i3,1+i2,1+i1)/2)
    inc3=2*i3+ic3,inc2=2*i2+ic2,inc1=2*i1+ic1
    do ic3=1+int((1-i3)/2),ic2=1+int((1-i2)/2),ic1=1&
       +int((1-i1)/2),2-int((1+i3,1+i2,1+i1)/2)
    inc3=2*i3+ic3,inc2=2*i2+ic2,inc1=2*i1+ic1
    do ic3=1+int((1-i3)/2),ic2=1+int((1-i2)/2),ic1=1&
       +int((1-i1)/2),2-int((1+i3,1+i2,1+i1)/2)
    inc3=2*i3+ic3,inc2=2*i2+ic2,inc1=2*i1+ic1
```

For me, such expressions are quite difficult to understand. Note also the lack
of spacing and indentation in the expanded versions.

I would like to try is to minimize the usage of *complicated preprocessor
patterns* in new code. (Of course, *complicated* is subjective here.) With human
intelligence, many complicated patterns can probably be avoided. In cases where
the power of the preprocessor is really needed, I suggest adding a comment
explaining what is done.

An important question is also whether we aim at general-dimension code or only
at 2D and 3D code, which could simplify matters. Some further minor points:

* Can we use min(^ND), max(^ND) instead of min1, min2, min3 etc.?
* Do we always need #ifdef, or could we just use flags at negligible performance
  penalty? (However, without other changes, this could lead to undefined
  routines)

## Output file formats

**Jannis**: In my experience, the VTK unstructured file format does not scale
well to large 3D simulations, see
[this link](http://scicomp.stackexchange.com/questions/23657/visualization-of-quadtree-octree-grids).
With *large* I mean for example ten million cells or more. Has this been a
problem for MPI-AMRVAC users already? In any case, the VTK unstructured format
is much more general than we need for our hyperoctree (or orthtree) meshes.

Somebody could look into alternatives, which will probably require collaboration
with Paraview/VTK developers.

# Minor improvements

## Commit messages, branches

**Jannis**: How do we organize our work? Do we use our own branches on the central
repository, which we at appropriate times merge into `master`?

And how do we communicate the changes we make to the developers and users of
MPI-AMRVAC? For the developers it makes sense to (at least) write informative
commit messages, so that we can see what we are working on and why somebody has
changed something.

For the user we could activate the mailing list again? We can also list
important changes on a webpage. There is already a changelog file included in
Doxygen.

## A code license

**Jannis**: As of now, it's not clear to me what the license of MPI-AMRVAC is, which could
be an issue for external users. The longer we wait, the harder it will be to get
all the contributors to agree on a license. Some suggestions:

  * GNU copyleft licenses: GPLv2 or GPLv3
  * MIT license (very liberal)
  * BSD license (have to choose clause, non-copyleft)

## Syntax of input files

**Jannis**: Perhaps we could switch to free-form input files with documentation included,
instead of the namelists that are currently used. I have some experience with
this, see https://github.com/jannisteunissen/fortran_config

## Mailing list

**Jannis**: The mailing list seems pretty inactive (I could find only 3 messages).
Do we want to put some effort into this? As of now, the mailing list is not
accessible for outside people it seems.

  * Pros: people receive news automatically, we can summarize the important
    things for the users
  * Cons: people hardly read mail, new people will not get old messages, what
    about the people not on the list?

## Optionally include boundary ghost cells in output

Suggested by **Utzi Utz**: Since numerical or physical 'problems' often start
near a boundary, it would be good to optionally include the boundary ghost cells
in the output, for debugging purposes.

## More flexible initialization

Suggested by **Utzi Utz**: It would be nice if the initialization was more
flexible, so that you could for example define one variable first, and then
compute other variables that depend on it.

**Jannis**: Perhaps this is already possible, but it would be good to think
about where we want to place the 'boundary' between modifying a user code and
modifying MPI-AMRVAC itself.
