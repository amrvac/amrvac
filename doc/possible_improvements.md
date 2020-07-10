# Suggestions and improvements

The suggestions and/or improvements below are collecting some views.
Please change/update this list according to your own ideas. 

# Improvements or suggestions

### Hybrid OpenMP/MPI

Scattered throughout the code there are OpenMP statements. Is the potential
performance gain worth it? And more importantly, do we have the debugging and
maintenance manpower to handle hybrid MPI/OpenMP?

Some literature on hybrid MPI/OpenMP vs pure MPI
performance:

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
performance.*

[link](https://computation.llnl.gov/casc/people/chow/pubs/hpaper.pdf)
*summary: hybrid not advantageous now, may change in future*

[link](https://www.epcc.ed.ac.uk/sites/default/files/PDF/mixedMode3.pdf)
*results indicate that performance portability of mixed OpenMP/MPI programs
is likely to be a significant issue on current systems: different
implementations have widely differing performance characteristics*


## Code standards

New modelers can learn best practices from the
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

## Limit the usage of the preprocessor 

The advantages of the **VACPP** preprocessor are outlined in
[this paper](http://www-personal.umich.edu/~gtoth/Papers/lasy.html):

* Compact code
* No separate 2D and 3D versions (which can lead to bugs)
* Efficiency (no complicated abstractions)

However, extensive usage of a preprocessor also has a major
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

Such expressions are quite difficult to understand. Note also the lack
of spacing and indentation in the expanded versions.

Therefore, we minimize the usage of *complicated preprocessor
patterns* in new code. (Of course, *complicated* is subjective here.) With human
intelligence, many complicated patterns can probably be avoided. In cases where
the power of the preprocessor is really needed, we suggest adding a comment
explaining what is done.

### Current usage of preprocessor

The following command sequence will show you the current usage of the preprocessor flags.
   
    cd src

    # grep options: -E extended regex, -o only matching, -h no file names
    # sort & uniq -c create a histogram
    grep -E -o -h "\^[A-Z]+" *.t */*.t -R | sort | uniq -c | sort -nr

## Output file formats

The VTK unstructured file format does not scale
well to large 3D simulations, see
[this link](http://scicomp.stackexchange.com/questions/23657/visualization-of-quadtree-octree-grids).
With *large* we mean for example ten million cells or more. 
In any case, the VTK unstructured format
is much more general than we need for our hyperoctree (or orthtree) meshes.

Somebody could look into alternatives, which will probably require collaboration
with Paraview/VTK developers.
