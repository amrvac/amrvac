# A list of possible improvements for MPI-AMRVAC

The improvements below have been divided in major and minor ones. This division
does not reflect their importance, but only how much work is involved.

# Major improvements

## Documentation

## Automated tests

## Modern code standards

indentation, namelists

## Limit the usage of the preprocessor (Jannis)

The advantages of the **VACPP**preprocessor are outlined in
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

```
    ! The braces have been omitted because they break the markdown formatting
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

### Suggested improvement

I would like to try to minimize the usage of *complicated preprocessor
patterns*, although *complicated* is of course very subjective here. My belief
is that with some human intelligence, most of the complicated patterns can
probably be avoided entirely. An important question is also whether we aim at
general-dimension code or only at 2D and 3D code, which could perhaps also
simplify matters.

Some minor points:
* Can we use min(^ND), max(^ND) instead of min1, min2, min3 etc.?
* Do we always need #ifdef, or could we just use flags at negligible performance penalty? (However, without other changes, this could lead to undefined routines)

# Minor improvements

## Commit messages, branches

## A code license
