# Source language

[TOC]

# Introduction {#source_intro}

This document describes the structure and the syntax of the **src/FILENAME.t**
source files, and also provides explanation to the most common expressions. A
brief summary at the end is provided as a quick reference guide.

A more general and possibly more enlightening description of the source language
can be found in a paper on the [LASY Preprocessor](http://www-personal.umich.edu/~gtoth/Papers/lasy.html).

# Structure {#structure}

MPI-AMRVAC does both the initialization as well as the advancing of the
variables in the governing PDEs. The program consists of several modules.
Modules are simply sets of subroutines and functions that belong together and
they are put in a single file.

The main program is in the file **amrvac.t** and all time-advancing actually
happens in **mod_advance.t**. Input and output routines are in the **amrvacio** directory which also contains the postprocess conversion part collected in **convert.t**.
Subroutines likely to be modified by the user are to be collected in the **mod_usr.t** module.
In this module, you must at the very minimum define all global parameters and
provide initial conditions for all (conserved) variables on the grid, by
following the template provided in **mod_usr_template.t**.

Independence of equations is ensured by putting all flow variables in a single
array **w**, and the variables are distinguished by their named indices, e.g.
**w(ix1,ix2,rho_)** would read **rho(ix1,ix2)** in a typical code designed for
a single equation. Here **rho_** is an integer parameter defined in the
**mod_rho_phys.t** module.

The physics modules have several variants, and they are found in the corresponding subdirectories called **rho** for advection, **hd** for hydro (plus dust) and **mhd** for MHD.
All physics modules are available in the library version of AMRVAC, and the user has to active the module of choice for the application in mind.

The spatial discretization of the equations algorithms are spread over the
modules **mod_finite_volume.t**, **mod_finite_difference.t**, **mod_tvd.t** where limiter implementations are obtained from **mod_mp5.t** **mod_ppm.t** **mod_limiter.t**.

The explicit temporal discretizations are in the main advancing module
**mod_advance.t**.

The info on the grid related quantities is all in the **geometry.t** module.

# Syntax {#syntax}

  1. Dimensions and vector components

Since we aim to solve equations independent of the number of dimensions,
arrays with unknown number of dimensions have to be declared and manipulated.
This is achieved by the extensive use of the preprocessor VACPP.

The basic idea of the Loop Annotation Syntax (LASY) is defining loops in the
source text. A first example may be an array with 3 indices:

        a({ix^D|,})   -->     a(ix1,ix2,ix3)

The string **ix^D** between the **{** and **|** characters is repeated NDIM=3
times with the pattern string **^D** replaced by **1**, **2** and **3**, and
the three resulting strings are separated by the string between the **|** and
**}** characters, in this case a single comma. As you may observe, some
special characters are used to help the preprocessor in recognizing the loop.
The loop is enclosed between curly brackets, the pattern consists of a **^**
character followed by uppercase letters (or **%**, or **&amp;**), and the
separator string is preceded by the **|** character. Thus the full syntax of
the loop is

        { text ^PATTERN text ^PATTERN ... | separator }

Fortunately some default values can be introduced to simplify the notation.
The default separator is the comma, thus the above example can be written as
**a({ix^D})**. Furthermore the preprocessor can expand patterns without
enclosing curly brackets. When it sees a pattern, it looks for _bounding
characters_ on both sides. These are one of _comma, space, newline, semicolon_
and _enclosing parentheses_. Since our **^D** pattern is enclosed by a left
and a right parenthesis, we may simply type

        a(ix^D)

which will expand to **a(ix1)**, **a(ix1,ix2)**, and **a(ix1,ix2,ix3)** for
the choices NDIM=1, 2, and 3, respectively.

The curly brackets are used only when the repetition should expand over some
bounding characters. The required separator is very often a single character,
such as in the calculation of sums or products:

        ix^D*    -->   ix1*ix2*ix3

Here the **^D** pattern is enclosed by a new line and a space character, and
the preprocessor checks the last character of the repeated string. If it is
one of comma, space, **+**, **-**, *****, **/**, **:**, **;**, or **\**, it is
taken to be a separator character. The **\** is replaced by a new line. A
simple use may be nested DO loops:

    {do ix^D=1,100 \ }
        a(ix^D)= ix^D*
    {enddo\ }

which translates to:

    do ix1=1,100
    do ix2=1,100
        a(ix1,ix2)= ix1*ix2
    enddo
    enddo

Here NDIM=2 was assumed. Note the need for curly brackets in the first line to
override the space and the comma, and also note the space in **a(ix^D)=
ix^D*** which is needed to bound the **ix^D*** loop. In the final
**{enddo\\ }** line the number of repetitions was assumed to be NDIM, as we
shall see this may not always be the case.

Suppose we want to assign the same value to NDIM variables:

        ix^D=1;     -->     ix1=1;ix2=1;

The semicolon is a bounding character, thus strictly speaking it does not
belong to the **ix^D=1** loop. The preprocessor, however, checks wether the
right bounding character is a semicolon, and if it is, it becomes the
separator. The original trailing semicolon is also preserved, which turns out
to be a useful feature.

Up to this point we used a single pattern **^D** only, which is replaced by
the numbers **1..NDIM**. During the code development it became obvious that
the preprocessor can be used for many other things than repeated indices for
the dimensions. The first new application is the components for vector
variables, which run from 1 to NDIR, and in general 1&lt;=NDIM&lt;=NDIR&lt;=3.
When more than one pattern is found in a loop,
the first one determines the number of repetitions, and only patterns with the
same first letter are replaced by their substitutes.

Look at the following typical case construct for NDIM=2:

    select case(iw)
      {case(mom(^D))
         mom(ixmin^D:ixmax^D)=w(ixmin^D:ixmax^D,mom(^D)) \ }
    end select

    -->

    select case(iw)
      case(mom(1))
         mom(ixmin^D:ixmax^D)=w(ixmin^D:ixmax^D,mom(1))
      case(mom(2))
         mom(ixmin^D:ixmax^D)=w(ixmin^D:ixmax^D,mom(2))
    end select

    -->

    select case(iw)
      case(mom(1))
         mom(ixmin1:ixmax1,ixmin2:ixmax2)=&
            w(ixmin1:ixmax1,ixmin2:ixmax2,mom(1))
      case(mom(2))
         mom(ixmin1:ixmax1,ixmin2:ixmax2)=&
            w(ixmin1:ixmax1,ixmin2:ixmax2,mom(2))
    end select

Notice that the preprocessor first expanded the loop for the first **^D**
pattern and substituted **^D** in the third index of the **w** array but the
other indices with the **^D** patterns were left alone. In the second
iteration the loops with **^D** patterns were expanded. The resulting lines
may become extremely long, thus the preprocessor breaks them into continuation
lines.

####

  2. Limits, segments, and nulpatterns

To further shorten the notation the often used **min^D** and **max^D** strings
were interpreted as a loop of the **^L** pattern (mnemonic: Limits). Yes, a
pattern may expand into other pattern(s) thus producing another loop. To
declare the limits of the array sections one may say

        integer:: ix^L  -->  integer:: ixmin^D,ixmax^D  -->

    integer:: ixmin1,ixmin2,ixmax1,ixmax2

The intermediate **^D** pattern becomes very important in the following
typical example:

        jx^L=ix^L+kr(idim,^D);                                       -->

    jxmin^D=ixmin^D+kr(idim,^D);jxmax^D=ixmax^D+kr(idim,^D);     -->

    jxmin1=ixmin1+kr(idim,1);jxmin2=ixmin2+kr(idim,2);
    jxmax1=ixmax1+kr(idim,1);jxmax2=ixmax2+kr(idim,2);

The purpose is to shift the limits by 1 in the idim direction. The **kr**
array is a Kronecker delta, it is 1 for the diagonal elements where the two
indices are the same, and 0 otherwise. The **j** in the **jx** is the next
letter in the alphabet after the **i**, thus **jx** is the mnemonic for
**ix****+****1**. Similarly **hx** is consistently used for **ix****-****1**.
The semicolon remains the separator for both loop expansions thanks to the
trailing semicolon in the intermediate step. Also note that the parentheses of
**kr(idim,^D)** do not bound the initial loop for **jx^L=...** because they do
not enclose it. Once we get used to the notation it becomes natural to think
of an **ix^L** type loop as a Pascal record or C structure consisting of 2, 4
or 6 integers. Therefore **^LADD**, **^LSUB** and **^LT** are introduced to do
operations and comparisons on them. **^LADD** (mnemonic: add to limits) is
replaced by **-** and **+** to decrease lower and increase upper limits,
**^LSUB** does the opposite by expanding to **+** and **-**, finally **^LT**
has substitutes **&gt;** and **&lt;** ensuring that the limits on the left are
within the limits on the right. Thus we may extend the **ix^L** limits by 1
with

        ix^L=ix^L^LADD1;  -->  ixmin^D=ixmin^D-1;ixmax^D=ixmax^D+1;

    ixmin1=ixmin1-1;ixmin2=ixmin2-1;ixmax1=ixmax1+1;ixmax2=ixmax2+1;

To check if the **ixI^L** input limits are not narrower than the **ix^L**
limits

        if(ixI^L^LTix^L|.or.|.or.)stop -->

    if(ixImin^D>ixmin^D|.or. .or. ixImax^D<ixmax^D|.or.)stop  -->

    if(ixImin1>ixmin1.or.ixImin2>ixmin2 .or. ixImax1<ixmax1.or.&
       ixImax2<ixmax2)stop

where note the repeated use of the **.or.** separator string.

It is possible to form array sections from the **^L** pattern by making use of
the **:** as a separator, but it turns out that introducing the new patterns
**^LIM** (mnemonic: three letter LIMits min and max) and **^S** (mnemonic:
Sections) is a better solution. The typically internally used **^LIM** pattern
expands to **min** and **max**, while **^S** expands to **^LIM1:..^LIMndim:**,
and therefore

        a(ix^S)  -->  a(ix^LIM1:,ix^LIM2:)  -->  a(ixmin1:ixmax1,ixmin2:ixmax2)

As you may suspect the whole exercise of introducing **^LIM** and **^S**
served the purpose of achieving this extremely compact notation for array
sections. In array declarations the **lo** and **hi** absolute limits are used
instead of the **min** and **max** actual limits. The **^LL** pattern
(mnemonic: Lowest/highest Limits) thus expands to **lo^D** and **hi^D**
similarly to **^L**, while **^LLIM** gives simply **lo** and **hi**, and
**^T** (mnemonic: total segment) is used in most array declarations:

        a(ix^T)  -->  a(ix^LLIM1:,ix^LLIM2:)  -->  a(ixlo1:ixhi1,ixlo2:ixhi2)

Now that we have many different patterns with different number of substitutes
it becomes useful and necessary to introduce patterns which are replaced by
nulstrings, but they determine the number and kind of substitutions. The
**^D&amp;** and **^C&amp;** patterns produce NDIM and NDIR repetitions
respectively. (The alternative names ^DLOOP and ^CLOOP can also be used to
avoid syntax problems at the end of a line, where the &amp; means continuation
according to FORTRAN 90.) A trivial application is the enddo-s at the end of
NDIM nested do loops:

        enddo^D&;     -->    enddo;enddo;

Obviously the **^D** pattern would not work here: **enddo^D; --&gt;
enddo1;enddo2;**. For scalar products of vector variables the **^D&amp;** is
used at the head of the loops to tell VACPP that first the **^D** patterns
should be expanded.

####

  3. Breaking the symmetry

The symmetry of the indices is sometimes broken. In the code it is mostly
related to the assumption of axial symmetry, when the first dimension becomes
special. The rest of the dimensions and sections are expanded from the **^DE**
(mnemonic: Dimensions Extra) and **^SE** (Sections Extra) patterns. The radial
distance may be related to the **x** coordinate array by

        r(ix^LIM1:)=x(ix^LIM1:,ixmin^DE,1)   -->

    r(ixmin1:ixmax1)=x(ixmin1:ixmax1,ixmin2,ixmin3,1)

then **r** may be a weight function for an array

        forall(ix= ix^LIM1:) a(ix,ix^SE)=r(ix)*a(ix,ix^SE)      -->

    forall(ix= ixmin1:ixmax1) a(ix,ixmin2:ixmax2)=r(ix)*a(ix,ixmin2:ixmax2)

If NDIM=1 the preprocessor removes the separator, in this case a comma in the
**ix,ix^SE** string, and the **ix^SE** loop is repeated 0 times, thus we get
**forall(ix= ixmin1:ixmax1) a(ix)=r(ix)*a(ix)** as expected.

Another type of symmetry breaking occurs when something is done for all the
indices separately, e.g. the boundary elements of an NDIM dimensional array
are filled up for each boundary separately. The **^D%** pattern is replaced by
**^%1..^%NDIM** patterns. The **^%N** pattern substitutes the text in front of
it in the N-th repetition, and the text behind in the other substitutions. In
general **head^%Ntail --&gt; tail,..,tail,head,tail,..,tail** with **head**
being at the N-th position. The number of repetitions is determined by the
first pattern in **head** and **tail**. Here is an application for boundary
conditions with **ixB^L** enclosing the boundary region and **ixMmin^D** being
the edge of the mesh:

        {case(^D)
       do ix= ixBmin^D,ixBmax^D
          a(ix^D%ixB^S)=a(ixMmin^D^D%ixB^S)
       end do \ }

    --> /loop for ^D/

    case(1)
       do ix= ixBmin1,ixBmax1
          a(ix^%1ixB^S)=a(ixMmin1^%1ixB^S)
       end do
    case(2)
       do ix= ixBmin2,ixBmax2
          a(ix^%2ixB^S)=a(ixMmin2^%2ixB^S)
       end do

    --> /loops for ^S taking ^%1 and ^%2 into account/

    case(1)
       do ix= ixBmin1,ixBmax1
          a(ix,ixB^LIM2:)=a(ixMmin1,ixB^LIM2:)
       end do
    case(2)
       do ix= ixBmin2,ixBmax2
          a(ixB^LIM1:,ix)=a(ixB^LIM1:,ixMmin2)
       end do

    --> /loops for ^LIM/

    case(1)
       do ix= ixBmin1,ixBmax1
          a(ix,ixBmin2:ixBmax2)=a(ixMmin1,ixBmin2:ixBmax2)
       end do
    case(2)
       do ix= ixBmin2,ixBmax2
          a(ixBmin1:ixBmax1,ix)=a(ixBmin1:ixBmax1,ixMmin2)
       end do

####

  4. Special features: dimensional switches

Depending on the actual use of MPI-AMRVAC some subroutines may or may not be
present, so something where you really need to use a 1D construct is best
embedded as

        {^IFONED   (something where x(ix1) occurs e.g.) }

This should be done with care, as the idea of the code is to use it for any
dimensionality. Still, certain subroutines are not used at all in 1D, or they
may be needed only in 2D, and the **^NOONED**, **^IFTWOD** etc. switches are
then used to make a part of the code conditional on the value of the NDIM
parameter.

## Summary

  1. VACPP Patterns

We assume NDIM=2, so you may try this interactively with **vacpp.pl
-d=2**. For a full list of defined patterns read the **&amp;patdef**
definitions in **vacpp.pl**.

        Pattern   Substitutes                              Mnemonic

    * * *

    ^D    --> 1,2                                      Dimensions
    ^DE   --> 2                                        Dimensions Extra
    ^D&   --> ,                                        NDIM repetitions
    ^DLOOP--> ,                                           NDIM repetitions
    ^DD   --> ^D,^D           --> 1,2,1,2              NDIM*NDIM dimensions
    ^D%   --> ^%1,^%2         --> head,tail,tail,head  NDIM*NDIM head-tail matrix
    ^C    --> 1,2,3                                    Components
    ^C&   --> ,,                                       NDIR repetitions
    ^CLOOP--> ,,                                           NDIR repetitions
    ^L    --> min^D,max^D     --> min1,min2,max1,max2  Limits in all D
    ^LIM  --> min,max                                  LIMits min and max
    ^LLIM --> lo,hi                                    Low and high LIMits
    ^LADD --> -,+                                      Add to ^L (extend)
    ^LSUB --> +,-                                      Subtract from ^L (shrink)
    ^LT   --> >,<                                      Less Than (compare ^L)
    ^LL   --> lo^D,hi^D       --> lo1,lo2,hi1,hi2      Lowest/highest Limits
    ^S    --> ^LIM1:,^LIM2:   --> min1:max1,min2:max2  Segments
    ^SE   --> ^LIM2:          --> min2:max2            Segments Extra
    ^T    --> ^LLIM1:,^LLIM2: --> lo1:hi1,lo2:hi2      Total segments

    * * *

####

  3. VACPP Expressions

        Source VACPP notation     Expanded Fortran 90 code

    * * *

    integer:: ix^D        --> integer:: ix1,ix2

    integer:: ix^L        --> integer:: ixmin1,ixmin2,ixmax1,ixmax2

    real:: w(ixG^T,nw)    --> real:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

    call subr(w,ix^L)     --> call subr(w,ixmin1,ixmin2,ixmax1,ixmax2)

    ixI^L=ix^L^LADD1;     --> ixImin1=ixmin1-1;ixImin2=ixmin2-1;
                              ixImax1=ixmax1+1;ixImax2=ixmax2+1;

    jx^L=ix^L+kr(2,^D);   --> jxmin1=ixmin1+kr(2,1);jxmin2=ixmin2+kr(2,2);
                              jxmax1=ixmax1+kr(2,1);jxmax2=ixmax2+kr(2,2);

    {^IFTWOD text}        --> text or '' (depending on ndim being equal to 2)

    {^NOONED text}        --> '' or text (depending on ndim being equal to 1}

    * * *
