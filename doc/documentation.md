# About the MPI-AMRVAC documentation

[TOC]

# Generating the documentation {#doc-gen}

The documentation of MPI-AMRVAC is generated using
[Doxygen](http://doxygen.org/). On every commit to the (owner-protected) master branch, documentation on http://amrvac.org
is automatically updated. You can also generate 
documentation locally if you have a recent version of Doxygen installed. The
latest binaries can be downloaded
[here](http://www.stack.nl/~dimitri/doxygen/download.html).

To locally generate the documentation, execute

    doxygen

in the documentation folder (`doc`). You can then open the documentation in
`html` folder with for example

    firefox html/index.html

# How to write documentation {#doc-howto}

Below some of the Doxygen basics are described, but for more information see the
online [manual](http://www.stack.nl/~dimitri/doxygen/manual/index.html). Having
a quick look at the documentation already present in MPI-AMRVAC (both in the
source code and in the `doc` folder) will also help get you started.

## Documenting source code {#doc-src}

You can write Doxygen comments almost in the same way as regular comments, using
the following syntax:

    ! The number of iterations (normal comment, ignored by Doxygen)
    integer :: x

    !> The number of iterations (Doxygen variant 1)
    integer :: bum_its

    integer :: x !< The number of iterations (Doxygen variant 2)

Note that `!>` describes the next statement, whereas `!<` describes the previous statement.
Multi-line comments can be formed in the following way:

    !> a long line
    !> of text
    !> it is really long

    !> a long line
    !! of text
    !! it is really long

You can document variables, functions, subroutines, modules, types and arguments.
Here are some examples to demonstrate the syntax:

    !> Compute the square of x
    subroutine square(x, x2)
        real, intent(in)  :: x  !< The number we will square
        real, intent(out) :: x2 !< The square of x

        ! This comment will not appear in the documentation
        x2 = x**2
    end subroutine square

    !> A module that contains nothing
    !>
    !> A longer description for the module that does nothing,
    !> although it seems hard to be very elaborate.
    module meaning_of_life
    end module meaning_of_life

    !> A 2D point coordinate
    type coordinate
        real :: x !< The x-coordinate
        real :: y !< The y-coordinate
    end type coordinate

Besides documenting your source code, you can also include other information
using special commands, preceded by a backslash (\). Examples are: **todo** to
list TODO items, **test** to describe test cases and **page** to create a
separate documentation page. The full list of commands is available
[here](http://www.stack.nl/~dimitri/doxygen/manual/commands.html).

## Documentation in markdown files {#doc-md}

It is also possible to write separate markdown files that will show up in the
generated documentation. In these markdown files Doxygen will automatically
generate links for names it recognizes. Here is a short example:

    # The title that doxygen will use for your page

    [TOC]

    # A section {#label-1}

    Please read the [subsection](@ref label-2).

    This is the general link syntax: [link name](link).

    ## A subsection {#label-2}

    You can also place LateX equations in a file, like this:

\f$ f(x) = \sin(x^2) \f$.

Note that in doxygen 1.8.11 you have to define section labels for the table of
contents to work, which appear at the location of the special `[TOC]` command.
