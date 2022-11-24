# Code style guide


[TOC]

# Introduction {#style-intro}

Below we present some guidelines for new code contributed to MPI-AMRVAC. The
goal of these guidelines is to make it easier to understand, modify and maintain
MPI-AMRVAC. Although some of the guidelines are simply a matter of preference
and convention, sticking to these preferences and conventions is still valuable.
The guide below has been inspired by a number of other Fortran style guides:
[1](http://www.fortran90.org/src/best-practices.html)
[2](http://www.fortran.com/Fortran_Style.pdf)
[3](http://research.metoffice.gov.uk/research/nwp/numerical/fortran90/f90_standards.html)
[4](https://github.com/Fortran-FOSS-Programmers/Best_Practices)

# Standards and language features {#style-standard}

Code should conform to a reasonably modern Fortran standard (e.g., Fortran 90,
95, 2003). Only make use of the latest language features if you really need
them, since this can cause problems for users with different compilers.

Do not use
[deprecated features](http://fortranwiki.org/fortran/show/Modernizing+Old+Fortran):

Deprecated feature | solution
---|---
common block | use a module
implicit typing | explicit typing
statement functions | define a real function
`double precision :: x` | `real(kind) :: x`
`real*N :: x` | `real(kind) :: x`
`integer*N :: x` | `integer(kind) :: x`
`character*N :: x` | `character(len=N) :: x`

# Naming conventions {#style-naming}

**Use self-explanatory names when possible**. Only variables that are locally
defined and have a simple meaning, such as loop indices, can have short names
such as `i` or `ix`. The farther away the definition of a variable is, the more
important it becomes that is has a self-explanatory name.

Here are some examples, inspired by looking at the current version of MPI-AMRVAC:

Description | Too short | Better
---|---|---
Time variable | `t` | `time`
Iteration counter | `it` | `iteration`
Maximum refinement level | `mxnest` | `max_refinement_level`
Subroutine to update AMR grid | `resettree` | `update_amr_grid`

A list of further naming guidelines:
* Use lowercase for Fortran constructs (`do`, `if` etc.)
* Use lowercase for the names you define (also of modules)
* Separate words with underscores (`read_snapshot` instead of `readSnapshot`)
* Give source files meaningful names
* Don't use Fortran keywords as names
* Avoid ambiguous names

# Code style conventions {#style-conventions}

* Use `>=`, `<=`, `/=` etc. instead of `.gt.`, `.lt.`, `.ne.`
* Include spaces around operators to improve readability: 
  ```{fortran}
      i = 5
      i = 5 + j
      i = 5/(a + b) + c
  ```
* Include spaces after commas to visually separate arguments:
  ```{fortran}
      call my_routine(a, b, c)
  ```
* Don't use semicolons to put multiple statements on one line
  ```{fortran}
      i = 1; j = i + 1 ! Bad
  ```
* Always use `::` in variable declarations
  ```{fortran}
      integer :: i
  ```

# Comments {#style-comments}

Write comments such that a reasonably experienced programmer can understand your
code as quickly as possible. In practice, this could mean writing the following
comments:

* **Modules**: describe their functionality
* **Subroutines** and **functions**: describe what they do
* **Difficult to read code**: describe what it does or simplify it
* **Important variables**: describe their meaning

When in doubt, you should probably add a comment: it will help the reader to
confirm his/her understanding. Do not write comments around trivial statements,
such as `i = i + 1`.

The [documentation](documentation.md) page explains how to write Doxygen
comments, which show up in the documentation of MPI-AMRVAC.

# Modules and programs {#style-mod}

A source file should contain either a module or a program. Try to use modules in
the following way:

```{fortran}
    module mod_name

      ! Selectively import when needed in the whole module
      use mod_first, only: ...

      implicit none ! Don't allow implicit types
      private       ! Entities are private by default

      public :: ... ! Explicitly state public entities

    contains ! Place subroutines and functions below

      subroutine test()
         use mod_second, only: ... ! Use local imports
         ...
      end subroutine test

    end module mod_name
```

The selective imports with `use ..., only:` help the reader with figuring out
where things come from. However, when you need a lot of functionality from
another module, it makes more sense to simply include everything with `use
module`.

# Functions and subroutines {#style-proc}

All functions and subroutines should at least have a brief comment describing
their functionality, as stated in @ref style-comments.

## Side-effects

Functions should not have *side-effects*, meaning that they do not change the
state of the program. Functions without such side-effects can be marked with
**pure**. Functions should also not change their arguments. Consider making pure
functions operating on scalars **elemental**, so that you can also use them on
arrays.

The result of a function should not be directly used in a `print` or `write`
statement, unless you are sure the function itself does no IO (input/output).
The reason for this is that recursive IO can lead to hard to debug hangs /
crashes.

Subroutines are ideally also **pure** (without side-effects), but they can change
their arguments.

## Arguments

All arguments should have their **intent** declared, using one of
**intent(in,out,inout)**.

Ideally, functions and subroutines get their information from their arguments.
If that requires a lot of arguments, perhaps the functions and subroutines can
be split in smaller pieces. Sometimes it can also help to define types with the
relevant information. There should not be more than about 4 or 5 arguments.

# Global and module variables {#style-global}

Don't use global variables. Instead place variables inside modules, and include
them where needed. In cases where this is not possible, gather the global
variables in a module and give them a recognizable prefix such as `GLOBAL_`.

To ensure that variables defined in a module are not accidentally changed, you can use the Fortran 2003 `protected` attribute:

```{fortran}
    integer, protected :: my_age
```

Such variables can only be changed from inside the defining module.

# Numbers {#style-num}

## Floating point

The usage of `double precision` is deprecated. Instead use the kind-parameter
`dp` to declare double precision numbers:

```{fortran}
    real(dp) :: x         ! good
    double precision :: x ! deprecated
```


Write floating point constants using the `_dp` suffix:
```{fortran}
    x = x * 1.1_dp ! good
    x = x * 1.1    ! bad: single precision constant
```

Integers are automatically converted to real numbers when they are mixed,
without a loss of precision. However, be careful with expressions such as:

    x = 0.5_dp * 5 / 6   ! Bad: is this (0.5_dp * 5)/6 or 0.5_dp * (5/6)?
    x = (0.5_dp * 5) / 6 ! Good

## Integers

Use the default integer type, unless there is a specific reason to use a larger
type. Be aware of the different ways of converting a real number `x` to an
integer:

* `int(x)`: round towards zero
* `floor(x)`: round down
* `ceiling(x)`: round up
* `mint(x)`: round to nearest integer

## Named constants vs magic numbers

Use named constants instead of *magic numbers*:

```{fortran}
    real(dp), parameter :: pi = 3.14_dp

    x = sin(pi * x)      ! good
    x = sin(3.14_dp * x) ! bad: magic number
```

However, don't use named constants for simple numbers such as `0`, or `0.5`:

```{fortran}
    x = 1.0_dp/pi         ! good
    x = one/pi            ! bad: what is one?
```

# Arrays {#style-array}

## Correct order for loops

Always loop over arrays in the *correct* order, meaning that the loop indices are ordered from right to left (opposite to loops in e.g., C or C++):

```{fortran}
    do j = 1, N
       do i = 1, N
          array(i, j) = 1
       end do
    end do
```

Because Fortran stores arrays in column-major order, the loop then follows the
'natural' memory order of the array. For the same reason, it pays off to think
about the ordering of your arrays, for example:

```{fortran}
    array(:, i, j) = 1 ! Fast
    array(i, j, :) = 1 ! Slower
```

## Allocatable arrays

Use allocatable arrays when:

* Arrays are 'large' (e.g., more than a few kilobytes)
* You would like to have a variable array size
* A pointer-based array is not strictly required

Allocatable arrays are automatically deallocated, so you cannot create memory
leaks. They are also automatically reallocated on assignment (with `gfortran`,
for `ifort` you have to enable this behavior).

## Passing arrays as arguments

Use **assumed-shape** arrays when performance is not critical, for example:

```{fortran}
    subroutine square(values, squares)
       integer, intent(in) :: values(:)
       integer, intent(out) :: squares(:)
       ...
    end subroutine square
```

The advantage of assumed-shape arrays is that they allow for run-time bounds
checking. When performance is critical, you can use explicit-shape arrays:

```{fortran}
    subroutine square(values, squares, n_values)
       integer, intent(in) :: n_values
       integer, intent(in) :: values(n_values)
       integer, intent(out) :: squares(n_values)
       ...
    end subroutine square
```

Another option is to embed the array in a type, and pass the type as argument.
See this
[guide](https://software.intel.com/en-us/articles/fortran-array-data-and-arguments-and-vectorization)
by Intel for more information about the performance differences.

## Array constants

In new code, use the `[...]` syntax to define array constants:

```{fortran}
    x = [1, 2, 3]   ! Fortran 2003 and later
    x = (/1, 2, 3/) ! Old style
```

# Indentation and whitespace {#style-indent}

Indent with spaces, using the following indentations:

Construct | number of spaces
---|---
`program`, `module`, `subroutine`, `function`, `associate` | 2
`type`, `interface` | 3
`do`, `if`, `select case`, `where`, `forall` | 3

Include empty lines to visually separate blocks of code. Include a newline
between functions and subroutines, after and before `do` loops, around `if`
blocks, after the variable declarations, etc.

# Warnings and runtime checks {#style-warning}

Enable all normal warnings during compilation, ideally using different
compilers. Fix all warnings, except those that are incorrect, and those for
which the fix does more harm than good.

Run the code with different inputs, perhaps using a different number of
processors, to detect the most obvious bugs. Enable all sensible run-time checks
to catch errors. Try to use a tool like `valgrind` to detect memory problems.

**Useful Gfortran flags**

Flag | meaning
---|---
-Wall | Enable all standard warnings
-Wextra | Enable additional warnings
-Wimplicit-interface | Warn about procedures with implicit interfaces
-fcheck=all | Enable all run-time checks
-ffpe-trap=... | Enable floating point traps. Recommended: invalid, zero, overflow
-finit-real=snan | Initialize real values with a signalling `NaN`
-finit-integer=-2147483648 | Initialize integers with this value
-g | Include debugging info into the executable
-pg | Enable support for profiling with a tool like `gprof`
-O0 | Disable optimizations (when debugging)
-fimplicit-none | Do not allow implicit types
-std=... | Set the Fortran standard to use (f95, f2003, etc.)


**Useful ifort flags**

Flag                 | meaning
---                  | ---
-traceback           | Generate a backtrace after a run-time error
-ftrapuv             | Initialize stack variables to unlikely values
-warn all            | Enable all standard warnings
-check all           | Enable all run-time checks
-check ...           | Useful are: bounds, uninit
-fpe[0,1,3]          | Lower values enable more floating point exceptions
-g                   | Include debugging info into the executable
-pg                  | Enable support for profiling with a tool like `gprof`
-O0                  | Disable optimizations (when debugging)
-implicitnone        | Do not allow implicit types
-stand ...           | Set the Fortran standard to use (f95, f2003, etc.)
