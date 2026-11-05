title: About the AGILE documentation

# Generating the documentation {: #doc-gen }
API documentation for AGILE is generated from the source code using
[FORD](https://forddocs.readthedocs.io/) (the Fortran Documenter). Unlike
upstream MPI-AMRVAC, this fork no longer uses Doxygen.

To generate it locally:

```bash
uv sync --extra docs
source .venv/bin/activate
doc/ford/build_docs.sh
```

This first expands `src/*.fpp` with `fypp` (FORD cannot parse fypp directives
directly), once per physics module (`doc/ford/docgen_{mhd,hd,ffhd,srhd,srmhd}.par`),
then runs FORD on the result. Open the output e.g. with

```bash
firefox doc/ford/html/index.html
```

Because AGILE's physics modules (`hd`, `mhd`, `ffhd`, `srhd`, `srmhd`) are
mutually-exclusive fypp branches, a single fypp pass with one `PHYS` value
only shows that physics module's variant of any file whose content depends
on it. `merge_physics_variants.py` auto-discovers (by comparison) which
files actually differ by physics module and keeps those as
`_hd`/`_mhd`/`_ffhd`/`_srhd`/`_srmhd`-suffixed variants; everything else is kept
once. See `doc/ford/project.md` and `doc/ford/build_docs.sh` for details.

# How to write documentation {: #doc-howto }
FORD understands the same `!>`/`!<` Doxygen-style comment syntax already
used throughout AGILE, so existing source comments don't need to change.
For the full set of features (note/warning/todo admonition blocks,
cross-references, LaTeX math, etc.) see the
[FORD user guide](https://forddocs.readthedocs.io/en/latest/user_guide/writing_documentation.html).

## Documenting source code {: #doc-src }
You can write documentation comments almost in the same way as regular
comments, using the following syntax:

```fortran
    ! The number of iterations (normal comment, ignored by FORD)
    integer :: x

    !> The number of iterations (variant 1)
    integer :: bum_its

    integer :: x !< The number of iterations (variant 2)
```

Note that `!>` describes the next statement, whereas `!<` describes the previous statement.
Multi-line comments can be formed in the following way:

```fortran
    !> a long line
    !> of text
    !> it is really long

    !> a long line
    !! of text
    !! it is really long
```

You can document variables, functions, subroutines, modules, types and arguments.
Here are some examples to demonstrate the syntax:

```fortran
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
```

## Documentation in markdown files {: #doc-md }
This page and the rest of `doc/*.md` are folded into the generated site as
plain FORD pages, via the `page_dir` setting -- see the
[FORD page documentation](https://forddocs.readthedocs.io/en/latest/user_guide/writing_documentation.html#writing-pages).
`doc/index.md` is the required root page; every other `doc/*.md` file is a
flat subpage of it. A few conventions carried over from the original
Doxygen-authored files:

* Every page needs a `title: ...` metadata line as its very first line.
* Section anchors use Python-Markdown's `attr_list` syntax,
  `## Heading {: #my-label }` (note the colon and spaces -- Doxygen's
  `{#my-label}` form isn't recognized), and are linked to with ordinary
  markdown links: `[text](#my-label)` on the same page, or
  `[text](otherpage.html#my-label)` from another page.
* Links to another page use its `.html` name, not `.md`
  (`[Getting started](getting_started.html)`, not `getting_started.md`).
* Inline LaTeX uses `\( ... \)`, and display equations use `\[ ... \]` or
  `$$ ... $$` (`mdx_math`'s defaults) -- not Doxygen's `\f$ ... \f$`,
  `\f[ ... \f]`, or `\f{env}{ ... \f}`.
