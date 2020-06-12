# Contributing to MPI-AMRVAC and its documentation

[TOC]

# Introduction {#contrib-intro}

This page describes how you can contribute code and documentation to MPI-AMRVAC.

# Testing your changes {#contrib-testing}

Having made some changes, the first thing you can check is whether the AMRVAC library still compiles. To compile the 1D, 2D and 3D version of the library in parallel with 8 jobs, you can run:

    cd lib
    make -j 8

If you have added new files, you might have to update dependencies in the makefiles, see @ref addmodule.md.

It is important to test your changes by running MPI-AMRVAC's test suite. This is done as
follows:

    cd tests
    make all

This will run programs in 1D, 2D and 3D and compare their output to previously
stored results, and report errors when differences are larger than some
threshold. Alternatively, you can run tests for individual physics modules:

    cd tests
    make rho
    make hd
    ...


## Adding tests {#contrib-adding-tests}

To add new tests, you need:

* A `mod_usr.t` or `mod_usr.f` file
* A suitable `.par` file so that your simulation runs quickly (in say 1 to 10 seconds)
* A file `test.make` following the template below

These files should be placed in the corresponding physics folder, e.g. `tests/hd`
for hydro problems. Assuming your `.par` file is called `my_test.par`, it should contain the following options in its filelist (see @ref par.md):

    base_filename='my_test'
    typefilelog='regression_test'

You can then use the following template for the `test.make` file:

    # Template for test.make files

    # Name of .par file you want to use
    PAR_FILE := my_test.par

    # Base file name as set in the .par file
    BASE_NAME := my_test

    # Change to -d=2 or -d=3 for 2D/3D problems
    SETUP_FLAGS := -d=1

    # This directory contains combinations of numerical schemes that you can use
    SCHEME_DIR := ../../schemes

    # Which of the schemes you want to use for your test
    SCHEMES := 2step_tvdlf_mm 3step_hll_cada 3step_hll_vl 3step_hll_ppm

    # Define the targets for the makefile, which are the log files produced by your
    # program.
    TESTS := $(SCHEMES:%=$(BASE_NAME)_%.log)

    # Include rules on how to compile and run the tests
    include ../../test_rules.make

    # Generate dependency rules for the tests, which are used to run them
    $(foreach s, $(SCHEMES),\
        $(eval $(s:%=$(BASE_NAME)_%.log): $(PAR_FILE) $(SCHEME_DIR)/$(s).par))

Of course, you can change `my_test` in the above files to something more
meaningful, e.g. `hd_dust_2d`. Now you can compile and run your tests with:

    make -f test.make

They should all fail the first time, but there should now be log files of the
"regression_test" type, which contain the volume averages of the variables and
the variables squared over time. If you are happy with the results, store them
in a folder called `correct_output`:

    mkdir correct_output
    cp *.log correct_output/

Now `make -f test.make` should pass. You can then add your test folder to the
corresponding directory in `tests/Makefile` file, so that the test is
automatically performed. Finally, commit and push your changes.

# Style guide {#contrib-style}

To steadily improve MPI-AMRVAC, it would be good if new contributions take into
account our [style guide](code_style_guide.md). Most importantly:

* Choose meaningful names for variables, functions etc.
* For non-trivial blocks of code or routines: describe what they are supposed to do

# Working with git {#contrib-git}

Git is a version control system that you can use for code and other 'line-based'
documents, such as LaTeX files. There are many tutorials on git, so it is
probably best to look at a couple of them and pick one you like, for example:

* https://help.github.com/articles/good-resources-for-learning-git-and-github/
* https://www.atlassian.com/git/tutorials/
* https://git-scm.com/docs/gittutorial
* http://rogerdudler.github.io/git-guide/

You can ask for access to the MPI-AMRVAC Github repository so that you can
`push` your changes to it. You can also clone the repository and/or file a merge
request.

# Contributing documentation {#contrib-doc}

The [documentation](documentation.md) page explains how to write documentation
for MPI-AMRVAC.
