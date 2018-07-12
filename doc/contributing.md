# Contributing to MPI-AMRVAC

[TOC]

# Introduction {#contrib-intro}

This page describes how you can contribute code and documentation to MPI-AMRVAC.

# Testing your changes {#contrib-testing}

Having made some changes, the first thing you can check is whether the AMRVAC library still compiles. To compile the 1D, 2D and 3D version of the library in parallel with 8 jobs, you can run:

    cd lib
    make -j 8

If you have added new files, you might have to update dependencies in the makefiles, see @ref addmodule.md.

It is important test your changes by running MPI-AMRVAC's test suite. This is done as
follows:

    cd tests
    make all

This will run programs in 1D, 2D and 3D and compare their output to previously
stored results, and report errors when differences are larger than some
threshold. Alternative, you can run tests for individual physics modules:

    cd tests
    make rho
    make hd
    ...


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
