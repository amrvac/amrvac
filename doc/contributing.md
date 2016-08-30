# Contributing to MPI-AMRVAC

[TOC]

# Introduction {#contrib-intro}

This page describes how you can contribute code and documentation to MPI-AMRVAC.

# Style guide {#contrib-style}

To steadily improve MPI-AMRVAC, we ask that new code contributions take into
account our [style guide](code_style_guide.md).

# git workflow {#contrib-git}

Please work in your own branch, and only include your changes in the master
branch when they are ready and tested. Suppose you have a branch called
`new-visualization-format` which you want to include in `master`, then a simple
way to include your changes would be to:

    # Checkout the master branch
    git checkout master

    # Fetch and merge with remote changes
    git pull

    # Merge with your branch
    git merge new-visualization-format

    # Perhaps resolve conflicts, and write a clear message describing what your
    # changes are.

    # Push your changes back to master
    git push

# Contributing documentation {#contrib-doc}

The [documentation](documentation.md) page explains how to write documentation
for MPI-AMRVAC.
