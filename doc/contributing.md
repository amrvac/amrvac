# Contributing to MPI-AMRVAC

[TOC]

# Introduction {#contrib-intro}

This page describes how you can contribute code and documentation to MPI-AMRVAC.

# Style guide {#contrib-style}

To steadily improve MPI-AMRVAC, we ask that new code contributions take into
account our [style guide](code_style_guide.md).

# Working with git {#contrib-git}

Git is a version control system that you can use for code and other 'line-based'
documents, such as LaTeX files. There are many tutorials on git, so it is
probably best to look at a couple of them and pick one you like, for example:

* https://help.github.com/articles/good-resources-for-learning-git-and-github/
* https://www.atlassian.com/git/tutorials/
* https://git-scm.com/docs/gittutorial
* http://rogerdudler.github.io/git-guide/

If you are not familiar with distributed version control systems, it can take
some time to learn how to use git. However, to contribute to MPI-AMRVAC, you
don't need to know some of the more advanced features. It can be as simple as:

    # Cloning the MPI-AMRVAC repository
    git clone https://gitlab.com/mpi-amrvac/amrvac.git
    cd amrvac

    # Create your own branch called 'new_feature'
    git branch new_feature

    # Switch to your branch
    git checkout new_feature

    # Work on file_a and file_b
    ...

    # Tell git to track the changes in file_a and file_b
    git add file_a file_b

    # Create a commit of your changes
    git commit

At this point, you can ask for access to the MPI-AMRVAC Gitlab repository so
that you can `push` your changes to it. You can also clone the repository and
file a merge request.

# Contributing documentation {#contrib-doc}

The [documentation](documentation.md) page explains how to write documentation
for MPI-AMRVAC.
