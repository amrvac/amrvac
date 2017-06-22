# MPI-AMRVAC

This is the development version of MPI-AMRVAC.

# Getting started

Read the 'Getting Started' page of the Doxygen documentation, also available via
doc/getting_started.md.

# Make a local copy of the repository

    git clone https://gitlab.com/mpi-amrvac/amrvac.git

# Update your local repository

To get the latest version of the branch you are on, perform:

    git pull

To merge the changes of a specific branch, you can use:

    git fetch
    git merge origin/branch_name

where `branch_name` could for example be be `master` or `physics_modules`.

If you have made local changes that you want to share, have a look at
doc/contributing.md.

# Running tests

Tests (which can be used as examples) are available in the `tests/` folder,
where they are sorted by their physics type. You can run many of them
automatically with the following command:

    cd tests
    make

# Generate local documentation

    cd doc
    doxygen

Then open doc/html/index.html with a web browser

