Work in progress...

This folder contains the migrated make files.

Inside the `make` directory there are a set of make files that should be sourced in numerical order.

- prelude: common definitions
- find-agile: make sure that `$(amrvac)` is defined from the `AMRVAC_DIR` variable.
- select-arch: the commandline should get a `arch=` argument, or architecture defaults to `gnu`. If there was a previous build, use the same value.
- read-config: there should be a `config.mk` located from where the make command is run.
- build-dir: prepares the build directory and symlinks to the latest build
- dependencies: runs the Python script `fortdepend` which scans dependencies and turns those into a `dependencies.mk` file
- fypp: configures fypp
- compile: defines compile rules
- link: defines linking rule(s)

