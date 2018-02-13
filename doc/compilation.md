# Compilation

# Compilation settings

Different predefined rules for compilation can be found in the `arch/` folder, for example:

* `default.arch`: The default settings, using GNU `gfortran`
* `debug.arch`: Debug flags, using GNU `gfortran`
* `intel.arch`: Using Intel `ifort`
* `inteldebug.arch`: Debug flags, using Intel `ifort`

These files set the compiler/linker and compiler flags. You can select the
default settings when you create a makefile, for example to use debug flags:

    $AMRVAC_DIR/setup.pl -d=2 -arch=debug

You can also edit the resulting `makefile` by hand later, or specify the compilation settings on the command line:

    make ARCH=debug

# Adding user libraries

If you make use of external libraries, you can add them to a file `local.make`.
This file will automatically be read (and included) by the AMRVAC `makefile`.
Suppose you want to include a library `libhdf5.so` (or `libhdf5.a`) which is in
a folder called `my/lib/folder`, then you could specify:

    LIBS += hdf5
    LIB_DIRS += my/lib/folder


