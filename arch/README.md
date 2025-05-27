Makefiles in this directory implement compiler and target specific options. These are sourced depending on the chosen compilation target.

- `compile` set to compiler command
- `f90_flags` common arguments between compiler and linker.
- `compile_flags` arguments added to compile command, defaults to `-c $(f90_flags)`.
- `link` set to linker command, defaults to compile command
- `link_flags` arguments added to link command, defaults to `$(f90_flags)`.
- `enabled` collects arguments for generating a unique hash.
