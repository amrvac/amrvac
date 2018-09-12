# Command line parameters

AMRVAC can be invoked with command line parameters, for example to do a restart
run or to convert `.dat` files. The documentation for the command line
parameters is available through:

    ./amrvac --help

# Examples

## Using multiple par files

    mpirun -np 4 ./amrvac -i file1.par file2.par

Note that settings from later par files override earlier ones. Par files can
contain only the namelists to which they make changes.

## Resume a previous simulation (which has not finished yet)

    mpirun -np 4 ./amrvac -i settings.par -resume

This will automatically look for the last snapshot (`.dat`) file and resume from
it.

## Restart a previous simulation

    mpirun -np 4 ./amrvac -i settings.par -if output/simulation0010.dat
