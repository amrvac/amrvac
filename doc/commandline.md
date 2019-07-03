# Command line parameters

AMRVAC can be invoked with command line parameters, for example to do a restart
run or to convert `.dat` files. The documentation for the command line
parameters is available through:

    ./amrvac --help

which currently prints:

    Usage example:
    mpirun -np 4 ./amrvac -i file.par [file2.par ...]
             (later .par files override earlier ones)

    Optional arguments:
    -convert             Convert snapshot files
    -if file0001.dat     Use this snapshot to restart from
                         (you can modify e.g. output names)
    -resume              Automatically resume previous run
                         (useful for long runs on HPC systems)
    -snapshotnext N      Manual index for next snapshot
    -slicenext N         Manual index for next slice output
    -collapsenext N      Manual index for next collapsed output

# Examples

## Using multiple par files

    mpirun -np 4 ./amrvac -i file1.par file2.par

Note that settings from later par files override earlier ones. There is one
exception to this: `base_filename`, for which the values are concatenated. Par
files do not need to contain all the namelists (see @ref par.md).

## Resume a previous simulation (which has not finished yet)

    mpirun -np 4 ./amrvac -i amrvac.par -resume

This will automatically look for the last available snapshot (`.dat`) file and
resume from it. Note that the indices for the next output files are
automatically set. If the simulation was ended normally reaching wall time limit
(see parameter `wall_time_max` to reach time limit of jobs on HPC systems), the
resumed run will keep the previous output (time) cadence and overwrite the 
snapshot file from which the run resumed, so that the whole data series keep
the same pace as if it was never stopped and resumed.

## Restart a previous simulation from a specific snapshot

    mpirun -np 4 ./amrvac -i amrvac.par -if output/sim0010.dat

To start the new output from a different index:

    mpirun -np 4 ./amrvac -i amrvac.par -if output/sim0010.dat -snapshotnext 100

## Convert a .dat file to another format

    mpirun -np 4 ./amrvac -i amrvac.par -if output/sim0010.dat -convert

See @ref par.md for a list of formats (which should be specified in the `.par` file)
