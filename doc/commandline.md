# Command line parameters

AMRVAC can be invoked with command line parameters, for example:

    amrvac --help

will print a welcome message:

```
    -----------------------------------------------------------------------------
    -----------------------------------------------------------------------------
    |         __  __ ____ ___        _    __  __ ______     ___    ____         |
    |        |  \/  |  _ \_ _|      / \  |  \/  |  _ \ \   / / \  / ___|        |
    |        | |\/| | |_) | |_____ / _ \ | |\/| | |_) \ \ / / _ \| |            |
    |        | |  | |  __/| |_____/ ___ \| |  | |  _ < \ V / ___ \ |___         |
    |        |_|  |_|_|  |___|   /_/   \_\_|  |_|_| \_\ \_/_/   \_\____|        |
    -----------------------------------------------------------------------------
    -----------------------------------------------------------------------------
    calling example:
    ./amrvac -i inifile -restart 100 -convert -if datamr/data
    default inifile is amrvac.par
    Note that inifile parameters overwrite the commandline
    -----------------------------------------------------------------------------
```

With the command line, it is possible to initiate a restart or convert. This
has the advantage that the par file does not need to be modified. To restart a
run, type

    amrvac -restart 10 -if datamr/data

This will restart from the 10 output of the datamr/data file. The next output
number is automatically set to 11. Note that the values assigned in the par-
file will overwrite the command line parameters.

In order to initiate a convert, consider the following example

    amrvac -restart 10 -if datamr/data -convert

With the option `-i inifile` it is also possible to prescribe an alternative
par-file, instead of the default `amrvac.par`. This comes in handy at the
convert stage.

An overview of the current command line parameters is given in the following
table.

Command | Sets variable | default
---|---|---
-i parfile | inifile | amrvac.par
-if filename | filenameini | 'unavailable'
-restart snapshot | snapshotini | -1
-convert | convert | F
-slice | slicenext | 0
-collapse | collapsenext | 0
--help | Prints the help message | F
