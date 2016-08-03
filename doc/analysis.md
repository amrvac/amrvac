# Writing a custom analysis subroutine

This is a subroutine which can be scheduled freely to derive quantities of
interest from the entire grid.

## Setup in par file

Similar to the other IO-mechanisms, analysis is scheduled within the savelist:

     &savelist;
            itsave(1,5)      = 0
            ditsave(5)       = 1
    /

The above example makes sure write_analysis is called every iteration.

## write_analysis subroutine

The file _$AMRVAC_DIR/src/amrvacio/analysis.t_ contains the subroutine
_write_analysis_. To use it, simply copy this file into the working directory
and modify. This routine is on the lowest abstraction level of _MPI-AMRVAC_
such that the user has the flexibility to loop over the the SFC and freely
perform integrals over the grid-functions, similar to the _printlog_special_.
A special file-unit is reserved for the IO: _unitanalysis_. This should be
used for any IO performed by this routine.

