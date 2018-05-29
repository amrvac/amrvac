# Writing a custom analysis subroutine

# Introduction

Users can write a subroutine to freely derive (and print or store) quantities of
interest from the entire grid. To activate this functionality, set the following
in your `usr_init` routine:

```{fortran}
subroutine usr_init()
  ...
  usr_write_analysis => my_analysis
  ...
end subroutine usr_init
```

# Setup in par file

Similar to the other IO-mechanisms, analysis is scheduled within the savelist:

     &savelist;
            dtsave_custom = 0.1 ! Save every 0.1 time units
            ditsave_custom = 10 ! Save every 10 iterations
    /

# write_analysis subroutine

The user-defined analysis subroutine takes no arguments. It operates on the
lowest abstraction level of MPI-AMRVAC so that the user has the flexibility to
loop over the the SFC and freely perform integrals over the grid-functions,
similar to the `printlog_special`. A special file-unit is reserved for the IO:
`unitanalysis`. This should be used for any IO performed by this routine.

# Example of analysis routine

```{fortran}
subroutine my_analysis()
  double precision :: tmp, wmax
  integer          :: iigrid, igrid, ierrmpi

  tmp = 0
  do iigrid = 1, igridstail
    igrid = igrids(iigrid)
    tmp = max(tmp, maxval(pw(igrid)%w(ixM^T, rho_)))
  end do

  call MPI_ALLREDUCE(tmp, wmax, 1, MPI_DOUBLE_PRECISION, &
       MPI_MAX, icomm, ierrmpi)

  print *, mype, global_time, "maximum density is", wmax
end subroutine my_analysis
```
