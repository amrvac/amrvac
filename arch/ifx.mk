arch := ifx

compile = mpiifx
f90_flags += -fpp -free -implicitnone

ifdef OPENMP
$(info Enabling OpenMP)
enabled += OPENMP
f90_flags += -qopenmp
endif

ifdef DEBUG
$(info Enable debugging symbols)
enabled += DEBUG
f90_flags += -g -O0
else
f90_flags += -O3
endif

