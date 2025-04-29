compile = mpif90
f90_flags += -ffree-form -fimplicit-none -Wall -cpp
f90_flags += -Wno-unused-dummy-argument	\
	     -Wno-unused-function -Wno-unused -Wno-uninitialized \
	     -Wno-zerotrip -Wno-target-lifetime

ifdef OPENMP
$(info Enabling OpenMP)
enabled += OPENMP
f90_flags += -fopenmp
endif

ifdef DEBUG
$(info Enable debugging symbols)
enabled += DEBUG
f90_flags += -g -O0
else
f90_flags += -O3
endif

