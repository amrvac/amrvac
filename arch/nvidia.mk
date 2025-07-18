arch := nvidia

compile = mpif90
f90_flags += -cpp -Mfree

ifdef NVTX
$(info Enabling NVTX)
enabled += NVTX
f90_flags += -DNVTX
link_flags += -lnvToolsExt
endif

ifdef OPENMP
$(info Enabling OpenMP)
enabled += OPENMP
f90_flags += -fopenmp
endif

ifdef OPENACC
$(info Enabling OpenACC)
f90_flags += -Wall -acc=gpu -Mvect=levels:5 -Minline
enabled += OPENACC
ifdef NOGPUDIRECT
$(info Disabling direct GPU-GPU copies)
f90_flags += -DNOGPUDIRECT
enabled += NOGPUDIRECT
endif
endif

ifdef INFO
f90_flags += -Minfo=all
endif

ifdef DEBUG
$(info Enable debugging symbols)
enabled += DEBUG
f90_flags += -g -O0
else
f90_flags += -O3 -fast
endif

link_flags += $(f90_flags)
