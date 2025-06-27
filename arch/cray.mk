arch := cray

compile = ftn
# disabled ftn messages: 878=multiple includes of same module
f90_flags += -ea -ef -eZ -ffree -eI -M878

ifdef OPENMP
$(info Enabling OpenMP)
enabled += OPENMP
f90_flags += -homp
else
f90_flags += -hnoomp
endif

ifdef OPENACC
$(info Enabling OpenACC)
f90_flags += -hacc
enabled += OPENACC
endif

ifdef DEBUG
$(info Enable debugging symbols)
enabled += DEBUG
f90_flags += -g -O0
else
f90_flags += -O3
endif

