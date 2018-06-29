# A custom compiler can be set with make F90=...
F90 ?= mpif90

ifeq ($(F90), mpif90)
# If no flags have been set, use these as the default
F90FLAGS ?= -O2 -std=f2008 -fopenmp -Wall -g -Wno-unused-dummy-argument	\
-Wno-unused-function

# For preprocessing
PPFLAGS := -cpp -DNDIM=$(NDIM)

# Add debugging flags when DEBUG = 1
ifeq ($(DEBUG), 1)
F90FLAGS += -fcheck=all -ffpe-trap=invalid,zero,overflow \
-pedantic -finit-real=snan
endif
endif

ifeq ($(F90), mpiifort)
# If no flags have been set, use these as the default
F90FLAGS ?= -warn all -O2 -stand f08 -assume realloc-lhs

# For preprocessing
PPFLAGS := -fpp -DNDIM=$(NDIM)
endif

# How to get .o object files from .f90 source files
%.o: %.f90
	$(F90) -c -o $@ $< $(F90FLAGS) $(PPFLAGS) $(addprefix -I,$(INCDIRS))

# How to get .mod files from .f90 source files (remake only if they have been
# removed, otherwise assume they are up to date)
%.mod: %.f90 %.o
	@test -f $@ || $(F90) -c -o $(@:.mod=.o) $< $(PPFLAGS) $(F90FLAGS) $(addprefix -I,$(INCDIRS))

# How to get executables from .o object files
%: %.o
	$(F90) -o $@ $^ $(F90FLAGS) $(PPFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

