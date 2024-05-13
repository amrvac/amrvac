.phony: all clean

all: laplace_openacc_ghostcell

clean:
	$(RM) laplace_openacc_ghostcell *.o

# Determine compiler brand
compiler_version = $(shell $(FC) --version)
compiler_brand = $(word 1, $(compiler_version))

ifeq ($(compiler_brand), GNU)
	FFLAGS ?= -O2 -fopenacc -Wall -g
	ifeq ($(DEBUG), 1)
		FFLAGS += -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan
	endif
else ifeq ($(compiler_brand), ifort)
	FFLAGS ?= -warn all -O2 -stand f08 -assume realloc-lhs
else ifeq ($(compiler_brand), nvfortran)
	FFLAGS ?= -Wall -O2 -acc -fast
endif

# How to get .o object files from .f90 source files
%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get .mod files from .f90 source files (remake only if they have been
# removed, otherwise assume they are up to date)
%.mod: %.f90 %.o
	@test -f $@ || $(FC) -c -o $(@:.mod=.o) $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get executables from .o object files
%: %.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
