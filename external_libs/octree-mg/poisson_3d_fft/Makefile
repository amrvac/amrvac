F90 ?= mpif90
F90FLAGS ?= -O2 -cpp

PROGRAMS := PSolver

all: $(PROGRAMS)

include definitions.make

clean:
	$(RM) $(OBJECTS_POIS3D) $(PROGRAMS) $(LIBPOIS3DFFT) *.mod

# Rules for compiling object files
%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

# Rule for compiling programs
%: %.f90
	$(F90) $(F90FLAGS) $< $(LIBPOIS3DFFT) -o $@

# Dependencies
$(PROGRAMS): $(LIBPOIS3DFFT)

