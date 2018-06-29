NDIM := 3
LIBOMG := libomg.a

.PHONY: all clean
all: $(LIBOMG)

vpath %.f90 ../src
vpath %.f90 ../poisson_3d_fft

include ../makerules.make
include ../src/definitions.make
include ../poisson_3d_fft/definitions.make

$(LIBOMG): $(OBJECTS)
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) $(OBJECTS) $(LIBOMG) *.mod
