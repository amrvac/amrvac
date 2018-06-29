# Create a static library
LIBPOIS3DFFT := libpois3dfft.a

OBJECTS_ABINIT := defs_basis.o defs_datatypes.o defs_xc.o drivexc.o xctetr.o	\
        xcwign.o xcxalp.o xchelu.o xcpzca.o xcspol.o xchcth.o xclb.o xcpbe.o	\
        invcb.o size_dvxc.o

OBJECTS_POIS3D := $(OBJECTS_ABINIT) poisson_solver.o timing.o dcopy.o

OBJECTS += $(OBJECTS_POIS3D)

$(LIBPOIS3DFFT): $(OBJECTS_POIS3D)
	$(RM) $@
	$(AR) rcs $@ $^

# Avoid circular dependency
$(filter-out defs_basis.o, $(OBJECTS_ABINIT)): \
	defs_basis.o

defs_datatypes.o: defs_basis.o

defs_xc.o: defs_datatypes.o

drivexc.o: defs_basis.o defs_xc.o

poisson_solver.o: psolver_main.f90 build_kernel.f90 psolver_base.f90 \
	xcenergy.f90 3dgradient.f90 fft3d.f90 scaling_function.f90
