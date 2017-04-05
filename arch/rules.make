VACPP := $(AMRVAC_DIR)/src/vacpp.pl

# Disable built-in make rules
.SUFFIXES:

# How to compile modules
mod_%.o mod_%.mod: mod_%.f
	$(F90) $(F90FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS))

# How to generate object files
%.o: %.f
	$(F90) $(F90FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS))

# How to translate .t source files to normal Fortran
%.f: %.t
	$(VACPP) $(PPFLAGS) -d=$(NDIM) $< > $(@)

# How to generate executables
%: %.o
	$(LINK) $(F90FLAGS) $^ -o $@ $(addprefix -L,$(LIB_DIRS)) $(addprefix -l,$(LIBS))


