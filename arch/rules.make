VACPP = $(AMRVAC_DIR)/src/vacpp.pl

%.o: %.f
	$(F90) $(F90FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS))

%.f: %.t
	$(VACPP) $(PPFLAGS) -d=$(NDIM)$(NDIM) $< > $(@)
