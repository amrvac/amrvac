ARCH          = huey.defs
SRC           = $(AMRVAC_DIR)/src

include $(AMRVAC_DIR)/arch/$(ARCH)

VPATH         =  $(SRC)/modules
INCLUDE_DIRS  += -I.

driver:

OBJ = driver.o  mod_odeint.o

############################################################################
%.f: %.t
	cp $< $@ 

%.o: %.f $(DEPS)
	$(F90) $(F90FLAGS) $(FFLAGS) $(INCLUDE_DIRS) -o $@ $<

driver.o: mod_odeint.o

driver: $(OBJ)
	$(LINK) $(OBJ) -o odeint 


#########################################################################
#  Make clean
#########################################################################

clean: ;        \
        /bin/rm -f *.o  *.f *.mod \
        odeint
