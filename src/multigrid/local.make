# Copy this file to a problem directory, it will automatically be included. Make
# sure MG_DIR is set to the right location.

MG_DIR := $(HOME)/git/octree-mg/lib_2d
LIBS += omg
INC_DIRS += $(MG_DIR)
LIB_DIRS += $(MG_DIR)

vpath %.t $(AMRVAC_DIR)/src/multigrid

# Dependencies
amrvac: mod_multigrid_coupling.o
mod_usr.o: mod_multigrid_coupling.mod
mod_multigrid_coupling.mod: $(MG_DIR)/libomg.a

# Clean multigrid relates files
.PHONY: clean_mg
clean: clean_mg
clean_mg:
	$(RM) mod_multigrid_coupling.o mod_multigrid_coupling.mod
