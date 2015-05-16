############################################################################
#
#                    AMRVAC Makefile
#
############################################################################
#
# For the first time read carefully everything and do selections.
# Later you will probably need to change the compiler flags only.
# You may delete compiler flags for other machines once everything works.
#
# Select compiling options by adding and removing the leading # characters.
#
############################################################################

# The Makefile should always use the sh shell.
SHELL = /bin/sh
amrvac:

ARCH          = huey.defs
SRC           = $(AMRVAC_DIR)/src

include $(AMRVAC_DIR)/arch/$(ARCH)

############## setup.pl: #######################
PHYSICS       = srmhd

VPATH         =  $(SRC)/$(PHYSICS):$(SRC)
INCLUDE_DIRS  += -I.
################################################

############### VACPreProcessor ################
#Collect data for preprocessor set by setup.pl:
ndim=2
ndir=3
phi=3
z=2
nf=1
ndust=0
eos=default
USER=nul
# SETVAC READS UP TO THIS POINT

PPFLAGS  = -d=$(ndim)$(ndir) -phi=$(phi)
PPFLAGS += -z=$(z) -p=$(PHYSICS) -nf=$(nf)
PPFLAGS += -ndust=$(ndust)
PPFLAGS += -eos=$(eos)

# ... For the compilation rules:
#
VACPP =$(AMRVAC_DIR)/src/vacpp.pl #

############################################################################
# EXTENSION NAMES FOR SOURCE AND OBJECT FILES
#
# $F is the suffix for the generic AMRVAC precompiled Fortran files
# $O is the suffix for the object files, usually .o
F=.f
O=.o
I=.i

##############################################################################
######                                                                 #######
######   YOU PROBABLY DO NOT NEED TO CHANGE ANYTHING BELOW THIS LINE   #######
######   BUT IT MAY BE WORTHWHILE TO READ IT                           #######
######                                                                 #######
##############################################################################

################# Definitions for source files #############################

INCLUDES += amrvacpar$F amrvacusrpar$F mod_indices$F mod_physicaldata$F \
mod_fix_conserve$F mod_connectivity$F \
amrvacdef$F amrvacsettings$F $(MAKE_)

FOBJECTS =	\
mod_physicaldata$F	mod_indices$F		initialize_vars$F	\
amrvac$F		advance$F		setdt$F			            \
amrini$F				                        \
amr_solution_node$F	                                        \
errest$F		amrgrid$F		amr_fix_conserve$F	        \
amr_ghostcells$F	amrvacphys$F		amrvacusr$F		    \
tvdlf$F			tvd$F			cd$F			            \
fd$F                                                                     \
load_balance$F		connectivity$F		mod_fix_conserve$F	\
mod_connectivity$F	boundary_conditions$F   comm_lib$F		\
forest$F		mod_forest$F		amr_neighbors$F		    \
amr_coarsen_refine$F	geometry$F		refine$F		    \
coarsen$F		set_B0$F                    \
selectgrids$F                       \

OBJECTS =	\
mod_physicaldata.o	mod_indices.o		initialize_vars.o	\
amrvac.o		advance.o		setdt.o			            \
amrini.o                                                               \
amr_solution_node.o                                                    \
errest.o		amrgrid.o		amr_fix_conserve.o	        \
amr_ghostcells.o	amrvacphys.o		                    \
amrvacusr.o		                                            \
tvdlf.o			tvd.o			cd.o			            \
fd.o                                                                     \
load_balance.o		connectivity.o		mod_fix_conserve.o	\
mod_connectivity.o      boundary_conditions.o	comm_lib.o	\
forest.o		mod_forest.o		amr_neighbors.o		    \
amr_coarsen_refine.o	geometry.o		refine.o		    \
coarsen.o		set_B0.o                    \
selectgrids.o                       \

################### Include the IO-related objects in amrvacio dir #########
VPATH += :$(SRC)/amrvacio
-include $(SRC)/amrvacio/makefile
################### Include additional objects in physics dir ##############
-include $(SRC)/$(PHYSICS)/makefile
################### Include additional modules in modules dir ##############
VPATH += :$(SRC)/modules
-include $(SRC)/modules/makefile
################### Include additional files in local dir     ##############
-include localmakefile
############################################################################

################### Translation and dependency rules #######################

PREFOR= #
PREPROC= $(VACPP_)

.SUFFIXES:
.SUFFIXES: .t $F $O

$(FOBJECTS)	: $(PREPROC)
$(OBJECTS)	: $(INCLUDES)

# General precompilation rule
.t$(F):
	$(VACPP) $(PPFLAGS) $< $(PREFOR) > $@

# General Fortran compilation rule
$(F)$(O): ; $(F90) $(FFLAGS) $(F90FLAGS) $(INCLUDE_DIRS) $<

amrvac: $(OBJECTS)
	$(LINK) $(OBJECTS) $(LIBS) -o amrvac
#########################################################################
#  Dependencies:
#########################################################################

mod_physicaldata.o : mod_indices.o
amr_fix_conserve.o : mod_fix_conserve.o
forest.o initialize_vars.o : mod_forest.o
amrvac.o initialize_vars.o : mod_connectivity.o mod_timing.o
#########################################################################
#  Make clean
#########################################################################

clean: ;	\
	/bin/rm -f *.o  *.f *.mod \
	amrvac
