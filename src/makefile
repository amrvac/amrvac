default:
	@echo '=== Call make or setup.pl from your work-directory, not src!   ==='
	@echo '=== Please read ../doc/getting_started.md                      ==='
################# Definitions for source files #############################

INCLUDES += amrvacpar$F amrvacusrpar$F mod_indices$F mod_physicaldata$F \
mod_fix_conserve$F mod_connectivity$F					\
mod_global_parameters$F amrvacsettings$F $(MAKE_)

FOBJECTS =	                                                        \
mod_physicaldata$F	mod_indices$F		initialize_vars$F	\
amrvac$F		advance$F		setdt$F	               	\
amrini$F				                        	\
amr_solution_node$F	                                        	\
errest$F		amrgrid$F		amr_fix_conserve$F      \
amr_ghostcells$F	amrvacphys$F		amrvacusr$F		\
tvdlf$F			tvd$F			cd$F			\
fd$F                                                                    \
load_balance$F		connectivity$F		mod_fix_conserve$F	\
mod_connectivity$F	boundary_conditions$F   comm_lib$F		\
forest$F		mod_forest$F		amr_neighbors$F		\
amr_coarsen_refine$F	refine$F                                        \
coarsen$F		set_B0$F             				\
selectgrids$F           geometry$F		mod_global_parameters$F	\

OBJECTS =								\
mod_physicaldata.o	mod_indices.o		initialize_vars.o	\
amrvac.o		advance.o		setdt.o			\
amrini.o                                                                \
amr_solution_node.o                                                     \
errest.o		amrgrid.o		amr_fix_conserve.o	\
amr_ghostcells.o	amrvacphys.o		                        \
amrvacusr.o		                                            	\
tvdlf.o			tvd.o			cd.o			\
fd.o                                                                    \
load_balance.o		connectivity.o		mod_fix_conserve.o	\
mod_connectivity.o      boundary_conditions.o	comm_lib.o	        \
forest.o		mod_forest.o		amr_neighbors.o		\
amr_coarsen_refine.o	refine.o                                        \
coarsen.o		set_B0.o                                        \
selectgrids.o           geometry.o		mod_global_parameters.o	\


############################################################################
#  Dependencies:
############################################################################
mod_global_parameters.o: mod_indices.o mod_physicaldata.o mod_connectivity.o
mod_physicaldata.o : mod_indices.o
amr_fix_conserve.o : mod_fix_conserve.o
forest.o initialize_vars.o : mod_forest.o mod_global_parameters.o
amrvac.o initialize_vars.o : mod_connectivity.o mod_timing.o
############################################################################
