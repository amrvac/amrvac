################# Definitions for source files #############################
INCLUDES = mod_physicaldata.t mod_fix_conserve.t mod_connectivity.t	\
mod_global_parameters.t mod_usr_methods.t mod_forest.t

FOBJECTS = initialize_vars.t advance.t setdt.t amrini.t		\
amr_solution_node.t errest.t amrgrid.t amr_fix_conserve.t tvdlf.t tvd.t cd.t	\
fd.t load_balance.t connectivity.t boundary_conditions.t comm_lib.t forest.t	\
amr_neighbors.t amr_coarsen_refine.t refine.t coarsen.t set_B0.t selectgrids.t	\
geometry.t mod_ghostcells_update.t mod_usr.t

############################################################################
#  Dependencies:
############################################################################
mod_global_parameters.o: mod_physicaldata.o mod_connectivity.o
amr_fix_conserve.o : mod_fix_conserve.o mod_ghostcells_update.o
forest.o initialize_vars.o: mod_global_parameters.o
initialize_vars.o : mod_connectivity.o mod_timing.o mod_input_output.o	\
	mod_physics.o mod_ghostcells_update.o
advance.o setdt.o: mod_thermal_conduction.o
tvdlf.o: mod_physics.o mod_ghostcells_update.o
tvd.o: mod_physics.o
cd.o: mod_physics.o
fd.o: mod_physics.o
comm_lib.o: mod_ghostcells_update.o
amr_coarsen_refine.o: mod_ghostcells_update.o
refine.o: mod_physics.o
coarsen.o: mod_physics.o
convert.o: mod_ghostcells_update.o
collapse.o: mod_physics.o
slice.o: mod_physics.o

amrgrid.o: mod_ghostcells_update.o

mod_input_output.o: mod_physics.o mod_usr_methods.o
mod_physics.o: mod_global_parameters.o mod_physics_roe.o mod_physics_hllc.o	\
	mod_physics_ppm.o
mod_ghostcells_update.o: mod_global_parameters.o mod_physics.o
mod_usr_methods.o: mod_global_parameters.o
############################################################################
