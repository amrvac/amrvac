################# Definitions for source files #############################
INCLUDES = mod_global_parameters.t mod_usr_methods.t mod_forest.t	\
mod_physicaldata.o mod_connectivity.t mod_constants.t mod_variables.t	\
mod_basic_types.t

FOBJECTS = mod_advance.t setdt.t amrini.t amr_solution_node.t errest.t		\
amrgrid.t mod_finite_volume.t mod_tvd.t mod_finite_difference.t load_balance.t	\
connectivity.t boundary_conditions.t comm_lib.t forest.t amr_neighbors.t	\
amr_coarsen_refine.t refine.t coarsen.t set_B0.t selectgrids.t mod_geometry.t	\
mod_ghostcells_update.t mod_small_values.t mod_ppm.t mod_mp5.t mod_limiter.t	\
mod_fix_conserve.t mod_source.t mod_initialize.t mod_amrvac.t mod_bc_data.t mod_lookup_table.t


%.o: module_deps.txt

module_deps.txt:
	list_module_deps.sh > module_deps.txt


############################################################################
#  Dependencies of INCLUDES
############################################################################

mod_global_parameters.o: mod_connectivity.mod
mod_global_parameters.o: mod_physicaldata.mod
mod_global_parameters.o: mod_constants.mod
mod_global_parameters.o: mod_variables.mod
mod_global_parameters.o: mod_basic_types.mod
mod_variables.o: mod_basic_types.mod
mod_usr_methods.o: mod_global_parameters.mod
