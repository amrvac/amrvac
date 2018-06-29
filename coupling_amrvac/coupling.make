SRC_F90 := m_data_structures.f90 m_build_tree.f90 m_load_balance.f90	\
m_ghost_cells.f90 m_allocate_storage.f90 m_mrgrnk.f90 m_restrict.f90	\
m_communication.f90 m_prolong.f90 m_multigrid.f90 m_octree_mg.f90	\
m_laplacian.f90 m_vlaplacian.f90 m_helmholtz.f90 m_vhelmholtz.f90	\
m_diffusion.f90 m_free_space.f90

ifneq ($(NDIM), 1)
OBJECTS += $(SRC_F90:%.f90=%.o) mod_multigrid_coupling.o
else
OBJECTS += mod_multigrid_coupling.o
endif

MY_FLAGS := -cpp -DNDIM=$(NDIM)

# Get .f90 files from the octree-mg/src folder
vpath %.f90 $(AMRVAC_DIR)/external_libs/octree-mg/src

# Get .t coupling files
vpath %.t $(AMRVAC_DIR)/external_libs/octree-mg/coupling_amrvac

# How to get .o object files from .f90 source files
%.o: %.f90
	$(F90) -c -o $@ $< $(F90FLAGS) $(MY_FLAGS) $(addprefix -I,$(INC_DIRS))

# How to get .mod files from .f90 source files (remake only if they have been
# removed, otherwise assume they are up to date)
m_%.mod: m_%.f90 m_%.o
	@test -f $@ || $(F90) -c -o $(@:.mod=.o) $< $(F90FLAGS) \
	$(MY_FLAGS) $(addprefix -I,$(INC_DIRS))

# AMRVAC dependencies
amrvac.o: mod_multigrid_coupling.mod
amr_coarsen_refine.o: mod_multigrid_coupling.mod

# Coupling dependency
ifneq ($(NDIM), 1)
mod_multigrid_coupling.o: m_octree_mg.mod
endif

m_free_space.o: INC_DIRS += $(HOME)/git/poisson_3d_fft
amrvac: LIBS += pois3dfft
amrvac: LIB_DIRS += $(HOME)/opt/sw/poisson_3d_fft

# Other dependencies
m_allocate_storage.o: m_data_structures.mod
m_allocate_storage.o: m_ghost_cells.mod
m_allocate_storage.o: m_prolong.mod
m_allocate_storage.o: m_restrict.mod
m_build_tree.o: m_data_structures.mod
m_communication.o: m_data_structures.mod
m_communication.o: m_mrgrnk.mod
m_diffusion.o: m_data_structures.mod
m_diffusion.o: m_helmholtz.mod
m_diffusion.o: m_multigrid.mod
m_diffusion.o: m_vhelmholtz.mod
m_free_space.o: m_data_structures.mod
m_free_space.o: m_multigrid.mod
m_free_space.o: m_restrict.mod
m_ghost_cells.o: m_communication.mod
m_ghost_cells.o: m_data_structures.mod
m_helmholtz.o: m_data_structures.mod
m_laplacian.o: m_data_structures.mod
m_load_balance.o: m_data_structures.mod
m_multigrid.o: m_data_structures.mod
m_multigrid.o: m_ghost_cells.mod
m_multigrid.o: m_helmholtz.mod
m_multigrid.o: m_laplacian.mod
m_multigrid.o: m_prolong.mod
m_multigrid.o: m_restrict.mod
m_multigrid.o: m_vhelmholtz.mod
m_multigrid.o: m_vlaplacian.mod
m_octree_mg.o: m_allocate_storage.mod
m_octree_mg.o: m_build_tree.mod
m_octree_mg.o: m_communication.mod
m_octree_mg.o: m_data_structures.mod
m_octree_mg.o: m_ghost_cells.mod
m_octree_mg.o: m_helmholtz.mod
m_octree_mg.o: m_load_balance.mod
m_octree_mg.o: m_multigrid.mod
m_octree_mg.o: m_prolong.mod
m_octree_mg.o: m_restrict.mod
m_octree_mg.o: m_vhelmholtz.mod
m_prolong.o: m_communication.mod
m_prolong.o: m_data_structures.mod
m_restrict.o: m_communication.mod
m_restrict.o: m_data_structures.mod
m_vhelmholtz.o: m_data_structures.mod
m_vlaplacian.o: m_data_structures.mod
