# This makefile is included in the ../lib_2d and ../lib_3d folder makefiles
ifndef NDIM
$(error NDIM is not set, type make in top or the lib folders)
endif

OBJECTS += m_data_structures.o m_build_tree.o m_load_balance.o m_ghost_cells.o	\
m_allocate_storage.o m_mrgrnk.o m_restrict.o m_communication.o m_prolong.o	\
m_multigrid.o m_octree_mg.o m_laplacian.o m_vlaplacian.o m_helmholtz.o		\
m_vhelmholtz.o m_diffusion.o m_free_space.o

ifeq ($(NDIM), 3)
m_free_space.o: poisson_solver.mod
endif

# Dependencies
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
m_free_space.o: m_ghost_cells.mod
m_free_space.o: m_multigrid.mod
m_free_space.o: m_prolong.mod
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
m_octree_mg.o: m_free_space.mod
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
