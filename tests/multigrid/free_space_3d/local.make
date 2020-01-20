# Instructions for using the free-space solver in MPI-AMRVAC
#
# 1. Clone the octree-mg repository somewhere. Below, it is assumed to be in
# ~/git/octree-mg
#
# 2. Type `make` in `the octree-mg/poisson_3d_fft` folder
#
# 3. Copy the file `octree-mg/single_module/m_free_space_3d.f90` to the folder
# with your MPI-AMRVAC user code, and rename it to `m_free_space_3d.t`
#
# 3. Create a file `local.make` with the following content

# Adjust this to where you have installed octree-mg
OCTREE_MG_DIR:=$${HOME}/git/octree-mg

INC_DIRS+=$(OCTREE_MG_DIR)/poisson_3d_fft
LIB_DIRS+=$(OCTREE_MG_DIR)/poisson_3d_fft
LIBS+=pois3dfft
mod_usr.o: m_free_space_3d.o
amrvac: m_free_space_3d.o
