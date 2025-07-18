#!/bin/bash
#SBATCH -p gpu_h100
#SBATCH -N 1
#SBATCH --gpus-per-node 4
#SBATCH --ntasks-per-node 4
#SBATCH --cpus-per-task 1
#SBATCH -t 00:10:00
#SBATCH -o mpi_local.out

module load 2023
module load NVHPC/24.5-CUDA-12.1.1 OpenMPI/4.1.5-NVHPC-24.5-CUDA-12.1.1 Nsight-Systems/2023.3.1 UCX-CUDA/1.14.1-GCCcore-12.3.0-CUDA-12.1.1

source $AMRVAC_DIR/.venv/bin/activate

#make clean
make -j8 arch=nvidia OPENACC=1

srun ./amrvac
