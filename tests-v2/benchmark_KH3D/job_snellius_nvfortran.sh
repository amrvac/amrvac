#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --gpus=1
#SBATCH --partition=gpu_a100
#SBATCH --mem=0
#SBATCH --job-name=NVHPC
#SBATCH --mail-type=ALL
#SBATCH --time=00-00:30:00
#SBATCH -o single_out-%j # STDOUT 
##SBATCH -e single_err-%j # STDERR

source $AMRVAC_DIR/.venv/bin/activate

module purge
module load 2023
module load OpenMPI/4.1.5-NVHPC-24.5-CUDA-12.1.1

export CUDA_LAUNCH_BLOCKING=1
export CUDA_MEMCHECK=1  # or use cuda-memcheck tool

##make
make clean-all
make -j16 arch=nvidia OPENACC=1

#export PGI_ACC_NOTIFY=16

srun ./amrvac

#cuda-memcheck ./amrvac

