#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --gpus=1
#SBATCH --partition=gpu_h100
#SBATCH --mem=0
#SBATCH --job-name=NVHPC
#SBATCH --mail-type=ALL
#SBATCH --time=00-00:30:00
#SBATCH -o out # STDOUT 
#SBATCH -e err # STDERR

module purge
module load 2023
module load OpenMPI/4.1.5-NVHPC-24.5-CUDA-12.1.1


$AMRVAC_DIR/setup.pl -d=3 -arch=nvidia_acc
#make allclean
make

#export PGI_ACC_NOTIFY=16

srun ./amrvac

