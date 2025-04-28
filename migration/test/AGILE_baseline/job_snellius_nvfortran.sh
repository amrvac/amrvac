#!/bin/bash
#SBATCH -p genoa
#SBATCH -N 16
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --job-name=NVHPC
#SBATCH --mail-type=ALL
#SBATCH --time=00-00:30:00
#SBATCH -o out # STDOUT 
#SBATCH -e err # STDERR

module load 2022
module load OpenMPI/4.1.4-NVHPC-22.7-CUDA-11.7.0

srun ./amrvac

