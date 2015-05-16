#!/bin/bash
#SBATCH --ntasks=4
#SBATCH --constraint=dual
#SBATCH --cpus-per-task=1
#SBATCH --job-name=amrvac3D
#SBATCH --mail-type=ALL
#SBATCH --partition=parallel
#SBATCH --time=00-02:00:00
# 
module load openmpi/intel-14.0.3/1.8.1
#
mpiexec -n 4 ./amrvac -i amrvac3D.par
