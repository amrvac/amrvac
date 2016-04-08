#!/bin/bash
#SBATCH --ntasks=192
#SBATCH --mem 64000
#SBATCH --constraint=dual
#SBATCH --cpus-per-task=1
#SBATCH --job-name=bench8N
#SBATCH --mail-type=ALL
#SBATCH --partition=test
#SBATCH --time=00-01:00:00
#SBATCH -o output/out8N # STDOUT
#SBATCH -e output/err8N # STDERR
# 
module load openmpi/intel-14.0.3/1.8.1
#
mpiexec -n 192 ./amrvac
