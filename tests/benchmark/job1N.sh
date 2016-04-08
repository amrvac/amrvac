#!/bin/bash
#SBATCH --ntasks=24
#SBATCH --mem 64000
#SBATCH --constraint=dual
#SBATCH --cpus-per-task=1
#SBATCH --job-name=bench1N
#SBATCH --mail-type=ALL
#SBATCH --partition=test
#SBATCH --time=00-01:00:00
#SBATCH -o output/out1N # STDOUT
#SBATCH -e output/err1N # STDERR
# 
module load openmpi/intel-14.0.3/1.8.1
#
mpiexec -n 24 ./amrvac
