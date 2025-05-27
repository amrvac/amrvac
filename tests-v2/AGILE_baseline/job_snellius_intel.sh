#!/bin/bash
#SBATCH -p genoa
#SBATCH -N 16
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --job-name=NVHPC
#SBATCH --mail-type=ALL
#SBATCH --time=00-01:10:00
#SBATCH -o out # STDOUT 
#SBATCH -e err # STDERR

module load 2022
module load impi/2021.6.0-intel-compilers-2022.1.0

srun ./amrvac

