#!/bin/bash
#SBATCH -p small
#SBATCH -N 1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=224G
#SBATCH --job-name=agile-baseline
#SBATCH --time=00:10:00
#SBATCH --output=%x-%j.out
#SBATCH --account=XXXXXXXXXX

module load LUMI/23.09 partition/C PrgEnv-cray

echo "Nodes: $SLURM_JOB_NUM_NODES"
srun ./amrvac
