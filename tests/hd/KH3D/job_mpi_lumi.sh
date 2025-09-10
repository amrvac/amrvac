#!/bin/bash
#SBATCH -p small-g
#SBATCH -N 1
#SBATCH --account <project number here>
#SBATCH --gpus-per-node 8
#SBATCH --ntasks-per-node 8
#SBATCH --cpus-per-task 1
#SBATCH -t 00:30:00
#SBATCH -o job_mpi_lumi.out

module load LUMI/24.03 partition/G PrgEnv-cray rocm/6.0.3

# See https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/lumig-job/

# Enable GPU-aware MPI support
export MPICH_GPU_SUPPORT_ENABLED=1

# Ensure proper GPU to CPU binding
cat << EOF > select_gpu
#!/bin/bash

export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
exec \$*
EOF

chmod +x ./select_gpu

# GPU 0 is bound to CPU 49, GPU 1 to CPU 57, etc.
CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"

srun --cpu-bind=${CPU_BIND} ./select_gpu ./amrvac
rm -f ./select_gpu
