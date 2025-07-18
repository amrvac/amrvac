#!/bin/bash
# 
# Use this wrapper for MPI tasks that see more than one GPU per task
# The wrappers sets one GPU to visible, equal to the local rank index
# (or counting down from the highest instance on the interactive node)

if [[ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]]; then
    rank=$OMPI_COMM_WORLD_LOCAL_RANK
elif [[ -n "$SLURM_LOCALID" ]]; then
    rank=$SLURM_LOCALID
else
    echo "Unknown MPI rank, did you run this script with mpirun/srun?"
    exit 1
fi


# Use GPU MIG on interactive GPU node
if [[ $(hostname) == int3.local.snellius.surf.nl ]]; then
    MIG=($(nvidia-smi -L | sed -nr "s|^.*UUID:\s*(MIG-[^)]+)\)|\1|p"))
    
    # 7 MIG instances per GPU, count down from the last one = 27
    let GPU=27-$rank
    export CUDA_VISIBLE_DEVICES=${MIG[$GPU]}
else
    export CUDA_VISIBLE_DEVICES=${rank}
fi

$@
