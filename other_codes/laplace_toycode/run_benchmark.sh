#!/bin/bash

n_gc=1

for ((grid_nx=256, n_it=16384; grid_nx<=8192; grid_nx*=2, n_it/=4)); do
    for ((block_nx=4; block_nx<=grid_nx; block_nx*=2)); do
        ./laplace $grid_nx $grid_nx $block_nx $block_nx $n_gc $n_it
    done
done
