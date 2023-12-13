#!/bin/bash

# Fixed values for dimensions, R_x, R_y, rho_0
Lx=10
Ly=10
Rx=5
Ry=5
rho_0=1

# Array of N values and thread counts
grids=(10 100 1000 10000 100000 1000000 10000000)
threads=(1 2 3 8 16 32 64 128)

# Compile the program
g++ -Wall -O3 -fopenmp -o parallel_2 parallel_2_noninteractive_itterative_potential_solver.cpp functions.h

# Run tests
for n in "${grids[@]}"; do
    for t in "${threads[@]}"; do
        echo "Running for N=$n with $t threads..."
        ./parallel_2 $Lx $Ly $Rx $Ry $rho_0 $n $t >> results.txt
    done
done
