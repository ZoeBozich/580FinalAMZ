#!/bin/bash
#SBATCH --job-name="potential_solver"
#SBATCH --output="potential_solver.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --account=sds154
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH -t 03:30:00

# Fixed values for dimensions, R_x, R_y, rho_0
Lx=10
Ly=10
Rx=5
Ry=5
rho_0=1

# Array of N values and thread counts
grids=(10 100 1000 10000 100000 1000000 10000000)
threads=(128 64 32 16 8 4 3 2 1)


# Compile the program
g++ -Wall -O3 -fopenmp -o parallel_3 parallel_2_noninteractive_itterative_potential_solver.cpp functions.h

# Run tests
for n in "${threads[@]}"; do
    for t in "${grids[@]}"; do
        echo "Running for $t threads with N=$n"
        ./parallel_3 $Lx $Ly $Rx $Ry $rho_0 $n $t >> results.txt
    done
done
