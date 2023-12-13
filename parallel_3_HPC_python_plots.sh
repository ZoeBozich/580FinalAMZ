#!/bin/bash
#SBATCH --job-name="potential_solver"
#SBATCH --output="potential_solver.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --account=sds154
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:03:00

# Fixed values for dimensions, R_x, R_y, rho_0
L_x=10
L_y=10
R_x=5
R_y=5
rho_0=1

module  load slurm

module load cpu/0.15.4 gcc/10.2.0


# Compile the program
echo "Compiling "
g++ -Wall -O3 -fopenmp -o parallel_3 parallel_3_noninteractive_itterative.cpp functions.h
chmod +x parallel_3


threads=1
grids=10
echo "Computing"
./parallel_3 $L_x $L_y $R_x $R_y $rho_0 $grids $threads

threads=1
grids=100
echo "Computing"
./parallel_3 $L_x $L_y $R_x $R_y $rho_0 $grids $threads

echo "Done"
