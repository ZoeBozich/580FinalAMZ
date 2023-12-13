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

module purge
module  load slurm
module  load cpu
module  load python/3.8.5

# Compile the program
echo "Compiling "
g++ -Wall -O3 -fopenmp -o parallel_2 parallel_2_noninteractive_itterative_potential_solver.cpp functions.h

threads=3
grids=1
echo "Computing"
./parallel_2 $L_x $L_y $R_x $R_y $rho_0 $grids $threads

threads=4
grids=10
echo "Computing"
./parallel_2 $L_x $L_y $R_x $R_y $rho_0 $grids $threads

echo "Opening the plotting file..."
python3 ./iterative_final_plotter.py
echo "Done"