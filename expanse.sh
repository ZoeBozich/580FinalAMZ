#!/bin/bash
#SBATCH --job-name="potential_solver"
#SBATCH --output="helloworld.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --account=sds154
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:01:00

L_x=10
L_y=10
R_x=5
R_y=5
rho_0=1
grid_points=1000  # Replace with required N

#This job runs with 1 nodes, 1 core per node and print hello world and date.
#Environment
module purge
module load slurm
module load cpu

OPTS="-Wall -O3"

g++ $OPTS -fopenmp -o potential_solver_parallel potential_solver_timer_parallel.cpp functions.h
./potential_solver_parallel $L_x $L_y $R_x $R_y $rho_0 $grid_points

echo "Opening the plotting file..."
python3 ./final_plotter.py
echo "Done"
