#!/bin/bash
#for local machine

L_x=10
L_y=10
R_x=5
R_y=5
rho_0=1

OPTS="-Wall -O2"
g++ $OPTS -fopenmp -o potential_solver_parallel potential_solver_edits_parallel.cpp functions.h

# Array of grid points to test
grid_points_array=(10 100 1000 10000 )

# Run tests for each grid point
for grid_points in "${grid_points_array[@]}"
do
    for num_threads in 2 4 6 # Add more thread counts if needed
    do
        echo "Running with grid_points=$grid_points and num_threads=$num_threads"
        ./potential_solver_parallel $L_x $L_y $R_x $R_y $rho_0 $grid_points $num_threads 
		
    done
done
