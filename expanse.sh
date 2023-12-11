#!/bin/bash

OPTS="-Wall -O3"

echo "This script builds the final project test program!"
echo "By: Zoe Bozich, Alec Clark, and Missael Flores "
g++ $OPTS -fopenmp -o potential_solver_parallel potential_solver_edits_parallel.cpp functions.h
./potential_solver
echo "Opening the plotting file..."
python3 ./final_plotter.py
