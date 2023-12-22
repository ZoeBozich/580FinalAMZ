#!/bin/bash

echo "This script builds the final project parallelized program!"
echo "By: Zoe Bozich, Alec Clark, and Missael Flores "

OPTS="-Wall -O3 -fopenmp"

g++  $OPTS -o parallel_solver parallel_solver.cpp

./parallel_solver
echo "Opening the plotting file..."
python3 ../Plotting/potential_plotter.py
