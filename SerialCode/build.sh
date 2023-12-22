#!/bin/bash

echo "This script builds the final project serial program!"
echo "By: Zoe Bozich, Alec Clark, and Missael Flores "

OPTS="-Wall -O3 -fopenmp"

g++  $OPTS -o potential_solver potential_solver.cpp

./potential_solver
echo "Opening the plotting file..."
python3 ../Plotting/potential_plotter.py
