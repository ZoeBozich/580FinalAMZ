#!/bin/bash

OPTS="-Wall -O3"

echo "This script builds the final project test program!"
echo "By: Zoe Bozich, Alec Clark, and Missael Flores "
g++  $OPTS -o potential_solver potential_solver.cpp

./potential_solver
echo "Opening the plotting file..."
python3 ./final_plotter.py