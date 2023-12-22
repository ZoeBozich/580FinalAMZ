#!/bin/bash

echo "This script builds the plots of computation time versus N!"
echo "By: Zoe Bozich, Alec Clark, and Missael Flores "

OPTS="-Wall -O3 -fopenmp"

g++  $OPTS -o thread_tester thread_tester.cpp

./thread_tester
echo "Opening the plotting file..."
python3 ../Plotting/thread_plotter.py
