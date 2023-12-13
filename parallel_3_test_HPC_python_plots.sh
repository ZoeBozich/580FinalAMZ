
L_x=10
L_y=10
R_x=5
R_y=5
rho_0=1


# Compile the program
echo "Compiling "
g++ -Wall -O3 -fopenmp -o parallel_3 parallel_3_noninteractive_itterative.cpp functions.h
chmod +x parallel_3

threads=2
grids=10
echo "Computing"
./parallel_3 $L_x $L_y $R_x $R_y $rho_0 $grids $threads

threads=2
grids=100
echo "Computing"
./parallel_3 $L_x $L_y $R_x $R_y $rho_0 $grids $threads

echo "Done"