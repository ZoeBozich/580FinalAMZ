# 580FinalAMZ
phys 580 final project 


############### parallel 1 ###############################

#Notes: Standard inputs via terminal & OPENMP Parallelization	

#Contents: parallel_1_basic_potential_solver.cpp	


g++ -Wall -O3 -fopenmp -o parallel_1 parallel_1_basic_potential_solver.cpp functions.h		

############### Parallel 2 ################################# 

#Notes: NEW shell script that loops for every core & gridd on the NEW cpp file based MODIFIED main parallelized cpp file 
Now has Unique output filename generation

#Contents: 

parallel_2_noninteractive_itterative.sh

parallel_2_noninteractive_itterative_potential_solver.cpp

############### Parallel 3 ################################# IN PROGRESS YALLS


#Notes: NEW shell script for HPC implementation of "Parallel 2" 

Its tricky to use python on hpc so this script will generate unique dat files for post-visualization
