#include "functions.h"
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>

int main(int argc, char *argv[]) {
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0] << " Lx Ly Rx Ry rho_0 N num_threads\n";
        return 1;
    }
	
    // Input validation
    std::cout << "Welcome to the Poisson Equation solver!" << std::endl;
    std::cout << "This program solves the potential for 2D gaussian charge density." << std::endl;
   
	Vec_DP L(2), R(2);
    L[0] = std::stod(argv[1]); // L_x
    L[1] = std::stod(argv[2]); // L_y
    R[0] = std::stod(argv[3]); // R_x
    R[1] = std::stod(argv[4]); // R_y
    double rho_0 = std::stod(argv[5]);
    int N = std::stoi(argv[6]);
    int num_threads = std::stoi(argv[7]);
	
	//Printing N, number of grids
	std::cout << "Number of grids: " << N << std::endl;

   
	//Set number of threads entered by the User
	omp_set_num_threads(num_threads);
	//Print maximum number of threads
	std::cout << "Maximum number of threads available: " << omp_get_max_threads() << std::endl;
	//Print the number of threads being used
	std::cout << "Number of threads used: " << num_threads << std::endl;

	// Generate a unique filename based on N and num_threads
    std::ostringstream filename;
    filename << "output_N" << N << "_threads" << num_threads << ".dat";
    std::string filename_str = filename.str();

    // Open the file with the new filename
    std::ofstream fp(filename_str);
    if (!fp) {
        std::cerr << "Error: Could not open file " << filename_str << std::endl;
        return 1;
    }


    // Creates a header to give relevant information to my python program
    // This is probably not a great approach to communicate with it, not sure how to improve!
    std::ofstream fp(DEFAULT_OUT);
    fp << KEY << std::endl;
    fp << create_header(L, R, rho_0, N) << std::endl;
    fp << KEY << std::endl;

    std::cout << "Calculating the potential..." << std::endl;

	//TIMER BEGINS
	double start_time = omp_get_wtime();

    // Obtains the potential
    Mat_DP V(N,N);
    double c_mn, rho_mn, x, y;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            V[i][j] = 0.0;
        }
    }

    // Time for a messy block of code... This definitely needs fixing.
    double h_x = L[0]/(N-1);
    double h_y = L[1]/(N-1);
    
    int temp = 1;//helping us skip the m&&n ==0 case to avoid dividing by zero.
    
	//PARALLELIZATION HERE 
	#pragma omp parallel for collapse(2)
	for (int m = 0; m < MAX_NM; m++){
        for (int n = temp; n < MAX_NM; n++){//int n=temp
			if (m == 0 && n == 0) continue; // Skip the m && n == 0 case
				// A lot of redundancy here compared to fourier_double_int routine.
				// Maybe can be merged into one big loop?
				double arg_x = m * M_PI / L[0];
				double arg_y = n * M_PI / L[1];
				rho_mn = zfourier(R, L, rho_0, N, arg_x, arg_y, h_x, h_y);// change function to accept arg_x, arg_y, hx, hy, don't need m, n
				c_mn = rho_mn / (arg_x*arg_x + arg_y*arg_y);
				for (int i = 0; i < N; i++)
				{
					x = i * h_x;
					for (int j = 0; j < N; j++)
					{
						y = j * h_y;
						// Note the 4/L_xL_y, to avoid using 2/sqrt(...) twice.
						V[i][j] += 4.0*c_mn*cos(arg_x*x)*cos(arg_y*y)/(L[0]*L[1]);
                }
            } 
        }
    }

    std::cout << "Done! Generating data file..." << std::endl;

    // Prints the results to file DEFAULT_OUT
    for (int i = 0; i < N; i++)
    {
        x = i * h_x;
        for (int j = 0; j < N; j++)
        {
            y = j * h_y;
            fp << std::setw(W) << x;
            fp << std::setw(W) << y;
            fp << std::setw(W) << rho(x,y,R,L,rho_0);
            fp << std::setw(W) << V[i][j];
            fp << std::endl;
        }
    }
    std::cout << "The data file has been generated!" << std::endl;
	
	//TIMER END
    double end_time = omp_get_wtime();
	//TIME CALCULATION
    double duration = end_time - start_time;
    //OUTPUT TIME
    std::cout << "Calculation completed in " << duration << " seconds." << std::endl;
	//Adding empty le=ine for readability
    std::cout << '\n';

	return 0;
}
