#include <omp.h>
#include <chrono>
#include <iostream>
#include "functions.h"

int main(int argc, char* argv[])
{
	// Check if the correct number of arguments is passed
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0] << " Lx Ly Rx Ry rho_0 N threads\n";
        return 1;
    }
	
    // Input validation
    // Parse command-line arguments
    double Lx = std::atof(argv[1]);
    double Ly = std::atof(argv[2]);
    double Rx = std::atof(argv[3]);
    double Ry = std::atof(argv[4]);
    double rho_0 = std::atof(argv[5]);
    int N = std::atoi(argv[6]) + 1;
    int threads = std::atoi(argv[7]);

    std::cout << "Number of grid points: " << N << std::endl;
	
    Vec_DP L(2), R(2), r_0(2);
    L[0] = Lx;
    L[1] = Ly;
    R[0] = Rx;
    R[1] = Ry;
	
    // Hardcoded x_0 and y_0 values
    double x_0 = 0.5 * L[0]; // Example: half of Lx
    double y_0 = 0.5 * L[1]; // Example: half of Ly
    r_0[0] = x_0;
    r_0[1] = y_0;
	
    //Set number of threads entered by the User
    omp_set_num_threads(threads);
	
    //Print maximum number of threads
    std::cout << "Maximum number of threads available: " << omp_get_max_threads() << std::endl;
    //Print the number of threads being used
    std::cout << "Number of threads used: " << threads << std::endl;
	
    //Unique output filename generation
    // Generate a unique FILENAME based on N and threads
    std::ostringstream filename;
    filename << "output_N" << N << "_threads" << threads << ".dat";
    std::string filename_str = filename.str();
    // Open the FILE with the new filename
    std::ofstream fp(filename_str);
    if (!fp) {
        std::cerr << "Error: Could not open file " << filename_str << std::endl;
        return 1;
    }
	
    fp << KEY << std::endl;
    fp << create_header(L, R, r_0, rho_0, N) << std::endl;
    fp << KEY << std::endl;

    std::cout << "Calculating the potential..." << std::endl;
    //TIMER BEGINS
    double start_time = omp_get_wtime();	
    
    // Obtains the potential
    double c_mn, rho_mn;
    Mat_DP V(N,N);
    Vec_DP x(N), y(N), h(2), args(2);
    
    h[0] = L[0]/(N-1);
    h[1] = L[1]/(N-1);
	
    for (int i = 0; i < N; i++)
    {
        x[i] = i * h[0];
        y[i] = i * h[1];
    }
    
    //int temp = 1;//helping us skip the m&&n ==0 case to avoid dividing by zero.
	
    //PARALLELIZATION
    #pragma omp parallel for collapse(2)
    for (int m = 0; m < N; m++){
		for (int n = 0; n < N; n++) {
			if (m == 0 && n == 0) continue; // Skip the case when both m and n are zero

            // A lot of redundancy here compared to fourier_double_int routine.
            // Maybe can be merged into one big loop?
            args[0] = m * M_PI / L[0];
            args[1] = n * M_PI / L[1];
            rho_mn = zfourier(R, L, r_0, rho_0, N, args, h, x, y);// change function to accept arg_x, arg_y, hx, hy, don't need m, n
            c_mn = rho_mn / (args[0]*args[0] + args[1]*args[1]);
            for (int i = 0; i < N; i++){
                for (int j = 0; j < N; j++){
                    // Note the 4/L_xL_y, to avoid using 2/sqrt(...) twice.
                    V[i][j] += 4.0*c_mn*cos(args[0]*x[i])*cos(args[1]*y[j])/(L[0]*L[1]);
                }
            }
        }
    }

    std::cout << "Done! Generating data file..." << std::endl;

    // Prints the results to file DEFAULT_OUT
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fp << std::setw(W) << x[i];
            fp << std::setw(W) << y[j];
            fp << std::setw(W) << rho(x[i],y[j],R,L,r_0,rho_0);
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
