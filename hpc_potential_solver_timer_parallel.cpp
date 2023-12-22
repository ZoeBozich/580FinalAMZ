#include "./Headers/functions.h"
#include <cstdlib> // For std::atof and std::atoi
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc < 8) {
        std::cerr << "Usage: " << argv[0] << " L_x L_y R_x R_y rho_0 N Threads" << std::endl;
        return 1;
    }

    Vec_DP L(2);
    L[0] = std::atof(argv[1]); // L_x
    L[1] = std::atof(argv[2]); // L_y

    Vec_DP R(2);
    R[0] = std::atof(argv[3]); // R_x
    R[1] = std::atof(argv[4]); // R_y

    double rho_0 = std::atof(argv[5]); // rho_0
    int N = std::atoi(argv[6]);        // N
    int num_threads = std::atoi(argv[7]); // Threads

    // Ensure num_threads does not exceed the maximum available
    num_threads = std::min(num_threads, omp_get_max_threads());
    omp_set_num_threads(num_threads);
    //Input validation
    std::cout << "Welcome to the Poisson Equation solver!" << std::endl;
    std::cout << "This program solves the potential for 2D gaussian charge density." << std::endl;

    //Creates a header to give relevant information to my python program
    //This is probably not a great approach to communicate with it, not sure how to improve!
    std::ofstream fp(DEFAULT_OUT);
    fp << KEY << std::endl;
    fp << create_header(L, R, r_0, rho_0, N, N_mn, M) << std::endl;
    fp << KEY << std::endl;

    std::cout << "Calculating the potential..." << std::endl;

    //Obtains the potential
    Mat_DP V(N,N);
    double c_mn, rho_mn, x, y;
	
	
	//For Missael's g14 laptop which has 8 cores and 16 threads
	//omp_set_num_threads(16); 
	double start_time = omp_get_wtime();

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            V[i][j] = 0.0;
        }
    }

    //Time for a messy block of code... This definitely needs fixing.
    double h_x = L[0]/(N-1);
    double h_y = L[1]/(N-1);
    
    int temp = 1;//helping us skip the m&&n ==0 case to avoid dividing by zero.
	
	#pragma omp parallel for
    for (int m = 0; m < MAX_NM; m++)
    {
        for (int n = temp; n < MAX_NM; n++)//int n=temp
        {
            // A lot of redundancy here compared to fourier_double_int routine.
            // Maybe can be merged into one big loop?
            double arg_x = m * M_PI / L[0];
            double arg_y = n * M_PI / L[1];
            rho_mn = fourier_double_int(R, L, rho_0, N, arg_x, arg_y, h_x, h_y);// change function to accept arg_x, arg_y, hx, hy, don't need m, n
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
        temp =0;
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
	double end_time = omp_get_wtime();
    std::cout << "Parallel section time taken: " << (end_time - start_time) << " seconds" << std::endl;

    return 0;
}
