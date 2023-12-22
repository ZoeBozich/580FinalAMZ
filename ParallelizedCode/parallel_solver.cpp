#include "../Headers/functions.h"

int main()
{
    // Input validation - fetch functions validate the inputs based on certain parameters
    std::cout << "Welcome to the Poisson Equation solver!" << std::endl;
    std::cout << "This program solves the potential for 2D gaussian charge density." << std::endl;
    
    Vec_DP L(2);
    std::cout << "Please enter the dimensions of the box, in meters:" << std::endl;
    L[0] = fetch_pos_double("L_x (m)");
    L[1] = fetch_pos_double("L_y (m)");

    Vec_DP R(2);
    std::cout << "Please enter the R_x and R_y parameters, in meters:" << std::endl;
    R[0] = fetch_bounded_double("R_x (m)", 0.5*L[0]);
    R[1] = fetch_bounded_double("R_y (m)", 0.5*L[1]);

    Vec_DP r_0(2);
    std::cout << "Please enter the x_0 and y_0 parameters, in meters:" << std::endl;
    r_0[0] = fetch_bounded_double("x_0 (m)", L[0]);
    r_0[1] = fetch_bounded_double("y_0 (m)", L[1]);

    std::cout << "Please enter the rho_0 parameter in Coulombs per meter^5:" << std::endl;
    double rho_0 = fetch_pos_double("rho_0 (C*m^5)");

    std::cout << "Please enter the number of grid points N for the potential matrix:" << std::endl;
    int N = fetch_pos_int("N",1);

    std::cout << "Please enter the maximum m and n values for the cosine terms, N_mn:" << std::endl;
    int N_mn = fetch_pos_int("N_mn",1);

    std::cout << "Please enter the number of subdivisions M for Boole's Method:" << std::endl;
    int M = fetch_pos_int("M",4) + 1;

	std::cout << "Please enter the number of threads to use (max " << omp_get_max_threads() << "): " << std::endl;
	int num_threads = fetch_pos_int("Threads", 1);
	if (num_threads > omp_get_max_threads())
    {
        std::cout << "   [Exceeds the max number of threads " << omp_get_max_threads() << ", will use this value.]"  << std::endl;
		num_threads = omp_get_max_threads();
	}
    
	omp_set_num_threads(num_threads);

    // V calculation, saved as a NR matrix
    std::cout << "Calculating the potential..." << std::endl;
    
	double start_time = omp_get_wtime();

    Mat_DP V(N,N);
    Vec_DP x(N), y(N), h(2);
    Vec_DP x_b(M), y_b(M), h_b(2); // For Boole's method
    
    h[0] = L[0]/(N-1);
    h[1] = L[1]/(N-1);

    h_b[0] = L[0]/(M-1);
    h_b[1] = L[1]/(M-1);

    for (int i = 0; i < N; i++)
    {
        x[i] = i * h[0];
        y[i] = i * h[1];
    }

    for (int i = 0; i < M; i++)
    {
        x_b[i] = i * h_b[0];
        y_b[i] = i * h_b[1];
    }
    
    
    // int temp = 1;// helping us skip the m&&n==0 case to avoid dividing by zero.
    // I had to remove this to produce good plots, since otherwise multiple threads would overwrite it
	#pragma omp parallel for
    for (int m = 0; m < N_mn; m++)
    {
        for (int n = 0; n < N_mn; n++)
        {
            if (m == 0 && n == 0 ) continue;
            Vec_DP args(2); // Moved these inside the loop so that threads do not conflict

            args[0] = m * M_PI / L[0];
            args[1] = n * M_PI / L[1];
            double rho_mn = fourier_double_int(R, L, r_0, rho_0, M, args, h, x_b, y_b);
            double c_mn = rho_mn / (args[0]*args[0] + args[1]*args[1]);
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    // Note the 4/L_xL_y, to avoid using 2/sqrt(...) twice.
                    V[i][j] += 4.0*c_mn*cos(args[0]*x[i])*cos(args[1]*y[j])/(L[0]*L[1]);
                }
            }
        }
    }

    std::cout << "Done! Generating data file..." << std::endl;

    // Creates a header to give relevant information to python program
    std::ofstream fp(DEFAULT_OUT);
    fp << KEY << std::endl;
    fp << create_header(L, R, r_0, rho_0, N, N_mn, M) << std::endl;
    fp << KEY << std::endl;
    fp << std::setw(W) << "x (m)";
    fp << std::setw(W) << "y (m)";
    fp << std::setw(W) << "rho (C/m^3)";
    fp << std::setw(W) << "V (C/m)";
    fp << std::endl;

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

    double end_time = omp_get_wtime();
    std::cout << "The data file has been generated! Elapsed time: " << end_time - start_time << " seconds." << std::endl;
    return 0;
}