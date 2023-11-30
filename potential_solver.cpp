#include "functions.h"


int main()
{
    // Input validation
    std::cout << "Welcome to the Poisson Equation solver!" << std::endl;
    std::cout << "This program solves the potential for 2D gaussian charge density." << std::endl;
    
    Vec_DP L(2);
    std::cout << "Please enter the dimensions of the box, in meters:" << std::endl;
    L[0] = fetch_pos_double("L_x (m)");// fetch functions validate the inputs based on certain parameters
    L[1] = fetch_pos_double("L_y (m)");

    Vec_DP R(2);
    std::cout << "Please enter the R_x and R_y parameters, in meters:" << std::endl;
    R[0] = fetch_bounded_double("R_x (m)", 0.5*L[0]);
    R[1] = fetch_bounded_double("R_y (m)", 0.5*L[1]);

    // TODO - Check/fix units
    std::cout << "Please enter the rho_0 parameter in Coulombs per meter^5:" << std::endl;
    double rho_0 = fetch_pos_double("rho_0 (C*m^5)");

    // TODO - Require N to be divisible by 4 for Boole's array approach
    std::cout << "Please enter the number of grid points N for the Fourier integrals approach:" << std::endl;
    int N = fetch_pos_int("N");

    // Creates a header to give relevant information to my python program
    // This is probably not a great approach to communicate with it, not sure how to improve!
    std::ofstream fp(DEFAULT_OUT);
    fp << KEY << std::endl;
    fp << create_header(L, R, rho_0, N) << std::endl;
    fp << KEY << std::endl;

    std::cout << "Calculating the potential..." << std::endl;

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
    return 0;
}
