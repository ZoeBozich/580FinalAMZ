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

    Vec_DP r_0(2);
    std::cout << "Please enter the x_0 and x_0 parameters, in meters:" << std::endl;
    r_0[0] = fetch_bounded_double("x_0 (m)", L[0]);
    r_0[1] = fetch_bounded_double("y_0 (m)", L[1]);

    // TODO - Check/fix units
    std::cout << "Please enter the rho_0 parameter in Coulombs per meter^5:" << std::endl;
    double rho_0 = fetch_pos_double("rho_0 (C*m^5)");

    // TODO - Require N to be divisible by 4 for Boole's array approach
    std::cout << "Please enter the number of grid points N for the Fourier integrals approach:" << std::endl;
    int N = fetch_pos_int("N") + 1;

    // Creates a header to give relevant information to my python program
    // This is probably not a great approach to communicate with it, not sure how to improve!
    std::ofstream fp(DEFAULT_OUT);
    fp << KEY << std::endl;
    fp << create_header(L, R, r_0, rho_0, N) << std::endl;
    fp << KEY << std::endl;

    std::cout << "Calculating the potential..." << std::endl;

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
    
    int temp = 1;//helping us skip the m&&n ==0 case to avoid dividing by zero.
    for (int m = 0; m < N; m++)
    {
        for (int n = temp; n < N; n++)//int n=temp
        {
            // A lot of redundancy here compared to fourier_double_int routine.
            // Maybe can be merged into one big loop?
            args[0] = m * M_PI / L[0];
            args[1] = n * M_PI / L[1];
            rho_mn = zfourier(R, L, r_0, rho_0, N, args, h, x, y);
            c_mn = rho_mn / (args[0]*args[0] + args[1]*args[1]);
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    // Note the 4/L_xL_y, to avoid using 2/sqrt(...) twice.
                    V[i][j] += 4.0*c_mn*cos(args[0]*x[i])*cos(args[1]*y[j])/(L[0]*L[1]);
                }
            }
        }
        temp = 0;
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
    return 0;
}
