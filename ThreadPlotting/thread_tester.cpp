#include "../Headers/functions.h"

int main()
{
    // This is hardcoded in unfortunately, feel free to change around
    int num_threads = 4; int num_Ns = 20;
    int threads[num_threads] = {1,2,4,8};
    int Ns[num_Ns] = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};

    // Input validation - fetch functions validate the inputs based on certain parameters
    std::cout << "Welcome to the thread count tester!" << std::endl;
    std::cout << "This program presents how the parallel program performs for different N and thread counts." << std::endl;
    
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

    std::cout << "Please enter the maximum m and n values for the cosine terms, N_mn:" << std::endl;
    int N_mn = fetch_pos_int("N_mn",1);

    std::cout << "Please enter the number of subdivisions M for Boole's Method:" << std::endl;
    int M = fetch_pos_int("M",4) + 1;

    // Creates a header to give relevant information to python program
    std::ofstream fp(DEFAULT_OUT);
    fp << KEY << std::endl;
    fp << create_header(L, R, r_0, rho_0, num_Ns, N_mn, M) << std::endl;
    fp << KEY << std::endl;
    fp << std::setw(W) << "num_threads";
    fp << std::setw(W) << "N";
    fp << std::setw(W) << "t (s)";
    fp << std::endl;
    
    std::cout << "Calculating the potential and recording computation time..." << std::endl;
    
    for (int i = 0; i < num_threads; i++)
    {
        for (int j = 0; j < num_Ns; j++)
        {
            fp << std::setw(W) << threads[i];
            fp << std::setw(W) << Ns[j];
            fp << std::setw(W) << timed_parallel_solver(L,R,r_0,rho_0,Ns[j],N_mn,M,threads[i]);
            fp << std::endl;
        }
    }

    std::cout << "The data file has been generated!" << std::endl;
    return 0;
}