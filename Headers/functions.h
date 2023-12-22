#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include "nr.h"
#define W 12 // Width per entry in the data file
#define KEY "#################################" // Defines file header
#define DEFAULT_OUT "output.dat" // Default file used to store data
#define MAX_NM 100 // The max number of m,n used

// Below is a collection of random routines from various programs throughout the course
// Feel free to review or change things if necessary - Alec

// Clears cin buffer after an error
// (I found this is the only way to do it consistently on my machine)
void cin_clear()
{
    std::cin.clear();
    while (std::cin.get() != '\n')
    {
        continue;
    }
}

// Requests parameter (paramIn), must be a positive integer divisible by div
int fetch_pos_int(std::string paramIn, int div)
{
    int paramOut = 0;

    while (true)
    {
        std::cout << " - " << paramIn << ": ";
        std::cin >> paramOut;
        if (std::cin.fail() || paramOut <= 0)
        {
            std::cout << "  [" << paramIn << " must be a positive integer.]" << std::endl;
            cin_clear();
        }
        else if (paramOut % div != 0)
        {
            std::cout << "  [" << paramIn << " must be divisible by " << div << ".]" << std::endl;
        }
        else
        {
            return paramOut;
        }
    }
}

// Requests parameter (paramIn), must be a positive double
double fetch_pos_double(std::string paramIn)
{
    double paramOut = 0;

    while (true)
    {
        std::cout << " - " << paramIn << ": ";
        std::cin >> paramOut;
        if (std::cin.fail() || paramOut <= 0)
        {
            std::cout << "  [" << paramIn << " must be a positive number.]" << std::endl;
            cin_clear();
        }
        else
        {
            return paramOut;
        }
    }
}

// Requests parameter (paramIn), must be a positive double less than upperBound
double fetch_bounded_double(std::string paramIn, double upperBound)
{
    double paramOut = 0;

    while (true)
    {
        std::cout << " - " << paramIn << ": ";
        std::cin >> paramOut;
        if (std::cin.fail() || paramOut <= 0.0)
        {
            std::cout << "  [" << paramIn << " must be a positive number.]" << std::endl;
            cin_clear();
        }
        else if (paramOut > upperBound)
        {
            std::cout << "  [" << paramIn << " must be less than " << upperBound << ".]" << std::endl;
            cin_clear();
        }
        else
        {
            return paramOut;
        }
    }
}


// Creates the header for the output data file, communicates relevant information to python script.
// If this is changed, python programming must be changed to reflect this!
std::string create_header(Vec_I_DP L, Vec_I_DP R, Vec_I_DP r_0, double rho_0, int N, int N_mn, int M)
{
    std::string out = "";
    out += std::to_string(L[0]) + ",";
    out += std::to_string(L[1]) + ",";
    out += std::to_string(R[0]) + ",";
    out += std::to_string(R[1]) + ",";
    out += std::to_string(r_0[0]) + ",";
    out += std::to_string(r_0[1]) + ",";
    out += std::to_string(rho_0) + ",";
    out += std::to_string(N) + ",";
    out += std::to_string(N_mn) + ",";
    out += std::to_string(M) + ",";
    out += std::to_string(W);
    
    return out;
}

// From project 2, this appears within the integration
double booles_array(double fx[], double h, int N)
{
    double accum = 0;

    for (int i = 4; i < N; i += 4)
    {
        accum += 7.0 * fx[i - 4] + 32.0 * fx[i - 3] + 12.0 * fx[i - 2] + 32.0 * fx[i - 1] + 7.0 * fx[i];
    }
    return (2 * h / 45.0) * accum;
}

// rho as defined in the project directions, x_0 and y_0 chosen as L_x/2 and L_y/2.
double rho(double x, double y, Vec_I_DP R, Vec_I_DP L, Vec_I_DP r_0, double rho_0)
{
    double output = exp(-((x-r_0[0])/R[0]*(x-r_0[0])/R[0]+(y-r_0[1])/R[1]*(y-r_0[1])/R[1]));
    return rho_0*R[0]*R[1]*output/M_PI;
}

// Accepts m and n as they appear in the Fourier integrals, all the parameters needed for rho, and number of grid points N.
double fourier_double_int(Vec_I_DP R, Vec_I_DP L, Vec_I_DP r_0, double rho_0, int N, Vec_I_DP args, Vec_I_DP h, Vec_I_DP x, Vec_I_DP y)
{
    double y_sq[N];   // Represent a bunch of very small squares in y...
    double x_rect[N]; // Stacks of squares form thin rectangles in x direction...
                      // Adding up all the rectangles gives the total area!

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            y_sq[j] = rho(x[i], y[j], R, L, r_0, rho_0)*cos(x[i]*args[0])*cos(y[j]*args[1]);
        }
        x_rect[i] = booles_array(y_sq, h[1], N);
    }
    return booles_array(x_rect, h[0], N);
}

//Zoe's trying some stuff here
double zboole(double fx[], double h, int M){
	double ends = 7*(fx[0]+fx[M-1]);
	double total = 0;
	for (int i=5; i<M; i+=4){
		total += 32 * fx[i - 4] + 12 * fx[i - 3] + 32 * fx[i - 2] + 14 * fx[i - 1];
	}
	return (2 * h / 45.0) * (total+ends);
}

double zfourier(Vec_I_DP R, Vec_I_DP L, Vec_I_DP r_0, double rho_0, int N, Vec_I_DP args, Vec_I_DP h, Vec_I_DP x, Vec_I_DP y)
{
    int M =40; //fixed number of subdivisions for the Boole's approximation of the double integral
    double xh = L[0]/(1.0*M);
    double yh = L[1]/(1.0*M);
    double y_sq[M];   // Represent a bunch of very small squares in y...
    double x_rect[M]; // Stacks of squares form thin rectangles in x direction...
                      // Adding up all the rectangles gives the total area!

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            y_sq[j] = rho(i*xh, j*yh, R, L, r_0, rho_0)*cos(i*xh*args[0])*cos(j*yh*args[1]);
        }
        x_rect[i] = zboole(y_sq, yh, M);
    }
    return zboole(x_rect, xh, M);
}


// Silently solves for the potential with these parameters, returns the elapsed time to solve in seconds.
double timed_parallel_solver(Vec_I_DP L, Vec_I_DP R, Vec_I_DP r_0, double rho_0, int N, int N_mn, int M, int num_threads)
{
	double start_time = omp_get_wtime();
	omp_set_num_threads(num_threads);

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
    
	#pragma omp parallel for
    for (int m = 0; m < N_mn; m++)
    {
        for (int n = 0; n < N_mn; n++)
        {
            if (m == 0 && n == 0 ) continue;
            Vec_DP args(2);

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
    // No printing of the potential required

    double end_time = omp_get_wtime();
    return end_time - start_time;
}