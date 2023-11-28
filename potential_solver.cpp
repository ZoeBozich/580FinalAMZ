#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
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

// Requests parameter (paramIn), must be a positive integer
int fetch_pos_int(std::string paramIn)
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

/* Probably not needed, but in case we wanted flexibility for an output file of choice

// Requests an output file name, if it already exists, requests approval to overwrite
std::ofstream fetch_file()
{
    std::string fileInput;
    std::cout << " - Enter file name: ";
    std::cin >> fileInput;

    std::ifstream in_fp(fileInput);
    if (in_fp.good())
    {
        char overwrite;
        std::cout << fileInput <<" already exists. Overwrite? (y/n): ";
        std::cin >> overwrite;

        if (overwrite == 'y' || overwrite == 'Y')
        {
            std::cout << "Overwriting " << fileInput << "..." << std::endl;
        }
        else
        {
            std::cout << "Restarting..." << std::endl;
            return fetch_file();
        }
    }
    in_fp.close();

    std::ofstream out_fp(fileInput.c_str());
    return out_fp;
}
*/

// Creates the header for the output data file, communicates relevant information to python script.
// If this is changed, python programming must be changed to reflect this!
std::string create_header(Vec_I_DP L, Vec_I_DP R, double rho_0, int N)
{
    std::string out = "";
    out += std::to_string(L[0]) + ",";
    out += std::to_string(L[1]) + ",";
    out += std::to_string(R[0]) + ",";
    out += std::to_string(R[1]) + ",";
    out += std::to_string(rho_0) + ",";
    out += std::to_string(N) + ",";
    out += std::to_string(W);
    
    return out;
}

// From project 2, this appears within the integration
double booles_array(double fx[], double h, int N)
{
    double accum = 0;

    for (int i = 4; i <= N; i += 4)
    {
        accum += 7.0 * fx[i - 4] + 32.0 * fx[i - 3] + 12.0 * fx[i - 2] + 32.0 * fx[i - 1] + 7.0 * fx[i];
    }
    return (2 * h / 45.0) * accum;
}

// rho as defined in the project directions, x_0 and y_0 chosen as L_x/2.
double rho(double x, double y, Vec_I_DP R, Vec_I_DP L, double rho_0)
{
    double output = exp(-((x-0.5*L[0])/R[0]*(x-0.5*L[0])/R[0]+(y-0.5*L[1])/R[1]*(y-0.5*L[1])/R[1]));
    return rho_0*R[0]*R[1]*output/M_PI;
}

// Accepts m and n as they appear in the Fourier integrals, all the parameters needed for rho,
// and number of grid points N. Right now this is a bit inefficient, needs fixes.
double fourier_double_int(int m, int n, Vec_I_DP R, Vec_I_DP L, double rho_0, int N)
{
    // This is lifted directly from my project 2, hence uses C++ arrays rather than NR ones
    // Should something be changed for consistency?
    double y_sq[N];   // Represent a bunch of very small squares in y...
    double x_rect[N]; // Stacks of blocks form thin rectangles in x direction...
                      // Adding up all the rectanges gives the total area!

    double h_x = L[0]/(N-1);
    double h_y = L[1]/(N-1);
    double arg_x = m * M_PI / L[0];
    double arg_y = n * M_PI / L[1];

    for (int i = 0; i <= N; i++)
    {
        double x_0 = i * h_x;
        for (int j = 0; j <= N; j++)
        {
            double y_0 = j * h_y;
            y_sq[j] = rho(x_0, y_0, R, L, rho_0)*cos(x_0*arg_x)*cos(y_0*arg_y);
        }
        x_rect[i] = booles_array(y_sq, h_y, N);
    }
    return booles_array(x_rect, h_x, N);
}

int main()
{
    // Input validation - AC
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

    // TODO - Check/fix units
    std::cout << "Please enter the rho_0 parameter in Coulombs per meter^5:" << std::endl;
    double rho_0 = fetch_pos_double("rho_0 (C*m^5)");

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

    for (int m = 1; m < MAX_NM; m++)
    {
        for (int n = 1; n < MAX_NM; n++)
        {
            // A lot of redundancy here compared to fourier_double_int routine.
            // Maybe can be merged into one big loop?
            double arg_x = m * M_PI / L[0];
            double arg_y = n * M_PI / L[1];
            rho_mn = fourier_double_int(m, n, R, L, rho_0, N);
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
    // NOTE: Despite what the project directions say, I think the plot turns out much better
    // if you replace the cosines above and in the fourier routine with sines... But still not perfect

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
            fp << std::setw(W) << V[i][j];
            fp << std::endl;
        }
    }
    std::cout << "The data file has been generated!" << std::endl;
    return 0;
}
