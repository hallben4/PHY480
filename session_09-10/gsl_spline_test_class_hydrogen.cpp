//  file: gsl_spline_test_class.cpp
// 
//  Test program for the gsl spline routines using the Spline class
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision history:
//      02/10/09 -- created from gsl_cubic_spline_test.cpp
//
//  Notes:
//   * uses the GSL interpolation functions (see online documentation) 
//
//*****************************************************************
// include files
#include <iostream>    // cout and cin
#include <iomanip>     // manipulators like setprecision
#include <fstream>
#include <cmath>
#include <string>     // C++ strings
#include <vector>
#include <list>
using namespace std;    
#include "GslSpline.h"  // Header file for the GSL Spline class

inline double sqr (double z) {return z*z;}  // inline function for z^2

// Inline functions for wavefunction and its derivative
inline double wavefunction (double x) {return 2*x*exp(-x);}

inline double wavefunction_prime (double x) {return 2*(1-x)*exp(-x);}

inline double wavefunction_integral (double x) {return -2*(1+x)*exp(-x);}

int
main (void)
{
    
    // Create file
    ofstream spline_out;
    spline_out.open("spline_cubic.txt");
    
    // Initialize file
    spline_out << " x        y_diff       y'_diff    " << endl;
    
    // Number of points and endpoints of r
    int n = 100000;
    int start = 1;
    int stop = 4;
    
    // Set x and y
    double x_values[n];
    double y_values[n];
    
    for (double i = 0; i < n; i++)
    {
        double x = start + i*(stop-start)/(n-1);
        x_values[int(i)] = x;
        y_values[int(i)] = wavefunction(x);
    }

  // Make the spline object
  string type_linear = "linear";
  string type_poly = "polynomial";
  string type_cubic = "cubic";
    
  Spline my_spline (x_values, y_values, n, type_cubic);
    
    int count = 0;
    double area = 0;
  // Loop through x
    for (double i = 0; i < 10*n; i++)
    {
        double x = start + i*(stop-start)/(10*n-1);
        
        // Evaluate the spline and derivatives
        double y = my_spline.y(x);
        double y_deriv = my_spline.yp(x);
        
        // Add to area
        area += y*(stop-start)/(10*n-1);
        
        spline_out << fixed << setprecision(6)
        << x << "  "
        << fabs(wavefunction(x)-y) << "  "
        << fabs(wavefunction_prime(x)-y_deriv)
        << endl;
        
        if (fabs(wavefunction_prime(x)-y_deriv) >= 0.000001)
        {
            count += 1;
        }

    }
    
//    area /= (stop-start);
    
    cout << "Integral estimation: " << area << endl;
    cout << "Exact integral: " << wavefunction_integral(stop) - wavefunction_integral(start) << endl;
    
    // Report the number of absolute errors greater than 10^(-6)
    cout << "There are " << count << " non-zero entry(s)" << endl;
    
    spline_out.close(); // Close file

  return (0);      // successful completion 
}
