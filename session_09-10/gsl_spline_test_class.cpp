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
using namespace std;    
#include "GslSpline.h"  // Header file for the GSL Spline class

inline double sqr (double z) {return z*z;}  // inline function for z^2

inline double bw (double x) {return 63900./(sqr(x-78.)+55.*55./4);}

int
main (void)
{
    
    // Create file
    ofstream spline_out;
    spline_out.open("spline_cubic.txt");
    
    // Initialize file
    spline_out << "    x     y_exact   y_spline   y'_spline  y''_spline"  << endl;
    
    // Set x and y
    double x_values[9] = {0., 25., 50., 75., 100., 125., 150., 175., 200.};
    double y_values[9] = {9.34, 17.9, 41.5, 83.5, 51.5, 21.5, 10.8, 6.29, 4.09};

  int npts = 9;

  // Make the spline object
  string type_linear = "linear";
  string type_poly = "polynomial";
  string type_cubic = "cubic";
  Spline my_spline (x_values, y_values, npts, type_cubic);

//  double x;
//  cout << "Enter x: ";
//  cin >> x;    // test point
    
  // Loop through x
    for (int x = 0; x <= 200; x = x+5)
    {
        
        // Evaluate the spline and derivatives
        double y = my_spline.y (x);
        double y_deriv = my_spline.yp (x);
        double y_deriv2 = my_spline.ypp (x);
        
        spline_out << fixed << setprecision(6)
        << x << "  " << bw(x) << "  " <<  y << "  " <<  y_deriv << "  " <<  y_deriv2 << "  " << endl;
    }
    
    spline_out.close(); // Close file

  return (0);      // successful completion 
}
