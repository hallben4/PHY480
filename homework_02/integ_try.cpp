//  file: integ_try.cpp
//
//  This is a test program for basic integration methods.               
//                                                                     
//  Programmer: Ben Hall
//
//  Revision history:
//      29-Feb-2020  original version, for 780.20 Computational Physics
//
//  Notes:
//   * define with floats to emphasize round-off error  
//   * compile with:  "g++ -Wall -c integ_test.cpp"
//   * adapted from: "Projects in Computational Physics" by Landau and Paez  
//             copyrighted by John Wiley and Sons, New York               
//             code copyrighted by RH Landau                           
//

// include files
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#include "integ_methods.h"	// prototypes for integration routines

float my_integrand (float x);

//const double ME = 2.7182818284590452354E0;    // Euler's number

int
main ()
{
  // set up the integration specifiction
  const int max_intervals = 501;	// maximum number of intervals
  const float lower = 0.0;	// lower limit of integration
  const float upper = 1.0;	// upper limit of integration

//  const double answer = 1. - 1. / ME;    // the "exact" answer for the test
  double result = 0.;  // approximate answer

  // open output file streams
  ofstream integ_hw2_out ("integ_hw2.dat");	// save data in integ_hw2.dat
    
  integ_hw2_out << "#  N     Simpsons     Milne        GSL          " << endl;
  integ_hw2_out << "#----------------------------------------------" << endl;
    

  // Milne's rule requires 4*i+1 points
  for (int i = 5; i <= max_intervals; i += 4)
  {
    integ_hw2_out << setw(4) << log10(i);

    result = simpsons_rule (i, lower, upper, &my_integrand);
//    integ_hw2_out << "  " << scientific << log10(fabs((result - answer)/answer));
    integ_hw2_out << "  " << scientific << setprecision(15) << log10(fabs(result));

    result = milnes_rule (i, lower, upper, &my_integrand);
//    integ_hw2_out << "  " << scientific << log10(fabs((result - answer)/answer));
    integ_hw2_out << "  " << scientific << setprecision(15) << log10(fabs(result));

    integ_hw2_out << endl;
  }

  cout << "data stored in integ_hw2.dat\n";
  integ_hw2_out.close ();

  return (0);
}

// the function we want to integrate 
float
my_integrand (float x)
{
  return sin(x);
}
