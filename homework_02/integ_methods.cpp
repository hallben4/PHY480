//  file: integ_methods.cpp
//
//  Routines for Simpson, Milne, and GSL routine
//                                                                     
//  Programmer:  Ben Hall
//
//  Revision history:
//      29-Feb-2020  original version
//
//  Notes:
//   * define with floats to emphasize round-off error  
//   * compile with:  "g++ -Wall -c integ_routines.cpp" or makefile
//   * adapted from: "Projects in Computational Physics" by Landau and Paez
//             copyrighted by John Wiley and Sons, New York               
//             code copyrighted by RH Landau  
//   * equation for interval h = (b-a)/(N-1) with x_min=a and x_max=b
// 
//************************************************************************

// include files
#include <cmath>
#include "integ_methods.h"   // integration routine prototypes

//************************************************************************


// Integration using Simpson's rule
double simpsons_rule ( int num_pts, float x_min, float x_max,
                      float (*integrand) (float x) )
{  
   double interval = ((x_max - x_min)/float(num_pts - 1));  // called h in notes
   double sum=  0.;  // initialize integration sum to zero
   
   for (int n=2; n<num_pts; n+=2)                // loop for odd points  
   {
     double x = x_min + interval * double(n-1);
     sum += (4./3.) * interval * integrand(x);
   }
   for (int n=3; n<num_pts; n+=2)                // loop for even points  
   {
     double x = x_min + interval * double(n-1);
     sum += (2./3.)*interval * integrand(x);
   }   
   // add in the endpoint contributions   
   sum +=  (interval/3.) * (integrand(x_min) + integrand(x_max));	
   
   return (sum);
}  

//************************************************************************

// Integration using Milne's rule
double milnes_rule ( int num_pts, float x_min, float x_max,
                     float (*integrand) (float x) )
{
    double interval = ((x_max - x_min)/double(num_pts - 1));  // called h in notes
    double sum=  0.;  // initialize integration sum to zero
    
    for (int n=2; n<num_pts; n+=2)                // loop for even points
    {
        double x = x_min + interval * float(n-1);
        sum += (64./45.) * interval * integrand(x);
    }
    
    for (int n=3; n<num_pts; n+=4)                // loop for 4*i-1
    {
        double x = x_min + interval * double(n-1);
        sum += (24./45.) * interval * integrand(x);
    }
    
    for (int n=5; n<num_pts; n+=4)                // loop for 4*i+1
    {
        double x = x_min + interval * double(n-1);
        sum += (28./45.) * interval * integrand(x);
    }

    // add in the endpoint contributions
    sum +=  (14./45.) * interval * (integrand(x_min) + integrand(x_max));
    
    return (sum);
}

//************************************************************************


