//  file: ode_test_class.cpp
//
//  C++ Program to test the ode differential equation solver from
//   the gsl numerical library using the Ode class wrapper for gsl.
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision history:
//     02/10/09  Original version based on ode_test.cpp
//
//  Notes:  
//   * Example taken from the GNU Scientific Library Reference Manual
//      Edition 1.1, for GSL Version 1.1 9 January 2002
//      URL: gsl/ref/gsl-ref_23.html#SEC364
//   * Compile and link with GslOde class files:
//       g++ -Wall -c ode_test_class.cpp 
//       g++ -Wall -c GslOde.cpp
//       g++ -Wall -o ode_test_class ode_test_class.o GslOde.o -lgsl -lgslcblas
//
//********************************************************************

// The following details are taken from the GSL documentation
// 
// The following program solves the second-order nonlinear 
//  Van der Pol oscillator equation (see background notes),
//
//     x"(t) + \mu x'(t) (x(t)^2 - 1) + x(t) = 0
//
// This can be converted into a first order system suitable for 
//  use with the library by introducing a separate variable for 
//  the velocity, v = x'(t).  We assign x --> y[0] and v --> y[1].
//  So the equations are:
// x' = v                  ==>  dy[0]/dt = f[0] = y[1]
// v' = -x + \mu v (1-x^2) ==>  dy[1]/dt = f[1] = -y[0] + mu*y[1]*(1-y[0]*y[0])
//
//*********************************************************************

// include files
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>    // C++ stringstream class (can omit iostream)
#include <cmath>
using namespace std;

#include "GslOde.h"   // Ode class for gsl (include Ode and Rhs)
   
// Class for the right side of the Van der Pol equation.
//  Derived from the Rhs class.
//  The VdP equation depends on one parameter mu.
class Rhs_VdP : public Rhs 
{
  public:
    Rhs_VdP (double mu_passed) {mu = mu_passed; num_eqs = 2;};
    ~Rhs_VdP () {};
    virtual int rhs (double t, const double y[], double f[]);
    virtual int jacobian (double t, const double y[], double *dfdy, 
                          double dfdt[]);
    double get_mu () {return mu;};                        
  private:
    double mu;    // Van der Pol parameter
};

// Function to evolve the differential equation and print the output
int evolve_and_print(Ode &vdp_ode,  Rhs_VdP &vdp_rhs,
                     const double x0, const double v0,
                     const double x1, const double v1,
                     const double x2, const double v2,
                     const double tmin, const double tmax, const double delta_t);


//*************************** main program ****************************
int
main ()
{
    double x0 = 0.0;          // initial conditions
    double v0 = 0.0;
    double x1 = 0.0;          // initial conditions
    double v1 = 0.0;
    double x2 = 0.0;          // initial conditions
    double v2 = 0.0;

  const double eps_abs = 1.e-8;    // absolute error requested 
  const double eps_rel = 1.e-10;   // relative error requested 

  double tmin = 0.;        // starting t value 
  double tmax = 100.;      // final t value 
  double delta_t = 0.01;   // step size in time
  
  double mu = 2;            // parameter for the diff-eq
  double nu = 3;            // parameter for the diff-eq
  Rhs_VdP vdp_rhs_1 (mu);   // set up the right side of the diff-eq
  Rhs_VdP vdp_rhs_2 (nu);   // set up a second instance of the right side of the diff-eq
    
  Ode vdp_ode_1 (vdp_rhs_2, eps_abs, eps_rel, "rk45");
  x0 = -1.5;
  v0 = 2.0;
  x1 = 0.1;
  v1 = 0.0;
  x2 = 1.0;
  v2 = 0.0;
  evolve_and_print(vdp_ode_1, vdp_rhs_2, x0, v0, x1, v1, x2, v2, tmin, tmax, delta_t);
  
  return 0;
}

//*******************************************************
int
evolve_and_print(Ode &vdp_ode, Rhs_VdP &vdp_rhs,
                 const double x0, const double v0,
                 const double x1, const double v1,
                 const double x2, const double v2,
                 const double tmin, const double tmax, const double delta_t)
{
  double y0[2];     // current solution vector
  y0[0] = x0;       // initial x value
  y0[1] = v0;       // initial v value
  double y1[2];     // current solution vector
  y1[0] = x1;       // initial x value
  y1[1] = v1;       // initial v value
  double y2[2];     // current solution vector
  y2[0] = x2;       // initial x value
  y2[1] = v2;       // initial v value

  double t0 = tmin;         // initialize t
  double t1 = tmin;         // initialize t
  double t2 = tmin;         // initialize t
    
  // Set up a file name with the initial values
  ostringstream my_stringstream0;  // declare a stringstream object
  my_stringstream0 << "ode_test_class" << "_mu_" << setprecision(2)
                  << vdp_rhs.get_mu()
                  << "_x0_" << setprecision(2) << y0[0]
                  << "_v0_" << setprecision(2) << y0[1] << ".dat";
    
    ostringstream my_stringstream1;  // declare a stringstream object
    my_stringstream1 << "ode_test_class" << "_mu_" << setprecision(2)
                    << vdp_rhs.get_mu()
                    << "_x0_" << setprecision(2) << y1[0]
                    << "_v0_" << setprecision(2) << y1[1] << ".dat";
    
    ostringstream my_stringstream2;  // declare a stringstream object
    my_stringstream2 << "ode_test_class" << "_mu_" << setprecision(2)
                    << vdp_rhs.get_mu()
                    << "_x0_" << setprecision(2) << y2[0]
                    << "_v0_" << setprecision(2) << y2[1] << ".dat";
    
    
    
  ofstream my_out0;    // now open a stream to a file for output
  my_out0.open(my_stringstream0.str().c_str());
    ofstream my_out1;    // now open a stream to a file for output
    my_out1.open(my_stringstream1.str().c_str());
    ofstream my_out2;    // now open a stream to a file for output
    my_out2.open(my_stringstream2.str().c_str());

  // print initial values and column headings
  my_out0 << "# Running ode_test with x0 = " << setprecision(2) << y0[0]
         << " and v0 = " << setprecision(2) << y0[1] << endl;
  my_out0 << "#     t            x           v   " << endl;
  my_out0 << scientific << setprecision (5) << setw (12) << t0 << " "
         << setw (12) << y0[0] << " " << setw (12) << y0[1] << endl;
    
    // print initial values and column headings
    my_out1 << "# Running ode_test with x0 = " << setprecision(2) << y1[0]
           << " and v0 = " << setprecision(2) << y1[1] << endl;
    my_out1 << "#     t            x           v   " << endl;
    my_out1 << scientific << setprecision (5) << setw (12) << t1 << " "
           << setw (12) << y1[0] << " " << setw (12) << y1[1] << endl;
    
    // print initial values and column headings
    my_out2 << "# Running ode_test with x0 = " << setprecision(2) << y2[0]
           << " and v0 = " << setprecision(2) << y2[1] << endl;
    my_out2 << "#     t            x           v   " << endl;
    my_out2 << scientific << setprecision (5) << setw (12) << t2 << " "
           << setw (12) << y2[0] << " " << setw (12) << y2[1] << endl;

  // step to tmax from tmin 
  double h = 1e-6;    // starting step size for ode solver 
  for (double t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
  {

    while (t0 < t_next)  // evolve from t to t_next
    {
        vdp_ode.evolve ( &t0, t_next, &h, y0 );
    }
    while (t1 < t_next)  // evolve from t to t_next
    {
        vdp_ode.evolve ( &t1, t_next, &h, y1 );
    }
    while (t2 < t_next)  // evolve from t to t_next
    {
        vdp_ode.evolve ( &t2, t_next, &h, y2 );
    }

    // print at t = t_next
    my_out0 << scientific << setprecision (5) << setw (12) << t0 << " "
            << setw (12) << y0[0] << " " << setw (12) << y0[1] << endl;
    my_out1 << scientific << setprecision (5) << setw (12) << t1 << " "
            << setw (12) << y1[0] << " " << setw (12) << y1[1] << endl;
    my_out2 << scientific << setprecision (5) << setw (12) << t2 << " "
            << setw (12) << y2[0] << " " << setw (12) << y2[1] << endl;
      
  }
    
  my_out0.close();
  my_out1.close();
  my_out2.close();

  return (0);    // successful completion 
}


//*******************************************************
//*******************************************************
// 
// Define the array of right-hand-side functions y[i] to be integrated.
//  The equations are:
// x' = v                  ==>  dy[0]/dt = f[0] = y[1]
// v' = -x + \mu v (1-x^2) ==>  dy[1]/dt = f[1] = -y[0] + mu*y[1]*(1-y[0]*y[0])
int
Rhs_VdP::rhs (double, const double y[], double f[])
{
  // std::cout << "rhs called with y[0] = " << y[0] << endl;
  // evaluate the right-hand-side functions at t 
  f[0] = y[1];
  f[1] = -y[0] + mu * y[1] * (1. - y[0] * y[0]);

  return 0;    // successful completion
}

//
// Define the Jacobian matrix of df_i/dy_j for i,j = {0,1}
int
Rhs_VdP::jacobian (double, const double y[], double *, double dfdt[])
{
  // fill the Jacobian matrix
  set_jacobian (0, 0, 0.0);                            // df[0]/dy[0] = 0 
  set_jacobian (0, 1, 1.0);                            // df[0]/dy[1] = 1 
  set_jacobian (1, 0, -2.0 * mu * y[0] * y[1] - 1.0);  // df[1]/dy[0] 
  set_jacobian (1, 1, -mu * (y[0] * y[0] - 1.0));      // df[1]/dy[1] 

  // set explicit t dependence of f[i] (none here)
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;

  return 0;    // successful completion
}
