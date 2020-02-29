//  file: integ_methods.h
// 
//  Header file for integ_methods.cpp
//
//
//  Programmer:  Ben Hall  furnstahl.1@osu.edu
//
//  Revision History:
//    29-Feb-2020 --- original version
//
//  To do:
//
//************************************************************************

//  begin: function prototypes 

extern double simpsons_rule ( int num_pts, float x_min, float x_max,
                       float (*integrand) (float x) );    // Simpson's rule
extern double milnes_rule ( int num_pts, float x_min, float x_max,
                            float (*integrand) (float x) );    // Milne's rule


//  end: function prototypes 
