//  file: area.cpp
//
//  This program calculates the area of a circle, given the radius.
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision history:
//      02-Jan-2004  original version, for 780.20 Computational Physics
//      01-Jan-2010  updates to "To do" wishlist
//      12-Jan-2016  comment out "using namespace std;"
//
//  Notes:  
//   * compile with:  "g++ -o area.x area.cpp"
//
//  To do:
//   1. output the answer with higher precision (more digits)
//   2. use a "predefined" value of pi or generate it with atan
//   3. define an inline square function
//   4. split the calculation off into a function (subroutine)
//   5. output to a file (and/or input from a file)
//   6. add checks of the input (e.g., for non-positive radii)
//   7. rewrite using a Circle class
//
//*********************************************************************// 

// include files
#include <iostream>	     // this has the cout, cin definitions
#include <iomanip>       // in order to use setprecision
#include <math.h>        // in order to use atan
#include <fstream>       // in order to write to a file
using namespace std;     // if omitted, then need std::cout, std::cin

//*********************************************************************//

double pi = 4 * atan (1);   // define pi in terms of atan

double square_function (double x);   // function prototype

// Inline function

//inline
//double
//square_function (double x)
//{
//    return (x * x);
//}

int
main ()
{
  double radius;    // every variable is declared as int or double or ...
    
  double precision = 10;     // Number of digits of area to display

  cout << "Enter the radius of a circle: ";	// ask for radius
  cin >> radius;

  double area = pi * square_function(radius);// standard area formula
    
  // Write radius to a file
    
  ofstream radius_file;
  radius_file.open("radius.txt");
  radius_file << "radius = " << radius << ",  area = " << setprecision(precision) << area;
  radius_file.close();
    
  // Check if radius is positive or zero
    
  if (radius >= 0)
  {
      cout << "radius is positive or zero" << endl;
  }

  if (radius < 0)
  {
      cout << "radius is negative" << endl;
  }

  return 0;			// "0" for successful completion
}

//************************** square_function **************************

double
square_function (double x)
{
    return (x * x);
}

//*********************************************************************//
