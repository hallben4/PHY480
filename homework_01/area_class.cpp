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
#include <math.h>        // in order to use atan
using namespace std;     // if omitted, then need std::cout, std::cin

//*********************************************************************//

double pi = 4 * atan (1);   // define pi in terms of atan

class Circle {
    int radius;
public:
    void set_values (double);
    double area() {return pi * radius * radius;}
};

void Circle::set_values (double x) {
    radius = x;}

int main () {
    
    double rad;
    
    cout << "Enter the radius of a circle: ";    // ask for radius
    cin >> rad;
    
    Circle circ;
    circ.set_values (rad);
    cout << "area: " << circ.area();
    return 0;
}

//*********************************************************************//
