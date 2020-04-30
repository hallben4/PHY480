//************************************************************************ 
//
//  Program to compare the difference between summing up and summing down
//
//  Programmer:  Ben Hall
//
//  Answer to part (c): For N < 3, the error between the two methods is extremely small (the calculations are basically equally precise) and so at this scale it doesn't matter to much whether you sum up or down. The, right after  N=4, the error starts to grow rapidly and so which method one uses matters more and more. It looks like the region of 4 <= log(N) <= 6 is a power law with approximate slope 2. This implies that in this region, the error scales qudratically with N. For N > 7, the error starts to asymptotically approach a constant, the slope goes to 1. Upward sum is more precise because with downward sum, one stars with a really large number and then adds smaller and smaller numbers so that by the end, one is adding very large numbers to very small numbers, which introduces error.
//
//*************************************************************************


// include files
#include <iostream>		// note that .h is omitted
#include <iomanip>		// note that .h is omitted
#include <fstream>
#include <cmath>
using namespace std;		// we need this when .h is omitted

//********************** begin main ******************************

// inline functions

inline float sum_up (int N)
{
    // initialize sum to zero
    float sum = 0;
    
    float num = 1;
    
    // for loop to do sum
    for (int n = 1; n <= N; n++)
    {
        sum += num/n;
    }
    
    return sum;
}

inline float sum_down (int N)
{
    // initialize sum to zero
    float sum = 0;
    
    float num = 1;
    
    // for loop to do sum
    for (int n = 1; n <= N; n++)
    {
        sum += num/(N-n+1);

    }
    
    return sum;
}

inline float sum_error (int N)
{
    
    return (0.5) * abs(sum_up(N) - sum_down(N)) / (abs(sum_up(N)) + abs(sum_down(N)));
}

int main ()
{
    // open the output file
    ofstream sum_error_file ("summing.dat");
    
    // define N_max
    double N_max = pow(10,9);
    
    // loop through N
    for (int N = 1; N <= N_max; N *= 2)
    {
        // Output to file
        sum_error_file << setprecision(18) << log10(N) << "   " << log10(sum_error(N)) << endl;
    }
    
    return (0);
}


