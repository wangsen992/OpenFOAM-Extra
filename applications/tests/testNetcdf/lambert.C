// Conversion test for lambert projection. 
// ref: https://mathworld.wolfram.com/LambertConformalConicProjection.html
#include <iostream>
#include <math.h>
#include "LambertConverter.H"
using namespace std;

int main(int argc, char *argv[])
{
    std::cout << "Hello World Yo" << std::endl;

    // ll to xy
    // Loading attributes
    // set target lat/lon
    float lat=41.0, lon=-89.0;

    mapProj::LambertConverter lam;

    float i, j;

    lam.get_xy(lat, lon, i, j);


    std::cout << "i, j : " << i << ", " << j << std::endl;
    std::cout << "d : " << sqrt(i*i + j*j) << std::endl;

    return 0;
}
