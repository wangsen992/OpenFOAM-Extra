#ifndef LambertConverter_H
#define LambertConverter_H

#include <iostream>
#include <math.h>
#include <utility>
#include <vector>
using namespace std;

namespace mapProj
{
  class LambertConverter
  {
      
      const float PI=3.1415927;
      const float RAD_PER_DEG=PI/180.0;
      const float DEG_PER_RAD=180.0/PI;
      
      const float stdlon;
      const float truelat1, truelat2;
      const float re_m;
      const int hemi;
      const float dx, dy;
      const float lat1, lon1;
      const float x1, y1;

  public:

      // Construct with every param
      LambertConverter
      (
        float stdlon, 
        float truelat1,
        float truelat2,
        float re_m,
        float hemi,
        float dx,
        float dy,
        float lat1,
        float lon1,
        float x1,
        float y1
      );

      // Default constructor
      LambertConverter
      ();

      ~LambertConverter(){};

      // Methods
        // Compute one pair of x,y
        void get_xy
        (
          const float lat,
          const float lon, 
          float& x,
          float& y
        );

        void get_xy
        (
          const vector<float> lat,
          const vector<float> lon,
          vector<float>& x,
          vector<float>& y
        );
  };

}

#endif
