
#include "LambertConverter.H"

mapProj::LambertConverter::LambertConverter
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
)
:
  stdlon(stdlon),
  truelat1(truelat1),
  truelat2(truelat2),
  re_m(re_m),
  hemi(hemi),
  dx(dx),
  dy(dy),
  lat1(lat1),
  lon1(lon1),
  x1(x1),
  y1(y1)
{}

mapProj::LambertConverter::LambertConverter
()
:
  stdlon(-87.6),
  truelat1(30.0),
  truelat2(60.0),
  re_m(6378137.0),
  hemi(1),
  dx(1.0),
  dy(1.0),
  lat1(40.9443),
  lon1(-89.3194),
  x1(0.0),
  y1(0.0)
{}


void mapProj::LambertConverter::get_xy
(
  const float lat,
  const float lon, 
  float& x,
  float& y
)
{
    float deltalon1;
    float tl1r, ctl1r;
    float rsw, deltalon;
    float rm, polei, polej;
    float rebydx;
    float arg, cone;
    
    // start conversion
    rebydx = re_m / dx;
    // Apply lambert
    
    // Ignored two checks (truelat1 != truelat2 and wide apart)
    std::cout << "Debug Line" << std::endl;
    std::cout << truelat1 << std::endl;
    std::cout << truelat2 << std::endl;
    std::cout << RAD_PER_DEG << std::endl;
    cone = (log(cos(truelat1*RAD_PER_DEG))-log(cos(truelat2*RAD_PER_DEG))) / 
          (log(tan((90.0 - abs(truelat1))*RAD_PER_DEG*0.5)) - 
           log(tan((90.0 - abs(truelat2))*RAD_PER_DEG*0.5)));

    // Compute lon differences and ensure we stay out of the forbidden cut zone
    deltalon1 = lon1 - stdlon;
    if (deltalon1 > 180.0)
    {
        deltalon1 = deltalon1 - 360.0;
    }
    else if (deltalon1 < -180.0) {
        deltalon1 = deltalon1 + 360.0;
    }

    // Convert truelat1 to radian and compute cos for later use
    tl1r = truelat1 * RAD_PER_DEG;
    ctl1r = cos(tl1r);

    // Compute the radius to our known lower-left corner
    rsw = rebydx * ctl1r/cone*
        pow
        (
          tan((90.0*hemi-lat1)*RAD_PER_DEG/2.0) / 
          tan((90.0*hemi - truelat1)*RAD_PER_DEG/2.0), 
          cone
        );

    // Find pole point
    arg = cone*(deltalon1*RAD_PER_DEG);
    polei = hemi*x1 - hemi*rsw*sin(arg);
    polej = hemi*y1 + rsw*cos(arg);
            
    // compute deltalon between known longitude and startdard lon and ensure
    // it is not in the cut zone
    deltalon = lon - stdlon;
    if (deltalon > 180.0)
    {
        deltalon = deltalon - 360.0;
    }
    else if (deltalon < -180.0) {
        deltalon = deltalon + 360.0;
    }

    // radius to desired point
    rm = rebydx*ctl1r/cone*
        pow
        (
          tan((90.0*hemi - lat)*RAD_PER_DEG/2.0) / 
          tan((90.0*hemi - truelat1)*RAD_PER_DEG/2.0),
          cone
        );
    arg = cone*(deltalon*RAD_PER_DEG);
    x = polei + hemi*rm*sin(arg);
    y = polej - rm*cos(arg);

    x = hemi * x;
    y = hemi * y;
  
}

void mapProj::LambertConverter::get_xy
(
  const vector<float> lat,
  const vector<float> lon,
  vector<float>& x,
  vector<float>& y
)
{
    for(int i = 0; i < lat.size(); i++)
    {
        get_xy(lat[i], lon[i], x[i], y[i]);
    }
}
