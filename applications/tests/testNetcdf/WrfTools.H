#ifndef WrfTools_H
#define WrfTools_H

#include "fvCFD.H"
#include <netcdf>
#include <string>
#include <vector>

struct WrfCaseInfo{
      int DateStrLen;
      int N_Time, N_west_east, N_west_east_stag;
      int N_south_north, N_south_north_stag;
      int N_bottom_top, N_bottom_top_stag;

      // Load attributes
      int MAP_PROJ;
      float CEN_LAT, CEN_LON, TRUELAT1, TRUELAT2, MOAD_CEN_LAT, STAND_LON;
      float POLE_LAT, POLE_LON, DX, DY;
      int Nx;
      int Ny;
      int Nz;
      int Ncellx;
      int Ncelly;
      int Ncellz;
      std::vector<std::string> timestamps;
};


autoPtr<polyMesh> meshFromNc(netCDF::NcFile&, Time&);

void testFunc(netCDF::NcFile&);

void readWrfCaseInfo(struct WrfCaseInfo *, netCDF::NcFile&);

#endif
