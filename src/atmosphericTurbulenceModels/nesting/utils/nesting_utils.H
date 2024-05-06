#ifndef UTILS_H
#define UTILS_H

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <netcdf>
#include <gdal/ogr_spatialref.h>

#include "polyMesh.H"
#include "fvMesh.H"

// Forward declarations to avoid namespace clash
// namespace Foam
// {
//   class polyMesh;
//   class scalar;
//   class string;
//   class labelList;
// }

// TODO: Provide a command-line parser class with checks on the arguments. 

typedef std::unordered_map<std::string, std::string> string_map;
typedef std::unordered_set<std::string> string_set;
typedef std::unordered_map<std::string, std::string> str_map;

struct WRF_PROJ_PARAMS
{
  int MAP_PROJ;
  float  TRUELAT1;
  float  TRUELAT2;
  float  MOAD_CEN_LAT;
  float  STAND_LON;
  float  POLE_LAT;
  float  POLE_LON;
};

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


// format parser for command line options
string_map parse_argv(int argc, char* argv[]);

// check if the input arguments are valid for the application
int check_argv(const string_map &arguments);

// Load proj string from file and produce a SRS
std::shared_ptr<OGRSpatialReference> ReadProj(const std::string& filename);

// Quick print of the SRS (gdal tutorial example)
void printSRS(std::shared_ptr<OGRSpatialReference> psrs);

//map with predefined keys is passed to function 
void getProjAtts(const netCDF::NcFile&, WRF_PROJ_PARAMS&);

// get crs from wrf params using gdal
std::shared_ptr<OGRSpatialReference> getCRS(const WRF_PROJ_PARAMS&);

// get GCS for lat lon representation
std::shared_ptr<OGRSpatialReference> getGCS(std::string gcs_str);

// ----------- WRF ---------------- //

void readWrfCaseInfo(struct WrfCaseInfo *, netCDF::NcFile&);

Foam::autoPtr<Foam::polyMesh> meshFromNc(netCDF::NcFile&, Foam::Time&, bool);

Foam::autoPtr<Foam::fvMesh> fvmeshFromNc(netCDF::NcFile&, Foam::Time&);

Foam::tmp<Foam::volVectorField> load_U(Foam::fvMesh&, netCDF::NcFile& datafile, size_t i);

Foam::tmp<Foam::volScalarField> load_var(Foam::fvMesh&, netCDF::NcFile& datafile, const std::string& varname, size_t i);
// ----------- OpenFOAM ---------------- //
// find cells from polyMesh patch within a distance
Foam::labelList getPatchCloseCells(const Foam::polyMesh&, Foam::string, Foam::scalar);

// ----------- Interpolations ------------ //
void interp3d
(
  const std::vector<double>& x0, 
  const std::vector<double>& y0, 
  const std::vector<double>& z0,
  const std::vector<double>& u0,
  double x1,
  double y1,
  double z1,
  std::vector<double>& u1
);

#endif
