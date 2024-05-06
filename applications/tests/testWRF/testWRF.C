#include <string>
#include <fstream>
#include "fvCFD.H"
#include "WRF.H"

using namespace std;
using namespace Foam;
int main(int argc, char **argv)
{
  // Load openfoam case
  std::string case_name("foamcase/system/controlDict"); 

  argList args(argc, argv);
  Time runTime(Time::controlDictName, args);

  fvMesh mesh
  (
    IOobject
    (
      polyMesh::defaultRegion,
      runTime.timeName(),
      runTime,
      IOobject::MUST_READ
    )
  );
  Info << "foamCase loaded." << endl;

  // Set up WRF case
  ifstream proj_file("proj4_neiu");
  std::string proj_string;
  std::getline(proj_file, proj_string);
  Foam::Time wrfTime("/app/gdal-test", "wrfCase");
  WRF wrf("ncfile", "./", "wrfCase", runTime, proj_string);
  Info << "WRF completed" << endl;

  Info << wrf.mesh().cells().size() << endl;

  return 0;
}
