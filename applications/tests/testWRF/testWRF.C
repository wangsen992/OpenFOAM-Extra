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

  // Testing get var function
  {
    Info<< average(wrf.U(1)) << endl;
    Info<< average(wrf.U(2)) << endl;
    Info<< average(wrf.U(3)) << endl;
  }
  {
    Info<< average(wrf.var("T", 1)) << endl;
    Info<< average(wrf.var("T", 2)) << endl;
    Info<< average(wrf.var("T", 3)) << endl;
  }

  // Testing interpolation of close-to-patch cells
  const pointField& foamCellCentres(mesh.cellCentres());
  labelList patchCellInd(getPatchCloseCells(mesh, "inlet", 20));
  scalarField Tinterp(patchCellInd.size());
  pointField patchCellCentres(patchCellInd.size());
  Info << "Running std::transform" << endl;
  std::transform
  (
    patchCellInd.cbegin(), 
    patchCellInd.cend(),
    patchCellCentres.begin(),
    [&](label i) { return wrf.transform(foamCellCentres[i]); }
  );
  Info << wrf.interpolate<Foam::vector>(patchCellCentres, wrf.U(2).ref()) << endl;;

  return 0;
}
