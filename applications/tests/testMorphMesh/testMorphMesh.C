#include <string>
#include <fstream>
#include "fvCFD.H"
#include "WRF.H"
#include <random>
#include <functional>

using namespace std;
using namespace Foam;
int main(int argc, char **argv)
{
  argList::addOption("test", "test name", "test only getVar");
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
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
    )
  );
  Info << "foamCase loaded." << endl;


  scalar amp(50);
  scalar zmax(max(mesh.cellCentres().component(2)));
  scalar zmin(min(mesh.cellCentres().component(2)));

  pointField pts(mesh.points());
  pointField moveVec(pts.size());

  std::default_random_engine random_engine;
  std::normal_distribution<scalar> dist(1, 0.2);
  auto my_dist=std::bind(dist, random_engine);
  while(runTime.loop())
  {
    runTime++;
    Info << "Time: " << runTime.timeName() << endl;
    scalar t = runTime.value();
    std::transform(pts.begin(), pts.end(), moveVec.begin(), 
                  [&](point pt)
                  {
                    return point(0, 0, 
                        ((zmax-pt.z())/ zmax)*amp
                       *pow(Foam::sin(6.28*(pt.y()+pt.x())/200 + 6.28*t/50),2));});

    mesh.movePoints(moveVec+pts);
    mesh.setInstance(runTime.timeName());
        
    mesh.write();
  }

  return 0;
}
