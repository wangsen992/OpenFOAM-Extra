#include <string>
#include <fstream>
#include "fvCFD.H"
#include "WRF.H"
#include <random>
#include <functional>
#include "indexedOctree.H"
#include "treeDataFace.H"
#include "boundBox.H"

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

  // Load WRF terrain
  std::string proj4(" +proj=lcc +lat_1=30 +lat_2=60 +lat_0=41.980352 +lon_0=-87.717729");
  WRF wrf("ncfile", "./", "wrfCase", runTime, proj4);
  fvMesh& wrfMesh(wrf.mesh());

  // transform the mesh terrain 
  runTime++;
  Info << "Time: " << runTime.timeName() << endl;
  Info << "wrf Bound: " << boundBox(wrf.mesh().points()) << endl;
  Info << boundBox(mesh.points()) << endl;
  const polyPatch& wrfGroundPatch(wrfMesh.boundaryMesh()["bottom"]);
  pointField groundPoints = wrfGroundPatch.points();

  indexedOctree<treeDataFace> wrfTree
  (
    treeDataFace(false, wrfGroundPatch),
    treeBoundBox(boundBox(wrfGroundPatch.points())),
    10,
    10,
    3
  );

  // Interp the ground points from openfoam to wrf bottom patch
  const polyPatch& foamGroundPatch(mesh.boundaryMesh()["ground"]);

  indexedOctree<treeDataFace> foamTree
  (
    treeDataFace(false, foamGroundPatch),
    treeBoundBox(boundBox(foamGroundPatch.points())),
    10,
    10,
    3
  );
  Foam::vector ll{0,0,3000};
  Foam::vector zvec{0,0,1};
  auto findVec = [&](const point& pt)
  {
    auto wrfHit =  wrfTree.findLine(pt-ll, pt+ll);
    point wrfPt = wrfHit.hitPoint(); 
    auto foamHit = foamTree.findLine(pt-ll, pt+ll);
    point foamPt;
    if (!foamHit.hit())
    {
      foamPt = foamTree.findNearest(pt, 1e9).hitPoint();
    }
    else
    {
      foamPt = foamHit.hitPoint(); 
    }

    return zvec * (wrfPt.z() - foamPt.z());
    // compress the high altitude regions
  };

  {
    pointField foamPts = mesh.points();
    vectorField vec(foamPts.size());
    std::transform
    (
      foamPts.cbegin(), 
      foamPts.cend(), 
      vec.begin(),
      findVec
    );

    // assuming foam mesh is always lower than wrf mesh
    scalar zmax = max(foamPts.component(2));
    scalar zmin = min(foamPts.component(2));
    scalar vec_zmax(max(vec.component(2)));
    scalar vec_zmin(min(vec.component(2)));
    Info << vec_zmax - vec_zmin << endl;

    vec = Foam::vector{0,0,vec_zmin} 
        +(
            (zmax - foamPts.component(2))/(zmax-zmin)
           *(vec - Foam::vector{0,0,vec_zmin})
         );

    mesh.movePoints(foamPts + vec);
    mesh.setInstance(runTime.timeName());
    mesh.write();
  }


  return 0;
}
