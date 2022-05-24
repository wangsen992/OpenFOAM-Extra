/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    testMesh

Description
    Refer to Troldborg et al., 2021 for the nomenclature

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "DynamicList.H"
#include "Ostream.H"
#include "fvCFD.H"
#include "fileName.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"
#include "searchableSurface.H"
#include "triSurfaceSearch.H"
#include "cellZone.H"
#include "cellZoneSet.H"
#include "cellClassification.H"
#include "meshSearch.H"
#include "pointIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// cellSet surfaceToCellSet
// (
//     fvMesh& mesh,
//     triSurface& surf
// )
// {
//     meshSearch queryMesh(mesh);
//     triSurfaceSearch querySurf(surf);
//     pointField outsidePoints(1);
//     outsidePoints[0] = mesh.points()[0];
//     cellClassification cellType
//     (
//       mesh,
//       queryMesh,
//       querySurf,
//       outsidePoints
//     );
// 
//     label cutCount = 0;
//     forAll(cellType, celli)
//     {
//         if (cellType[celli] == cellClassification::CUT)
//         {
//             cutCount += 1;
//         }
//     }
//     Info << "Total number of cutted cells: " 
//       << cutCount << endl;
// 
//     // Obtain the cells that intersects with the surface to prepare for
//     // immersed boundary method implementation
//     labelHashSet cutCells(cutCount);
//     label icut = 0;
//     forAll(cellType, celli)
//     {
//         if (cellType[celli] == cellClassification::CUT)
//         {
//            cutCells.insert(celli);
//            icut += 1;
//         }
//     }
//     cellSet forcingCells(mesh, "forcingCells", cutCells, IOobject::AUTO_WRITE);
// };

//- Extract forcing point positions from the mesh and surface using an
//arbitrary method by finding the surface projection of the average of all
//projections of cell vertices
void calcForcingPointAndNormal
(
  fvMesh& mesh, 
  triSurface& surf, 
  labelHashSet& forcingCells,
  pointIOField& forcingPoints,
  vectorIOField& forcingPointNormals
)
{
  // Create a pointField to store points
  triSurfaceSearch querySurf(surf);

  // Iterate over the cells to compute forcing points
  List<PointIndexHit<point>> pointHits(8);
  List<point> plist(8);
  label celli;
  face faceHit;
  PointIndexHit<point> pointHit;

  forAll(forcingCells.toc(), i)
  {
    celli = forcingCells.toc()[i];
    querySurf.findNearest
    (
       mesh.cells()[celli].points(mesh.faces(), mesh.points()),
       scalarField(List<scalar>(label(8), scalar(1))),
       pointHits
    );
    forAll(pointHits, hiti)
    {
        plist[hiti] = pointHits[hiti].hitPoint();
    }
    pointHit = querySurf.nearest(average(plist), vector(1,1,1));
    forcingPoints[i] = pointHit.hitPoint();
    faceHit = surf.faces()[pointHit.index()];

    forcingPointNormals[i] = normalised(face::area<pointField>(faceHit.points(surf.points())));
  }
}

int main(int argc, char *argv[])
{
  
    argList::addOption("filepath", "/path/to/file");
    argList args(argc, argv);
    Time runTime(Time::controlDictName, args);
    fvMesh mesh
    (
        IOobject
        (
          fvMesh::defaultRegion,
          runTime.constant(),
          runTime,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
        )
    );
    Info << "mesh info" << nl
      << "mesh points number: " << mesh.points().size() << endl;

    fileName f(args.optionRead<string>("filepath"));
    Info << "filepath: " << f << endl;

    // Info << "Loading the triSurface file..." << endl;
    // Load triSurface
    // triSurface surf(f);
    triSurfaceMesh surf
    (
        IOobject
        (
            f.name(),
            runTime.constant()/"geometry",
            mesh,
            IOobject::MUST_READ
        )
    );
    Info << "Instantiation done" << endl;
    Info << "triSurface loaded: " << nl
          << "Number of points: " << surf.points()().size() << endl;
    

    // Use meshSearch and surfSearch and cellClassification to identify cells
    // that are cut by the surface. 
    pointField outsidePoints(1);
    outsidePoints[0] = point( 1, 1, 1);
    meshSearch queryMesh(mesh);
    triSurfaceSearch querySurf(surf.surface());
    cellClassification cellType(mesh, queryMesh, querySurf, outsidePoints);

    label cutCount = 0;
    forAll(cellType, celli)
    {
        if (cellType[celli] == cellClassification::CUT)
        {
            cutCount += 1;
        }
    }
    Info << "Total number of cutted cells: " 
      << cutCount << endl;

    // Obtain the cells that intersects with the surface to prepare for
    // immersed boundary method implementation
    labelHashSet cutCells(cutCount);
    label icut = 0;
    forAll(cellType, celli)
    {
        if (cellType[celli] == cellClassification::CUT)
        {
           cutCells.insert(celli);
           icut += 1;
        }
    }
    cellSet forcingCells(mesh, "forcingCells", cutCells, IOobject::AUTO_WRITE);
    // cellZoneSet cutCellZoneSet(mesh, "cutCells", forcingCells, IOobject::AUTO_WRITE);
    
    // Working with meshCellZones
    meshCellZones& cellZones(mesh.cellZones());
    label sz = cellZones.size();
    cellZones.setSize(sz+1);
    cellZones.set
    (
      sz,
      new cellZone
      (
          "forcingCells",
          forcingCells.toc(),
          sz,
          cellZones
      )
    );
    cellZones.writeOpt() = IOobject::AUTO_WRITE;
    cellZones.instance() = mesh.facesInstance();
    
    pointIOField sc
    (
      IOobject
      (
          "sc",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      surf.surface().faceCentres()
    );
    pointIOField surfFaceCtrsNorm
    (
      IOobject
      (
          "nsc",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      surf.surface().faceNormals()
    );

    pointIOField forcingPts
    (
      IOobject
      (
          "fp",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      pointField(List<point>(cutCells.toc().size()))
    );

    vectorIOField forcingPtNormal
    (
      IOobject
      (
          "nfp",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      vectorField(List<vector>(cutCells.toc().size()))
    );
    
    calcForcingPointAndNormal(mesh, surf, cutCells, forcingPts, forcingPtNormal);

    scalar dpp(1);
    vectorIOField probePoints
    (
      IOobject
      (
          "pp",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      forcingPts + dpp * forcingPtNormal
    );

    runTime++;
    runTime.write();

    Info<< "End\n" << endl;
  
  return 0;
}
