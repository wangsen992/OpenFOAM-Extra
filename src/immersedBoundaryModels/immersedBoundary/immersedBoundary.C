/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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

Class
    Foam::fv::immersedBoundary

\*---------------------------------------------------------------------------*/

#include "immersedBoundary.H"
#include "cellClassification.H"
#include "interpolationCellPoint.H"
// * * * * * * * * * * * * Static Member Functions  * * * * * * * * * * * * * //
Foam::labelHashSet Foam::immersedBoundary::calcForcingCells
(
  const polyMesh& mesh,
  const triSurfaceMesh& surf
)
{
    pointField outsidePoints(1);
    outsidePoints[0] = mesh.points()[0];
    meshSearch queryMesh(mesh);
    triSurfaceSearch querySurf(surf.surface());
    cellClassification cellType(mesh, queryMesh, querySurf, outsidePoints);
    
    DynamicList<label> dynList(0);
    forAll(cellType, celli)
    {
        if (cellType[celli] == cellClassification::CUT)
        {
            dynList.append(celli);
        }
    }
    // Obtain the cells that intersects with the surface to prepare for
    // immersed boundary method implementation
    labelHashSet forcingCellsHashSet(dynList.size());
    forAll(dynList, i)
    {
         forcingCellsHashSet.insert(dynList[i]);
    }
    
    return forcingCellsHashSet;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
void Foam::immersedBoundary::calcForcingPointAndNormal
(
  const polyMesh& mesh, 
  const triSurface& surf, 
  const labelHashSet& forcingCells,
  pointField& forcingPoints,
  vectorField& forcingPointNormals
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
    // Forcing point here is computed using a crude averaging method
    // over the projection of cell vertices onto the surface instead 
    // of the actual center of gravity
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

void Foam::immersedBoundary::calc()
{
    // Use meshSearch and surfSearch and cellClassification to identify cells
    // that are cut by the surface. 
    calcForcingPointAndNormal
    (
      mesh_, 
      surf_, 
      forcingCells_,
      forcingPoints_, 
      forcingPointNormals_
    );

    scalar dpp(0.5);
    probePoints_ = forcingPoints_ + dpp * forcingPointNormals_;
    
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::immersedBoundary::immersedBoundary
(
    const polyMesh& mesh,
    const triSurfaceMesh& surf
)
:
    mesh_(mesh),
    surf_(surf),
    forcingCells_
    (
      mesh_,
      "forcingCells",
      calcForcingCells(mesh_, surf_),
      IOobject::NO_WRITE
    ),
    surfaceFaceCenters_(surf_.surface().faceCentres()),
    surfaceFaceNormals_(surf_.surface().faceNormals()),
    forcingPoints_(forcingCells_.size()),
    forcingPointNormals_(forcingCells_.size()),
    probePoints_(forcingCells_.size())
{
    calc();
}
    

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::volVectorField> fFp()
{}

Foam::tmp<Foam::volVectorField> Foam::immersedBoundary::fnFp(Foam::volVectorField& U)
{
    // Interpolate velocity to the forcing point
    // Try "interpolationCellPoint.H"
    interpolationCellPoint<Foam::vector> interpCell(U);
    labelField forcingCellLabels(forcingCells_.toc());
    vectorField Ufp
    (
      interpCell.interpolate
      (
        forcingPoints_,
        forcingCellLabels
      ) 
      & forcingPointNormals_ * forcingPointNormals_ // normalize to normal dir
    );

    // Compute the difference between normal velocity and desired value
    // Note there is a requirement to store previous values
    dimensionedVector Utarget(dimVelocity, vector(0, 0, 0));

}

Foam::tmp<Foam::volVectorField> fpFp()
{}
