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
#include "surfaceMeshTools.H"
#include "dimensionSet.H"

// * * * * * * * * * * * * Static Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
void Foam::immersedBoundary::calcForcingPointAndNormal
(
  const polyMesh& mesh, 
  const triSurface& surf, 
  const cellSet& forcingCells,
  vectorCellSet& forcingPoints,
  vectorCellSet& forcingPointNormals
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
  labelList forcingCellList(forcingCells.sortedToc());

  forAll(forcingCellList, i)
  {
    celli = forcingCellList[i];
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
    forcingPoints[celli] = pointHit.hitPoint();
    faceHit = surf.faces()[pointHit.index()];

    forcingPointNormals[celli] = normalised(face::area<pointField>(faceHit.points(surf.points())));
  }
}

void Foam::immersedBoundary::readCoeffs()
{}

void Foam::immersedBoundary::init()
{
    // Init surface face properties
    forAll(surfaceFaceCenters_, i)
    {
        surfaceFaceCenters_[i] = surf_.surface().faceCentres()[i];
        surfaceFaceNormals_[i] = surf_.surface().faceNormals()[i];
    }
    
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
    label celli;
    labelList forcingCellsList(forcingCells_.sortedToc());
    forAll(forcingCellsList, i)
    {
        celli = forcingCellsList[i];
        probePoints_[celli] = forcingPoints_[celli] + dpp * forcingPointNormals_[celli];
    }

    readCoeffs();
    
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
    dpp_(0.5),
    forcingCells_
    (
      mesh_,
      "forcingCells",
      surfaceMeshTools::findSurfaceCutCells(mesh_, surf_),
      IOobject::NO_WRITE
    ),
    surfaceFaceCenters_(surf_.surface().faceCentres().size()),
    surfaceFaceNormals_(surf_.surface().faceNormals().size()),
    forcingPoints_(forcingCells_.size()),
    forcingPointAreas_(forcingCells_.size()), // This is not calculated
    forcingPointNormals_(forcingCells_.size()),
    probePoints_(forcingCells_.size()),
    Kp_(dimensionSet(1, 0, -1, 0, 0), 10),
    KI_(dimensionSet(1, 0, -1, 0, 0), 10),
    ustar_(dimVelocity, 0.1),
    fU_(forcingCells_.toc().size())
{
    init();
}
    

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<dimensionedVectorCellSet> Foam::immersedBoundary::fnFp(Foam::volVectorField& U)
{
    // Interpolate velocity to the forcing point
    // Try "interpolationCellPoint.H"
    interpolationCellPoint<Foam::vector> interpCell(U);
    labelField forcingCellLabels(forcingCells_.sortedToc());
    pointField forcingPointField = vectorCellSet2pointField(forcingPoints_);
    pointField forcingPointNormalField = vectorCellSet2pointField(forcingPointNormals_);
    
    tmp<vectorField> tUfp
    (
      interpCell.interpolate
      (
        forcingPointField,
        forcingCellLabels
      ) 
      & forcingPointNormalField * forcingPointNormalField // normalize to normal dir
    );
    

    // Compute the difference between normal velocity and desired value
    // Note there is a requirement to store previous values
    // Current version has no previus iteration values
    vector Utarget(0, 0, 0);

    tmp<dimensionedVectorCellSet> tfnfp
    (
        new dimensionedVectorCellSet(forcingCells_.toc().size())
    );

    label celli;
    forAll(forcingCellLabels, i)
    {
        celli = forcingCellLabels[i];
        tfnfp.ref()[celli] = Kp_ * dimensionedVector(dimVelocity, (tUfp()[i] - Utarget));
    }

    return tfnfp;
}

tmp<dimensionedVectorCellSet> Foam::immersedBoundary::fFp(Foam::volVectorField& U)
{

}
