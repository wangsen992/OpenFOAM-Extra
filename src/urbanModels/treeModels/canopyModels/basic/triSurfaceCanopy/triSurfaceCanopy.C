/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "triSurfaceCanopy.H"
#include "meshSearch.H"
#include "cellClassification.H"
#include "surfaceMeshTools.H"

namespace Foam
{
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicCanopy>
dimensionedScalarCellSet triSurfaceCanopy<BasicCanopy>::calcLAD
(
    const polyMesh& mesh,
    triSurfaceMesh& surface,
    labelHashSet& cellsIndex
)
{
    // Compute the total area of leaves within a mesh cell
    const pointField& surfacePoints = surface.points();
    const faceList& faces = surface.surface().faces();
    HashTable<scalar, label> cellLeafArea(cellsIndex.size());

    forAll(cellsIndex, celli)
    {
        // Info << "cellIndex[celli] = " << cellsIndex.toc()[celli] << endl;
        cellLeafArea.set(cellsIndex.toc()[celli], 0);
    }

    Info << "calculating cellLeafArea..." << endl;
    label celli;
    point faceCenter;
    forAll(faces, i)
    {
        
        faceCenter = faces[i].centre(surfacePoints);
        celli = mesh.findCell(faceCenter);
        cellLeafArea[celli] += mag(surface.surface().faceAreas()[i]);
    } 
    
    Info << "Computing leafAreaDensity..." << endl;

    HashTable<dimensionedScalar, label> leafAreaDensity(cellsIndex.size());

    labelList cellList = cellLeafArea.toc();
    forAll(cellList, i)
    {
        leafAreaDensity.set
        (
          cellList[i], 
          dimensionedScalar
          (
            dimArea/dimVolume, 
            cellLeafArea[cellList[i]] / mesh.cellVolumes()[cellList[i]]
          )
        );
    }

    return leafAreaDensity;
}

template<class BasicCanopy>
dictionary triSurfaceCanopy<BasicCanopy>::coeffs()
{
    word subDictName = "triSurfaceCanopyCoeffs";
    const dictionary& canopyDict(BasicCanopy::properties());
    dictionary triSurfaceCanopyDict(canopyDict.subDict(subDictName));
    return triSurfaceCanopyDict;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class BasicCanopy>
triSurfaceCanopy<BasicCanopy>::triSurfaceCanopy
(
    const fvMesh& mesh
)
:
    BasicCanopy(mesh),
    coeffs_(coeffs()),
    surface_
    (
        IOobject
        (
            coeffs_.lookup("surface"),
            mesh.time().constant()/"geometry",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    canopyCells_
    (
      surfaceMeshTools::findSurfaceCutCells(mesh, surface_)
    ),
    lad_(calcLAD(mesh, surface_, canopyCells_)),
    Cd_(canopyCells_.toc().size())
{
    // Arbitrary setting Cd_
    scalar Cd = coeffs_.lookupOrAddDefault("Cd", 0.5);
    labelList cells = canopyCells_.sortedToc();
    forAll(cells, i)
    {
        Cd_.set(cells[i], Cd);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

}
// ************************************************************************* //

