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
Class
    Foam::canopyGeometryModel

\*---------------------------------------------------------------------------*/

#include "canopySurfaceModel.H"
#include "dimensionedScalarFwd.H"
#include "meshSearch.H"
#include "cellClassification.H"
#include "triSurfaceSearch.H"

namespace Foam
{

labelHashSet canopySurfaceModel::findSurfaceCutCells
(
    polyMesh& mesh,
    triSurfaceMesh& surface
)
{
    
    pointField outsidePoints(1);
    outsidePoints[0] = mesh.points()[0]; // Highly unlikely this is inside
    meshSearch queryMesh(mesh);
    triSurfaceSearch querySurf(surface.surface());
    cellClassification cellType(mesh, queryMesh, querySurf, outsidePoints);
    
    DynamicList<label> dynList(0);
    forAll(cellType, celli)
    {
        if (cellType[celli] == cellClassification::CUT)
        {
            dynList.append(celli);
        }
    }
    labelHashSet cutCellsIndex(dynList);

    return cutCellsIndex;
}

HashTable<dimensionedScalar, label> canopySurfaceModel::calcLAD
(
    polyMesh& mesh,
    triSurfaceMesh& surface,
    labelHashSet& cellsIndex
)
{
    // Compute the total area of leaves within a mesh cell
    const pointField& meshPoints = mesh.points();
    const faceList& faces = surface.surface().faces();
    HashTable<scalar, label> cellLeafArea(cellsIndex.size());

    forAll(cellsIndex, celli)
    {
        cellLeafArea.set(cellsIndex[celli], 0);
    }

    label celli;
    forAll(faces, i)
    {
        celli = mesh.findCell(meshPoints[faces[i][0]]);
        cellLeafArea[celli] += mag(surface.surface().faceAreas()[i]);
    } 

    HashTable<dimensionedScalar, label> leafAreaDensity(cellsIndex.size());

    labelList cellList = cellLeafArea.toc();
    forAll(cellList, i)
    {
        leafAreaDensity.set
        (
          cellList[i], 
          dimensionedScalar
          (
            dimArea, 
            cellLeafArea[cellList[i]] 
          )
        );
        leafAreaDensity[cellList[i]] 
          = leafAreaDensity[cellList[i]] / 
                 dimensionedScalar
                 (
                    dimVolume,
                    mesh.cellVolumes()[cellList[i]]
                 );

    }

    return leafAreaDensity;
}

canopySurfaceModel::canopySurfaceModel
(
    polyMesh& mesh,
    triSurfaceMesh& surface
)
:
    mesh_(mesh),
    surface_(surface),
    canopyCellsIndex_(findSurfaceCutCells(mesh_, surface_)),
    lad_(calcLAD(mesh_, surface_, canopyCellsIndex_))
{
    
};
} 
