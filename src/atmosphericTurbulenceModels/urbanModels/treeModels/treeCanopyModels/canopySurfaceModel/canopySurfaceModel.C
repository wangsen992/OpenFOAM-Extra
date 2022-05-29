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
    Info << "running findSurfaceCutCells" << endl; 
    DynamicList<label> dynList(0);

    // Compute the total area of leaves within a mesh cell
    const pointField& surfacePoints = surface.points();
    const faceList& faces = surface.surface().faces();

    label celli;
    point faceCenter;
    forAll(faces, i)
    {
        
        faceCenter = faces[i].centre(surfacePoints);
        celli = mesh.findCell(faceCenter);
        dynList.append(celli);
    } 

    return labelHashSet(dynList);
}

HashTable<dimensionedScalar, label> canopySurfaceModel::calcLAD
(
    polyMesh& mesh,
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
