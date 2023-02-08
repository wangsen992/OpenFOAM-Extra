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

#include "fileName.H"
#include "canopyTriSurfaceModel.H"
#include "dimensionedScalarFwd.H"
#include "triSurfaceSearch.H"
#include "surfaceMeshTools.H"

namespace Foam
{

template<class BaseCanopyModel>
dimensionedVectorCellSet canopyTriSurfaceModel<BaseCanopyModel>::calcLAD
(
    const polyMesh& mesh,
    const triSurface& surface,
    const labelHashSet& cellsIndex
)
{
    // Compute the total area of leaves within a mesh cell
    const pointField& surfacePoints = surface.points();
    const faceList& faces = surface.faces();
    HashTable<vector, label> cellLeafArea(cellsIndex.size());

    forAll(cellsIndex, celli)
    {
        cellLeafArea.set(cellsIndex.toc()[celli], vector(0,0,0));
    }

    Info << "calculating cellLeafArea..." << endl;
    label celli;
    point faceCenter;
    vector faceArea;
    forAll(faces, i)
    {
        faceCenter = faces[i].centre(surfacePoints);
        celli = mesh.findCell(faceCenter);
        faceArea = surface.faceAreas()[i];
        cellLeafArea[celli] += vector
                              (
                                sqrt(sqr(faceArea.x())),
                                sqrt(sqr(faceArea.y())),
                                sqrt(sqr(faceArea.z()))
                              );
    } 
    Info << "Computing leafAreaDensity..." << endl;
    HashTable<dimensionedVector, label> leafAreaDensity(cellsIndex.size());

    labelList cellList = cellLeafArea.toc();
    forAll(cellList, i)
    {
        leafAreaDensity.set
        (
          cellList[i], 
          dimensionedVector
          (
            dimArea/dimVolume, 
            cellLeafArea[cellList[i]] / mesh.cellVolumes()[cellList[i]]
          )
        );
    }

    return leafAreaDensity;
}

template<class BaseCanopyModel>
canopyTriSurfaceModel<BaseCanopyModel>::canopyTriSurfaceModel
(
    const treeModel& tree
)
:
    canopySurfaceModel<BaseCanopyModel>(tree),
    surface_(fileName(this->surfaceModelDict().lookup("file")))
{
    
    // Correct the intersected cells by the leafy surface
    labelHashSet canopyCellsIndex = surfaceMeshTools::findSurfaceCutCells(tree.mesh(), surface_);
    this->canopyCells() = canopyCellsIndex;
    
    // Compute the leaf area density
    dimensionedVectorCellSet ladCalc
    (
      calcLAD(tree.mesh(), surface_, this->canopyCells())
    );
    // set leafy area density to LAD
    this->lad() = ladCalc;
};

} 

