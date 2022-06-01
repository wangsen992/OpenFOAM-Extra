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

\*---------------------------------------------------------------------------*/

#include "surfaceMeshTools.H"

namespace Foam
{

labelHashSet surfaceMeshTools::findSurfaceCutCells
(
    const polyMesh &mesh, 
    const triSurfaceMesh &surface
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
}
