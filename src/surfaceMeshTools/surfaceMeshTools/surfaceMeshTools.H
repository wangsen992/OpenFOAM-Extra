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

Description
    Utility functions to help making surface-mesh type of operation easier.

SourceFiles
    surfaceMeshTools.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceMeshTools_H
#define surfaceMeshTools_H

#include "fvCFD.H"
#include "polyMesh.H"
#include "triSurfaceMesh.H"

namespace Foam
{
namespace surfaceMeshTools
{
    
    // Find mesh cells that are cut by surface
    labelHashSet findSurfaceCutCells
    (
        const polyMesh& mesh,
        const triSurfaceMesh& surface
    );

    // Write the triSurface representation of the surface
    void writeObjFile(OFstream, triSurface);
}
}

#endif
