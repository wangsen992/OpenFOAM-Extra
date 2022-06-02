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
    Foam::meshSurfaceCutter

Description
    Helper class that cuts a surface (comprises of faces) by mesh cells 
    & faces. 

SourceFiles
    meshSurfaceCutter.C
\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceCutter_H
#define meshSurfaceCutter_H

#include "polyMesh.H"
#include "triSurfaceMesh.H"

namespace Foam
{

namespace surfaceMeshTools
{

class meshSurfaceCutter
{
    // Private Member Data

        //- Reference to mesh and surface
        const polyMesh&     mesh_;
        const triSurfaceMesh&     surf_;

public:
    
    // Constructors

        //- Construt from components
        meshSurfaceCutter
        (
            const polyMesh& mesh, 
            const triSurfaceMesh& surf
        );

    // Destructor
        
        ~meshSurfaceCutter(){};

    // Public Member Functions
        
        //- Cut surface with mesh
        tmp<triSurfaceMesh> cutSurf();

};
}
}

#endif
