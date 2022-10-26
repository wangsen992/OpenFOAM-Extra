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
    Foam::trunkSurfaceModel

Description
    Base class for the handling the tree trunk surface geometry information
    
SourceFiles
    trunkSurfaceModel.C
\*---------------------------------------------------------------------------*/

#ifndef trunkSurfaceModel_H
#define trunkSurfaceModel_H

#include "fvCFD.H"
#include "polyMesh.H"
#include "triSurfaceMesh.H"

namespace Foam
{

class trunkSurfaceModel
{
    // Private Member Data
    
        //- Base mesh reference
        polyMesh&   mesh_;

        //- Surface file representing the canopy
        triSurfaceMesh&   surface_;

        //- cell list influenced by the canopy
        labelHashSet     canopyCellsIndex_;
    // Private Member Functions

        //- Compute cells that are cut by the surface
        
public:
    
    // Constructors
        
        //- Construct from component
        trunkSurfaceModel
        (
          polyMesh& mesh,
          triSurfaceMesh& surface
        );

    // Destructor
    

    // Member Functions
        // Access

        // Parameters

};

}

#endif
