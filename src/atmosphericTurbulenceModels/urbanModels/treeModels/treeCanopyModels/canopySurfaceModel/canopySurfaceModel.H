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

Description
    Provides an interface for providing useful geoemtry information for
    the tree canopy. Offers both fine grained and bulk information based on
    different type of input. 

SourceFiles
    canopyGeometryModel.C
\*---------------------------------------------------------------------------*/

#ifndef canopySurfaceModel_H
#define canopySurfaceModel_H

#include "fvCFD.H"
#include "polyMesh.H"
#include "triSurfaceMesh.H"

namespace Foam
{

class canopySurfaceModel
{
    
    // Private Member Data
        
        //- Base mesh reference
        polyMesh&   mesh_;

        //- Surface file representing the canopy
        triSurfaceMesh&   surface_;

        //- cell list influenced by the canopy
        labelHashSet     canopyCellsIndex_;

        //- Leaf area density
        HashTable<dimensionedScalar, label> lad_;
         
    // Private Member Functions
        
        static labelHashSet findSurfaceCutCells
        (
            polyMesh& mesh,
            triSurfaceMesh& surface
        );

        static HashTable<dimensionedScalar, label> calcLAD
        (
            polyMesh& mesh,
            triSurfaceMesh& surface,
            labelHashSet& cellsSet
        );


public: 
    
    // Constructors
    
        // Construct from surface file and mesh
        canopySurfaceModel(polyMesh& mesh, triSurfaceMesh& surface);

    // Destructor
        ~canopySurfaceModel(){};

    // Member Functions

        // Access
            //- reference to mesh
            const polyMesh& mesh(){return mesh_;};

            //- reference to surfaceMesh
            const triSurfaceMesh& surface(){return surface_;};

            //- reference to canopy cells
            const labelHashSet& canopyCells() const {return canopyCellsIndex_;};

        // Parameters for all types of analysis

            //- reference to leaf area density
            const HashTable<dimensionedScalar, label>& lad() const {return lad_;};

};
}

#endif 
