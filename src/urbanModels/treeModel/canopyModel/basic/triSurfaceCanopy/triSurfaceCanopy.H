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

Class
    Foam::triSurfaceCanopy

Description
    Surface engine to produce geometry-related parameters for canopy model

SourceFiles
    triSurfaceCanopy.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceCanopy_H
#define triSurfaceCanopy_H

#include "volFields.H"
#include "uniformGeometricFields.H"
#include "triSurfaceMesh.H"
#include "cellSetHashTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class triSurfaceCanopy Declaration
\*---------------------------------------------------------------------------*/

template<class BasicCanopy>
class triSurfaceCanopy
:
    public BasicCanopy
{
protected:

    // Protected Data

        //- Coefficients
        dictionary coeffs_;

        //- Surface file representing the canopy
        triSurfaceMesh   surface_;
    
        //- Cell Set of Canopy Influence
        labelHashSet canopyCells_;

        // Cell Set Fields

            //- Leaf Area Density [m^-1]
            dimensionedScalarCellSet lad_;

            //- Drag Coefficients
            scalarCellSet   Cd_;

    // Protected Member Functions

        static dimensionedScalarCellSet calcLAD
        (
            const polyMesh& mesh,
            triSurfaceMesh& surface,
            labelHashSet& cellsSet
        );

        dictionary coeffs();

public:

    // Constructors
    
        //- Construct from mesh
        triSurfaceCanopy
        (
          const fvMesh& mesh
        );
    
    //- Destructor
    virtual ~triSurfaceCanopy(){};
    
    // Member Functions
    
        //- Access to canopy form-related variables

            //- Cells
            virtual const labelHashSet& canopyCells() const
            {
                return canopyCells_;
            }

            //- Leaf Area Density [m^-1]
            virtual const dimensionedScalarCellSet& lad() const
            {
                return lad_;
            }

            //- Drag Coefficients
            virtual const scalarCellSet& Cd() const
            {
                return Cd_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "triSurfaceCanopy.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
