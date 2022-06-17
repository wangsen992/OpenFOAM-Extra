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
    Foam::dragCanopy

Description
    Include a basic drag model interface to the canopy base

SourceFiles
    dragCanopy.C

\*---------------------------------------------------------------------------*/

#ifndef dragCanopy_H
#define dragCanopy_H

#include "basicCanopy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class heThermo Declaration
\*---------------------------------------------------------------------------*/

class dragCanopy
:
    virtual public basicCanopy
{
public:
    
    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the implementation class
        class composite;

        //- Runtime type information
        TypeName("dragCanopy");

        //- Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            dragCanopy,
            fvMesh,
            (const fvMesh& mesh),
            (mesh)
        );
        
    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<dragCanopy> New
        (
            const fvMesh&
        );
        
    //- Destructor
    virtual ~dragCanopy(){};

    // Member Functions
    
        //- Canopy source terms to the transport equations

            //- Source term to momentum equation
            virtual dimensionedVectorCellSet& fU() = 0;

         //    //- Source term to temperature equation
         //    virtual dimensionedScalarCellSet& fT() = 0;

         //    //- Source term to species equation
         //    virtual dimensionedScalarCellSet& fY(const label specie) = 0;
         //    virtual dimensionedScalarCellSet& fY(const word& specieName) = 0;
};

/*---------------------------------------------------------------------------*\
                 Class dragCanopy::implementation Declaration
\*---------------------------------------------------------------------------*/

class dragCanopy::implementation
:
    virtual public dragCanopy
{
protected:
    
    // Protected Data 

        // Momentum Source
        dimensionedVectorCellSet fU_;

public:

    // Constructors

        //- Construct from mesh
        implementation(const fvMesh&);
    
      
    //- Destructor
    virtual ~implementation(){};
    
    // Member Functions
    
          //- Source term to momentum equation
          virtual dimensionedVectorCellSet& fU(){return fU_;};

    // Member Operators

          //- Disallow default bitwise assignment
          void operator=(const implementation&) = delete;
};

/*---------------------------------------------------------------------------*\
                 Class dragCanopy::composite Declaration
\*---------------------------------------------------------------------------*/

class dragCanopy::composite
:
    public basicCanopy::implementation,
    public dragCanopy::implementation
{
public:

    // Constructors

        //- Construct from mesh
        composite
        (
            const fvMesh& mesh
        )
        :
            basicCanopy::implementation(mesh),
            dragCanopy::implementation(mesh)
        {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
