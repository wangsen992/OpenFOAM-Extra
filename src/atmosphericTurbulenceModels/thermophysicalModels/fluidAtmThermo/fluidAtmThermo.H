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
    Foam::fluidAtmThermo

Description
    A fluid thermo extended for atm (supporting multi-component mixture). 
    This branch of thermo module enables the use of potential
    temperature (or entropy conserved temperature variable) for energy transport 
    computation by providing variables attributes and interface relevant for 
    atmospheric science problems. 

    The final application can still be written in either with potential 
    temperature or with internal energy form. 

    pRef and rhoRef are assigned after object construction and hydrostatic
    initialization, and used for buoyancy definition. 

    The key difference for atmThermo is that due to the presence of gravity,
    there is a hydrostatic balance so definition of buoyancy becomes important. 

\*---------------------------------------------------------------------------*/

#ifndef fluidAtmThermo_H
#define fluidAtmThermo_H

#include "fvCFD.H"
#include "fluidReactionThermo.H"
#include "speciesTable.H"

typedef fluidReactionThermo BaseFluidThermo;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class fluidReactionThermo Declaration
\*---------------------------------------------------------------------------*/
class fluidAtmThermo
:
    virtual public BaseFluidThermo
{
public:
    
    // Public Classes
    
        //- Forward declare the implementation class
        class implementation;
    
    //- Runtime type information
    TypeName("fluidAtmThermo");

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidAtmThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<fluidAtmThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );

    //- Destructor
    virtual ~fluidAtmThermo();
    
    // Member Functions

        // Access to thermodynamic state variables

            //- Reference Pressure [Pa]
            virtual volScalarField& pRef ()= 0;
            
            //- Reference Pressure [Pa]
            virtual const volScalarField& pRef() const = 0;
            
            //- Reference Density [Pa]
            //  Non-const access allowed for assignment
            virtual volScalarField& rhoRef() = 0;
            
            //- Reference Density [Pa]
            virtual const volScalarField& rhoRef() const = 0;
            
        // Access to Mixture Specie Information
        
        // Derived variables from available references (see above) 

            //- Buoyancy [kg m-2 s-2]
            // Note: Buoyancy is defined in its most fundamental form
            // as the product of density variation with gravity
            virtual const volVectorField& b() const =  0;
            virtual volVectorField& b() =  0;
            
};
/*---------------------------------------------------------------------------*\
                     Class fluidReactionThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidAtmThermo::implementation
:
    virtual public fluidAtmThermo
{
protected:
    // Protected data

        // Model Parameters

            //- Reference Pressure [Pa]
            volScalarField  pRef_;

            //- Gravitational Acceleration
            uniformDimensionedVectorField g_;

public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const fvMesh&, const word& phaseName);

    //- Destructor
    virtual ~implementation();
    
    // Member Functions

        // Access to thermodynamic state variables

            //- Reference Pressure [K]
            //  Non-const access allowed for transport equations
            virtual volScalarField& pRef();

            //- Reference Pressure [K]
            virtual const volScalarField& pRef() const;
            
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
};
#endif

