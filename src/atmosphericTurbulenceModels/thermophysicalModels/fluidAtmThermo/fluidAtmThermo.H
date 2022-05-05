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
    A fluid thermo extended for atm (supporting multi-component mixture)

\*---------------------------------------------------------------------------*/

#ifndef fluidAtmThermo_H
#define fluidAtmThermo_H

#include "fvCFD.H"
#include "fluidReactionThermo.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class fluidReactionThermo Declaration
\*---------------------------------------------------------------------------*/
class fluidAtmThermo
:
    virtual public fluidReactionThermo
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

            //- Potential Temperature [K]
            //  Non-const access allowed for transport equations
            virtual volScalarField& theta() = 0;

            //- Potential Temperature [K]
            virtual const volScalarField& theta() const = 0;
            
            //- Specific Humidity [kg m-3]
            //  Non-const access allowed for transport equations
            virtual volScalarField& q() = 0;

            //- Specific Humidity [kg m-3]
            // Const access allowed for transport equations
            virtual const volScalarField& q() const = 0;
            
            //- Liquid Water Content (LWC) [kg m-3]
            //  Non-const access allowed for transport equations
            virtual volScalarField& lwc() = 0;

            //- Liquid Water Content (LWC) [kg m-3]
            virtual const volScalarField& lwc() const = 0;

            //- Reference Pressure [Pa]
            virtual const dimensionedScalar p0() const = 0;
            
             //- Virtual Potential Temperature
             virtual tmp<volScalarField> theta_v() const = 0;
             
             //- Mixing Ratio []
             virtual tmp<volScalarField> r() const = 0;
 
             //- Liquid Mixing Ratio
             virtual tmp<volScalarField> rl() const = 0;

        // Access to Mixture Specie Information

            //- Return the table of species
            virtual const speciesTable& species() const = 0;

            //- Molecular Weight of Given Specie
            virtual scalar Wi(const label speciei) const = 0;

            //- Return the const mass-fraction field for a specie given by index
            virtual const volScalarField& Y(const label specie) const = 0;
            
            //- Return the const mass-fraction field for a specie given by name
            virtual const volScalarField& Y(const word& specieName) const = 0;
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

        // Fields

            //- Potential Temperature [K]
            volScalarField& theta_;

            //- Moisture (specific humidity) [kg m-3]
            // Replaced with Y("H2O")
            // volScalarField& q_;

            //- Liquid Water Content (LWC) [kg m-3]
            volScalarField& lwc_;

            //- Reference Pressure [Pa]
            dimensionedScalar   p0_;

public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const fvMesh&, const word& phaseName);


    //- Destructor
    virtual ~implementation();
    
    // Member Functions

        // Access to thermodynamic state variables

            //- Potential Temperature [K]
            //  Non-const access allowed for transport equations
            virtual volScalarField& theta();

            //- Potential Temperature [K]
            virtual const volScalarField& theta() const;
            
            //- Liquid Water Content (LWC) [kg m-3]
            //  Non-const access allowed for transport equations
            virtual volScalarField& lwc();

            //- Liquid Water Content (LWC) [kg m-3]
            virtual const volScalarField& lwc() const;

            //- Reference Pressure [K]
            virtual const dimensionedScalar p0() const;


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
};
#endif

