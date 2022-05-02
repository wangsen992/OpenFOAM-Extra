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

#include "fluidReactionThermo.H"

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
            
            //- Moisture (specific humidity) [kg m-3]
            //  Non-const access allowed for transport equations
            virtual volScalarField& q() = 0;

            //- Moisture (specific humidity) [kg m-3]
            virtual const volScalarField& q() const = 0;
            
            //- Liquid Water Content (LWC) [kg m-3]
            //  Non-const access allowed for transport equations
            virtual volScalarField& lwc() = 0;

            //- Liquid Water Content (LWC) [kg m-3]
            virtual const volScalarField& lwc() const = 0;
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
            volScalarField& q_;

            //- Liquid Water Content (LWC) [kg m-3]
            volScalarField& lwc_;
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
            
            //- Moisture (specific humidity) [kg m-3]
            //  Non-const access allowed for transport equations
            virtual volScalarField& q();

            //- Moisture (specific humidity) [kg m-3]
            virtual const volScalarField& q() const;
            
            //- Liquid Water Content (LWC) [kg m-3]
            //  Non-const access allowed for transport equations
            virtual volScalarField& lwc();

            //- Liquid Water Content (LWC) [kg m-3]
            virtual const volScalarField& lwc() const;
};


}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
