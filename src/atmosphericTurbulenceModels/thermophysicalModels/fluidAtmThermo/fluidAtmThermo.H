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

    // Static Function
    static scalar exner
    (
        const scalar p, 
        const scalar p0, 
        const scalar gamma
    );

    static tmp<volScalarField> exner
    (
        const volScalarField& p, 
        const dimensionedScalar& p0,
        const volScalarField& gamma
    );

    static tmp<volScalarField> exner
    (
        const volScalarField& p, 
        const dimensionedScalar& p0,
        const scalar gamma
    );
    
    // Member Functions

        // Access to thermodynamic state variables

            //- Potential Temperature [K]
            //  Non-const access allowed for transport equations
            virtual volScalarField& theta() = 0;

            //- Potential Temperature [K]
            virtual const volScalarField& theta() const = 0;
           
            //- Specific Liquid Water Content
            //  Non-const access allowed for transport equations
            virtual volScalarField& ql() = 0;

            //- Specific Liquid Water Content
            virtual const volScalarField& ql() const = 0;
            
            //- Reference Pressure [Pa]
            virtual const dimensionedScalar p0() const = 0;
            
            //- Specific Latent Heat of Water [J/kg]
            virtual const dimensionedScalar Lv() const = 0;
            

        // Access to Mixture Specie Information
        
            //- Specific Humidity [kg/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& q();

            //- Specific Humidity [kg/kg]
            // Const access allowed for transport equations
            virtual const volScalarField& q() const;
            
        
        // Derived variables from available references (see above) 
            //- Virtual Potential Temperature [K]
            //  Non-const access allowed for transport equations
            virtual tmp<volScalarField> theta_v() const;

            //- Liquid Water Potential Temperature [K]
            //  Non-const access allowed for transport equations
            virtual tmp<volScalarField> theta_l() const;
            
            //- Total Specific Humidity
            //  Non-const access allowed for transport equations
            virtual tmp<volScalarField> qw() const;

            //- Liquid Water Content (LWC) [kg m-3]
            // implemented in rhoAtmThermo for rho access
            virtual tmp<volScalarField> lwc() const = 0;
       
            //- Mixing Ratio []
            virtual tmp<volScalarField> r() const;

            //- Liquid Mixing Ratio []
            // Mass of liquid over Mass of Dry Air (similar to mixing ratio) 
            virtual tmp<volScalarField> rl() const;

            //- Return the partial pressure for a specie given by index
            virtual tmp<volScalarField> pp(const label specie) const;

            //- Return the partial pressure for a specie given by name
            virtual tmp<volScalarField> pp(const word& specieName) const;
            //- Return the saturation vapor pressure
            virtual tmp<volScalarField> es() const;
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
            dimensionedScalar  p0_;

            //- Specific Latent Heat of Vaporization [J/kg]
            dimensionedScalar Lv_;
        // Fields

            //- Potential Temperature [K]
            volScalarField theta_;
            
            //- Specific Liquid Water Content
            // [Note] ql is initialized to zero assuming non-saturated
            // moist air. For saturated air, liquid water content is updated
            // from application level, either by solving theta_l-ql relation 
            // or by multi-phase flow system (PBM or LCM)
            volScalarField ql_;
            

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
            
            //- Specific Liquid Water Content
            //  Non-const access allowed for transport equations
            virtual volScalarField& ql();

            //- Specific Liquid Water Content
            virtual const volScalarField& ql() const;

            //- Reference Pressure [K]
            virtual const dimensionedScalar p0() const;

            //- Specific Latent Heat of Water [J/kg]
            virtual const dimensionedScalar Lv() const;
            
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
};
#endif

