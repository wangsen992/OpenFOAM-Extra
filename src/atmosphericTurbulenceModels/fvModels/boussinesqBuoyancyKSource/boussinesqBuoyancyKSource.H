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
    Foam::fv::boussinesqBuoyancyKSource

Description
    Calculates and applies the buoyancy contribution to TKE based on Boussinesq 
    Approximation to momentum equation

Usage
    Example usage:
    \verbatim
    boussinesqBuoyancyKSource
    {
        type        boussinesqBuoyancyKSource;
    }
    \endverbatim

SourceFiles
    boussinesqBuoyancyKSource.C

\*---------------------------------------------------------------------------*/

#ifndef boussinesqBuoyancyKSource_H
#define boussinesqBuoyancyKSource_H

#include "fvModel.H"
#include "uniformDimensionedFields.H"
#include "dimensionedVector.H"
#include "fluidAtmThermo.H"
#include "fluidAtmThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class buoyancyForce Declaration
\*---------------------------------------------------------------------------*/
 
class boussinesqBuoyancyKSource
:
    public fvModel
{
    // Private Data

        //- Optional phase name
        word phaseName_;

        //- Reference Temperature
        const uniformDimensionedScalarField& theta0_;

        //- Gravitational Acceleration
        const uniformDimensionedVectorField& g_;

        //- thermo model
        const fluidAtmThermo& thermo_;
        
        //- Turbulent diffusivity of heat K_H
        const volScalarField& alphat_;



    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();
        

public:

    //- Runtime type information
    TypeName("boussinesqBuoyancyKSource");

    // Constructors

        //- Construct from explicit source name and mesh
        boussinesqBuoyancyKSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        boussinesqBuoyancyKSource(const boussinesqBuoyancyKSource&) = delete;

    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluate

            //- Add explicit contribution to incompressible momentum equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const boussinesqBuoyancyKSource&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
