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
    Foam::fv::geostrohpicForce

Description
    Calculates and applies the geostrohpic force $f \cross (U - Ug)$ to the momentum equation
    corresponding to the specified velocity field.

Usage
    Example usage:
    \verbatim
    buoyancyForce1
    {
        type        geostrohpicForce;

        U           U;            // Name of the velocity field
        Ug          (10 0 0);     // Value of geostrophic velocity
        f           1e-4;         // Value of Coriolis Parameter
    }
    \endverbatim

SourceFiles
    geostrophicForce.C

\*---------------------------------------------------------------------------*/

#ifndef geostrophicForce_H
#define geostrophicForce_H

#include "fvModel.H"
#include "uniformDimensionedFields.H"
#include "dimensionedVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class buoyancyForce Declaration
\*---------------------------------------------------------------------------*/
 
class geostrophicForce
:
    public fvModel
{
    // Private Data

        //- Optional phase name
        word phaseName_;

        //- Name of the velocity field
        word UName_;

        //- Velocity Field
        volVectorField& U_;

        //- Geostrohpic Velocity
        dimensionedVector Ug_;

        //- Coriolis Parameter
        //- Note: For larger domains, f_ can be non-uniform.
        dimensionedVector f_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

public:

    //- Runtime type information
    TypeName("geostrophicForce");

    // Constructors

        //- Construct from explicit source name and mesh
        geostrophicForce
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        geostrophicForce(const geostrophicForce&) = delete;

    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluate

            //- Add explicit contribution to incompressible momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const geostrophicForce&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
