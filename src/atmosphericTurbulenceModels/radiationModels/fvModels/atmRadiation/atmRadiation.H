/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
    Foam::fv::atmRadiation

Description
    Calculates and applies the atmRadiation source to the energy equation.

Usage
    Example usage:
    \verbatim
    radiationCoeffs
    {
        libs    ("libradiationModels.so");
        type    atmRadiation;
    }
    \endverbatim

SourceFiles
    atmRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_H
#define radiation_H

#include "fvModel.H"
#include "uniformDimensionedFields.H"
#include "atmRadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class atmRadiation Declaration
\*---------------------------------------------------------------------------*/

class atmRadiation
:
    public fvModel
{
    // Private Data

        //- The atmRadiation model pointer
        word phase_;
        mutable autoPtr<atmRadiationModel> radiation_;


public:

    //- Runtime type information
    TypeName("atmRadiation");


    // Constructors

        //- Construct from explicit source name and mesh
        atmRadiation
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        atmRadiation(const atmRadiation&) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluate

            //- Add explicit contribution to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const atmRadiation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
