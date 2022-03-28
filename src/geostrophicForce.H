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
    Foam::geostrophicForce

Description

SourceFiles
    geostrophicForceI.H
    geostrophicForce.C
    geostrophicForceIO.C

\*---------------------------------------------------------------------------*/

#ifndef geostrophicForce_H
#define geostrophicForce_H

#include "fvModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class geostrophicForce;
Istream& operator>>(Istream&, geostrophicForce&);
Ostream& operator<<(Ostream&, const geostrophicForce&);


/*---------------------------------------------------------------------------*\
                          Class geostrophicForce Declaration
\*---------------------------------------------------------------------------*/

class geostrophicForce
:
    public fvModel
{
    // Private Data

        //- Description of data_

        word            UName_;

    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    // Static Data Members

        //- Static data staticData

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

        // Access

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

        // Write


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const geostrophicForce&) = delete;

    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, geostrophicForce&);
        friend Ostream& operator<<(Ostream&, const geostrophicForce&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
