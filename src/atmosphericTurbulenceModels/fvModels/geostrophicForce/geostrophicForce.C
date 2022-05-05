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

\*---------------------------------------------------------------------------*/

#include "geostrophicForce.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(geostrophicForce, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        geostrophicForce,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::geostrophicForce::readCoeffs()
{
    phaseName_ = coeffs().lookupOrDefault<word>("phase", word::null);

    UName_ =
        coeffs().lookupOrDefault<word>
        (
            "U",
            IOobject::groupName("U", phaseName_)
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::geostrophicForce::geostrophicForce
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    phaseName_(word::null),
    UName_(word::null),
    U_
    (
      this->mesh().lookupObjectRef<volVectorField>(UName_)
    ),
    Ug_
    (
      "Ug",
      dimVelocity,
      vector(coeffs().lookup<vector>("Ug"))
    ),
    f_
    (
      "f",
      dimTime / dimTime / dimTime,
      vector(coeffs().lookup<vector>("f"))
    )
{
    readCoeffs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::geostrophicForce::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::geostrophicForce::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    eqn += f_ ^ (U_ - Ug_);
}


void Foam::fv::geostrophicForce::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    eqn += rho*g_;
}


void Foam::fv::geostrophicForce::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    eqn += alpha*rho*g_;
}


bool Foam::fv::geostrophicForce::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
