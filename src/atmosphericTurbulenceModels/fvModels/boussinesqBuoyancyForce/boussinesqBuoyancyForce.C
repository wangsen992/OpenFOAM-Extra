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

#include "boussinesqBuoyancyForce.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(boussinesqBuoyancyForce, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        boussinesqBuoyancyForce,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::boussinesqBuoyancyForce::readCoeffs()
{
    phaseName_ = coeffs().lookupOrDefault<word>("phase", word::null);

    UName_ =
        coeffs().lookupOrDefault<word>
        (
            "UName",
            IOobject::groupName("U", phaseName_)
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::boussinesqBuoyancyForce::boussinesqBuoyancyForce
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    phaseName_(word::null),
    UName_
    (
        coeffs().lookupOrDefault<word>
        (
            "UName",
            IOobject::groupName("U", phaseName_)
        )
    ),
    U_
    (
      mesh.lookupObjectRef<volVectorField>(UName_)
    ),
    theta0_
    (
        "theta0", 
        dimTemperature, 
        scalar(coeffs().lookup<scalar>("theta0"))
    ),
    g_
    (
        mesh.lookupObjectRef<uniformDimensionedVectorField>("g")
    ),
    pthermo_
    (
      fluidAtmThermo::New(mesh)
    )
{
    readCoeffs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::boussinesqBuoyancyForce::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::boussinesqBuoyancyForce::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    eqn -= g_ * (pthermo_->theta_v() - theta0_) / theta0_;
}


void Foam::fv::boussinesqBuoyancyForce::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    eqn -= rho * g_ * (pthermo_->theta_v() - theta0_) / theta0_;
}


void Foam::fv::boussinesqBuoyancyForce::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    eqn -= alpha*rho * g_ * (pthermo_->theta_v() - theta0_) / theta0_;
}


bool Foam::fv::boussinesqBuoyancyForce::read(const dictionary& dict)
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
