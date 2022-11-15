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

#include "IOobject.H"
#include "boussinesqBuoyancyKSource.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(boussinesqBuoyancyKSource, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        boussinesqBuoyancyKSource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::boussinesqBuoyancyKSource::readCoeffs()
{
    phaseName_ = coeffs().lookupOrDefault<word>("phase", word::null);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::boussinesqBuoyancyKSource::boussinesqBuoyancyKSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    phaseName_(word::null),
    theta0_
    (
      mesh.lookupObjectRef<uniformDimensionedScalarField>("theta0")
    ),
    g_
    (
      mesh.lookupObjectRef<uniformDimensionedVectorField>("g")
    ),
    thermo_
    (
      mesh.lookupObjectRef<fluidAtmThermo>("thermophysicalProperties")
    ),
    alphat_
    (
      mesh.lookupObjectRef<volScalarField>("alphat")
    )
{
    readCoeffs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::boussinesqBuoyancyKSource::addSupFields() const
{
    return wordList(1, "k");
}


void Foam::fv::boussinesqBuoyancyKSource::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    tmp<volVectorField> tgradTheta(fvc::grad(thermo_.theta_v()));
    const volVectorField& gradTheta(tgradTheta.ref());
    // Warning: Doesn't work in incompressible mode since alphat has rho in it
    eqn += alphat_ * g_ & gradTheta / theta0_;
}


void Foam::fv::boussinesqBuoyancyKSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    Info << "[boussinesqBuoyancyKSource.C] runing addSup with rho and theta" << endl;
    tmp<volVectorField> tgradTheta(fvc::grad(thermo_.theta()));
    const volVectorField& gradTheta(tgradTheta.ref());

    eqn += alphat_ * g_ & gradTheta / theta0_;
}


void Foam::fv::boussinesqBuoyancyKSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    tmp<volVectorField> tgradTheta(fvc::grad(thermo_.theta_v()));
    const volVectorField& gradTheta(tgradTheta.ref());
    eqn += alpha * alphat_ * g_ & gradTheta / theta0_;
}


bool Foam::fv::boussinesqBuoyancyKSource::read(const dictionary& dict)
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
