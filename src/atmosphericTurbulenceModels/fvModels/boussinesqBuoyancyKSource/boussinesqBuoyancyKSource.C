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
    phaseName_(coeffs().lookupOrDefault<word>("phase", word::null)),
    thermo_
    (
      mesh.lookupObjectRef<fluidAtmThermo>
      (
        IOobject::groupName
        (
          "thermophysicalProperties",
          phaseName_
        )
      )
    ),
    alphat_
    (
      mesh.lookupObjectRef<volScalarField>
      (
        IOobject::groupName
        (
          "alphat",
          phaseName_
        )
      )
    )
{
    readCoeffs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::boussinesqBuoyancyKSource::addSupFields() const
{
    return wordList(1, IOobject::groupName("k", phaseName_));
}


void Foam::fv::boussinesqBuoyancyKSource::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    tmp<volScalarField> tgradb(fvc::grad(thermo_.b())->component(8));
    const volScalarField& gradb(tgradb.ref());
    // Warning: Doesn't work in incompressible mode since alphat has rho in it
    eqn += - alphat_ * gradb;
}


void Foam::fv::boussinesqBuoyancyKSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    eqn += fvc::reconstruct
           (
              ( 
                fvc::snGrad(thermo_.b().component(2))
              ) * thermo_.b().mesh().magSf()
            )->component(2) *  (-alphat_);
}

void Foam::fv::boussinesqBuoyancyKSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    // [Note] Not sure exactly why, but the commented block will cause
    // oscillation which travels downwards. From local 2D stationary flow test, 
    // adding KSource directly to eqn.source() seems to solve this problem.

    // This problem reappeared
    // Also, need to define the buoyancy flux source properly now, as 
    // there is a discrepancy of g2/c2 at the setting of neutral
    tmp<volScalarField> tKSource = 
        (
           fvc::reconstruct
             (
                ( 
                  fvc::snGrad(thermo_.b().component(2))
                ) * thermo_.b().mesh().magSf()
             )->component(2) 
        )

           * (- alpha * alphat_)
           ;
    const volScalarField& KSource = tKSource();

    forAll(eqn.source(), i)
    {
        eqn.source()[i] += KSource[i];
    }
    // eqn -= fvc::reconstruct
    //        (
    //           ( 
    //             fvc::snGrad(thermo_.b().component(2))
    //           ) * thermo_.b().mesh().magSf()
    //        )->component(2) * (- alpha * alphat_);
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
