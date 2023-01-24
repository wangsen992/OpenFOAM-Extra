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

\*---------------------------------------------------------------------------*/

#include "rhoAtmThermo.H"
#include "word.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rhoAtmThermo, 0);
    defineRunTimeSelectionTable(rhoAtmThermo, fvMesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoAtmThermo::implementation::implementation
(
    const fvMesh& mesh,
    const word& phaseName
)
:
  rhoRef_
  (
    IOobject
    (
       phasePropertyName("rhoref", phaseName),
       mesh.time().constant(),
       mesh,
       IOobject::READ_IF_PRESENT,
       IOobject::AUTO_WRITE
    ),
    mesh,
    dimDensity
  ),
  b_
  (
    IOobject
    (
       phasePropertyName("b", phaseName),
       mesh.time().timeName(),
       mesh,
       IOobject::READ_IF_PRESENT,
       IOobject::AUTO_WRITE
    ),
    mesh,
    dimDensity * dimAcceleration
  )
{
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::rhoAtmThermo> Foam::rhoAtmThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<rhoAtmThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rhoAtmThermo::~rhoAtmThermo()
{}


Foam::rhoAtmThermo::implementation::~implementation()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::rhoAtmThermo::implementation::rhoRef()
{
    return rhoRef_;
}


const Foam::volScalarField& Foam::rhoAtmThermo::implementation::rhoRef() const
{
    return rhoRef_;
}

Foam::volVectorField& Foam::rhoAtmThermo::implementation::b()
{
    return b_;
}

const Foam::volVectorField& Foam::rhoAtmThermo::implementation::b() const
{
    return b_;
}

Foam::tmp<Foam::volVectorField> Foam::rhoAtmThermo::implementation::bByRho() const
{
    tmp<volVectorField> bByRho
    (
        new volVectorField
        (
            IOobject
            (
                phasePropertyName("bByRho", phaseName()),
                b_.mesh().time().timeName(),
                b_.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            b_.mesh(),
            dimAcceleration
        )
    );

    return bByRho;
}

// ************************************************************************* //
