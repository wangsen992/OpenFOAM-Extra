/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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

#include "fluidAtmThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidAtmThermo, 0);
    defineRunTimeSelectionTable(fluidAtmThermo, fvMesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidAtmThermo::implementation::implementation
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    theta_(lookupOrConstruct(mesh, "theta")),
    lwc_(lookupOrConstruct(mesh, "lwc")),
    p0_("p0", dimPressure, pow(10,5))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidAtmThermo> Foam::fluidAtmThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<fluidAtmThermo>(mesh, phaseName);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidAtmThermo::~fluidAtmThermo()
{}


Foam::fluidAtmThermo::implementation::~implementation()
{}

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::fluidAtmThermo::exner
(
    const scalar p,
    const scalar p0,
    const scalar gamma
)
{
    return Foam::pow((p/p0), gamma);
}

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::exner
(
    const volScalarField& p, 
    const dimensionedScalar& p0,
    const volScalarField& gamma
)
{
    return Foam::pow((p/p0), gamma);
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::volScalarField& Foam::fluidAtmThermo::implementation::theta()
{
    return theta_;
}


const Foam::volScalarField& Foam::fluidAtmThermo::implementation::theta() const
{
    return theta_;
}

Foam::volScalarField& Foam::fluidAtmThermo::implementation::lwc()
{
    return lwc_;
}


const Foam::volScalarField& Foam::fluidAtmThermo::implementation::lwc() const
{
    return lwc_;
}

const Foam::dimensionedScalar Foam::fluidAtmThermo::implementation::p0() const
{
    return p0_;
}

