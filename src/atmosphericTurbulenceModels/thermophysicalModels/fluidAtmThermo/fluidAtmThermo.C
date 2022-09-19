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

Foam::volScalarField& Foam::fluidAtmThermo::implementation::q()
{
    return const_cast<volScalarField&>(this->composition().Y("H2O"));
}

const Foam::volScalarField& Foam::fluidAtmThermo::implementation::q() const
{
    return this->composition().Y("H2O");
}

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::implementation::theta_v() const
{
    
    Foam::tmp<Foam::volScalarField> ttheta_v
    (
        Foam::volScalarField::New
        (
            "theta_v",
            this->theta_ 
            * (
                1 + 0.608 * this->q()
              //  - this->lwc_ // This is currently disabled, need to be
              //  careful when implementing
              )
        )
    );
    return ttheta_v;
}

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::implementation::r() const
{
    Foam::tmp<Foam::volScalarField> tr
    (
        Foam::volScalarField::New
        (
            "r",
            this->composition().Y("H2O") / this->composition().Y("dryAir")
        )
    );
    return tr;
}

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::implementation::rl() const
{
    Info << "Returning theta for testing rl()" << endl;
    return this->theta_;
}


Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::implementation::pp
(
  const label speciei
) const
{
    return 
    ( 
      this->composition().Y(speciei) * this->composition().Wi(speciei) // R 
    * this->p()
    );
}

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::implementation::pp
(
  const word& specieName
) const
{
  const speciesTable& speciesTbl(this->composition().species());
  label speciei(speciesTbl[specieName]);
  return this->pp(speciei);
}

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::implementation::es () const
{
    Foam::tmp<Foam::volScalarField> Tcelsius
    (
        this->T() - dimensionedScalar(dimTemperature, 273.15)
    );

    Foam::tmp<Foam::volScalarField> tes
    (
        6.112 * 
        exp( 
             (17.67 * Tcelsius) 
           / (Tcelsius + dimensionedScalar(dimTemperature, 243.5)                    )
           )
    );
    return tes * dimensionedScalar(dimPressure, 1);
}
