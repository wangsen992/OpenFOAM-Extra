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
    p0_("p0", dimPressure, pow(10,5)),
    Lv_("Lv", dimEnergy/dimMass, 2.26*pow(10,6)),
    theta_
    (
      IOobject
      (
        "theta",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimTemperature
    ),
    ql_
    (
      IOobject
      (
        "ql",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimless
    )
{
}

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
    const scalar poConst
)
{
    return Foam::pow((p/p0), poConst);
}

static tmp<scalarField> exner
(
    const scalarField p, 
    const scalar p0, 
    const scalarField poConst
)
{
    return Foam::pow((p/p0), poConst);
};

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::exner
(
    const volScalarField& p, 
    const dimensionedScalar& p0,
    const volScalarField& poConst
)
{
    return Foam::pow((p/p0), poConst);
}

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::exner
(
    const volScalarField& p, 
    const dimensionedScalar& p0,
    const scalar poConst
)
{
    return Foam::pow((p/p0), poConst);
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

Foam::volScalarField& Foam::fluidAtmThermo::implementation::ql()
{
    
    return ql_;
}

const Foam::volScalarField& Foam::fluidAtmThermo::implementation::ql() const
{
    
    return ql_;
}
const Foam::dimensionedScalar Foam::fluidAtmThermo::implementation::p0() const
{
    return p0_;
}

const Foam::dimensionedScalar Foam::fluidAtmThermo::implementation::Lv() const
{
    return Lv_;
}

Foam::volScalarField& Foam::fluidAtmThermo::q()
{
    return const_cast<volScalarField&>(this->composition().Y("H2O"));
}

const Foam::volScalarField& Foam::fluidAtmThermo::q() const
{
    return this->composition().Y("H2O");
}

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::poConst() const
{
    Foam::tmp<Foam::volScalarField> tpoConst
    (
        Foam::volScalarField::New
        (
            "poConst",
            (Cp() - Cv()) / Cp()
        )
    );
    return tpoConst;
    
}

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::theta_v() const
{
    Foam::tmp<Foam::volScalarField> ttheta_v
    (
        Foam::volScalarField::New
        (
            "theta_v",
            theta() * (1 + 0.608 * q() - ql())
        )
    );
    return ttheta_v;
}

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::theta_l() const
{
    Foam::tmp<Foam::volScalarField> ttheta_l
    (
        Foam::volScalarField::New
        (
            "theta_l",
            theta() - (Lv() / Cp() * theta() / T()) * ql()
        )
    );
    return ttheta_l;
}

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::qw() const
{
    Foam::tmp<Foam::volScalarField> tqw
    (
        Foam::volScalarField::New
        (
            "qw",
            ql() + q()        
        )
    );
    return tqw;
}

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::r() const
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

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::rl() const
{
    Foam::tmp<Foam::volScalarField> trl
    (
        Foam::volScalarField::New
        (
            "rl",
            ql() / this->composition().Y("dryAir")
        )
    );
    return trl;
}
Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::pp
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

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::pp
(
  const word& specieName
) const
{
  const speciesTable& speciesTbl(this->composition().species());
  label speciei(speciesTbl[specieName]);
  return this->pp(speciei);
}

Foam::tmp<Foam::volScalarField> Foam::fluidAtmThermo::es () const
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
