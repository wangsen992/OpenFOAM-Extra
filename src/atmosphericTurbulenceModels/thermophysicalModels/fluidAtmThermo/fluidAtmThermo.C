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
    theta_v_(lookupOrConstruct(mesh, "theta_v")),
    theta_l_(lookupOrConstruct(mesh, "theta_l")),
    ql_(lookupOrConstruct(mesh, "ql")),
    qw_(lookupOrConstruct(mesh, "qw")),
    p0_("p0", dimPressure, pow(10,5)),
    Lv_("Lv", dimEnergy/dimMass, 2.26*pow(10,6))
{
    // Pure virtual functions are being called at this point (e.g. Cp()) 
    // which can be compiled but can cause unwanted behaviour. Try moving 
    // updates below to the final thermodynamics compute engine 
    // (e.g. thetaRhoAtmThermo)
    theta_v_ = theta_ * ( 1 + 0.608 * q() - ql());

    theta_l_ = theta_ - (Lv_ / Cp() * theta_ / T()) * ql_;

    Info << "Check default value of ql: " << nl
         << average(ql_) << endl;

    qw_ = q() + ql();
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

Foam::volScalarField& Foam::fluidAtmThermo::implementation::theta_v()
{
    
    return theta_v_;
}

const Foam::volScalarField& Foam::fluidAtmThermo::implementation::theta_v() const
{
    
    return theta_v_;
}

Foam::volScalarField& Foam::fluidAtmThermo::implementation::ql()
{
    
    return ql_;
}

const Foam::volScalarField& Foam::fluidAtmThermo::implementation::ql() const
{
    
    return ql_;
}

Foam::volScalarField& Foam::fluidAtmThermo::implementation::qw()
{
    
    return qw_;
}

const Foam::volScalarField& Foam::fluidAtmThermo::implementation::qw() const
{
    
    return qw_;
}

const Foam::dimensionedScalar Foam::fluidAtmThermo::implementation::p0() const
{
    return p0_;
}

const Foam::dimensionedScalar Foam::fluidAtmThermo::implementation::Lv() const
{
    return Lv_;
}

Foam::volScalarField& Foam::fluidAtmThermo::implementation::q()
{
    return const_cast<volScalarField&>(this->composition().Y("H2O"));
}

const Foam::volScalarField& Foam::fluidAtmThermo::implementation::q() const
{
    return this->composition().Y("H2O");
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
    Foam::tmp<Foam::volScalarField> trl
    (
        Foam::volScalarField::New
        (
            "rl",
            ql_ / this->composition().Y("dryAir")
        )
    );
    return trl;
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
