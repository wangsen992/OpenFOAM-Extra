/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2021 OpenFOAM Foundation
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

#include "dropletActivation.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvmDdt.H"
#include "fvcDdt.H"

#include "fluidAtmThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(dropletActivation, 0);
    addToRunTimeSelectionTable
    (
        nucleationModel,
        dropletActivation,
        dictionary
    );
}
}
}

using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::nucleationModels::dropletActivation::
dropletActivation
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    nucleationModel(popBal, dict),
    dNuc_("nucleationDiameter", dimLength, dict),
    velGroup_
    (
        refCast<const velocityGroup>
        (
            popBal.mesh().lookupObject<phaseModel>
            (
                IOobject::groupName
                (
                    "alpha",
                    dict.lookup("velocityGroup")
                )
            ).dPtr()()
        )
    ),
    nucPhase_
    (
        popBal_.mesh().lookupObject<phaseModel>
        (
            IOobject::groupName("alpha", dict.lookup("nucPhase"))
        )
    ),
    pair_
    (
        popBal_.fluid().phasePairs()
        [
            phasePair(velGroup_.phase(), nucPhase_)
        ]
    ),
    specieName_(dict.lookup("specie")),
    e_
    (
        IOobject
        (
          "dropletActivation.e",
          nucPhase_.time().timeName(),
          nucPhase_.mesh(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        nucPhase_.mesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    es_
    (
        IOobject
        (
          "dropletActivation.es",
          nucPhase_.time().timeName(),
          nucPhase_.mesh(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        nucPhase_.mesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    S_
    (
        IOobject
        (
          "dropletActivation.S",
          nucPhase_.time().timeName(),
          nucPhase_.mesh(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        nucPhase_.mesh(),
        dimensionedScalar(dimless, Zero)
    )
{
    if
    (
        dNuc_.value() < velGroup_.sizeGroups().first().dSph().value()
     || dNuc_.value() > velGroup_.sizeGroups().last().dSph().value()
    )
    {
        FatalIOErrorInFunction(dict)
            << "Nucleation diameter " << dNuc_.value() << "m outside of range ["
            << velGroup_.sizeGroups().first().dSph().value() << ", "
            << velGroup_.sizeGroups().last().dSph().value() << "]." << nl
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::diameterModels::nucleationModels::dropletActivation::addToNucleationRate
(
    volScalarField& nucleationRate,
    const label i
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const volScalarField& rho = fi.phase().rho();

    const volScalarField& q = nucPhase_.Y()[1];
    Info << "[dropletActivation] q name : " << q.name() << endl;
    const volScalarField& p = fi.phase().thermo().p();
    volScalarField Tc
    (
        fi.phase().thermo().T() 
      - dimensionedScalar(dimTemperature, 273.15)
    );

    // Prototyping model for Twomney's formula
    // It's easy to implement without modeling the transport for CCN
    forAll(nucleationRate, i)
    {
        e_[i] = q[i] * p[i] / (0.622 + 0.378 * q[i]);
        es_[i] = 611.21 * exp(Tc[i] * (18.678 - Tc[i]/234.5)/(257.14 + Tc[i]));
        S_[i] = max(0, (e_[i] - es_[i]) / es_[i]);
    }
    
    const scalar dmidtfSign =
        velGroup_.phase().name() == pair_.first() ? +1 : -1;

    nucleationRate +=
        popBal_.eta(i, pi/6*pow3(dNuc_))*dmidtfSign
        * dimensionedScalar(inv(dimVol * dimTime), 50 * pow(10, 6)) 
        * pow(S_, 0.5);
}


// ************************************************************************* //
