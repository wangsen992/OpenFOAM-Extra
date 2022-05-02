/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "className.H"
#include "coefficientMultiComponentMixture.H"
#include "coefficientWilkeMultiComponentMixture.H"
#include "constTransport.H"
#include "valueMultiComponentMixture.H"
#include "singleComponentMixture.H"
#include "SpecieMixture.H"

#include "rhoThermo.H"
#include "rhoAtmThermo.H"
#include "heRhoAtmThermo.H"

#include "forGases.H"
#include "forLiquids.H"
#include "forTabulated.H"
#include "makeAtmThermo.H"

#include "heRhoAtmThermo.C"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeRhoAtmThermos(Mixture, ThermoPhysics)                         \
    makeAtmThermos                                                        \
    (                                                                          \
        rhoThermo,                                                             \
        rhoAtmThermo,                                                     \
        heRhoAtmThermo,                                                           \
        Mixture,                                                               \
        ThermoPhysics                                                          \
    )

#define makeRhoAtmThermo(Mixture, ThermoPhysics)                          \
    makeAtmThermo                                                         \
    (                                                                          \
        rhoAtmThermo,                                                     \
        heRhoAtmThermo,                                                           \
        Mixture,                                                               \
        ThermoPhysics                                                          \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCoeffGases(makeRhoAtmThermos, coefficientMultiComponentMixture);

    forCoeffGases
    (
        makeRhoAtmThermos,
        coefficientWilkeMultiComponentMixture
    );

    forGases(makeRhoAtmThermo, singleComponentMixture);

    // forCoeffLiquids(makeRhoAtmThermos, coefficientMultiComponentMixture);
    // forLiquids(makeRhoAtmThermos, valueMultiComponentMixture);
    // forLiquids(makeRhoAtmThermo, singleComponentMixture);

    // forTabulated(makeRhoAtmThermos, valueMultiComponentMixture);
    // forTabulated(makeRhoAtmThermo, singleComponentMixture);
}

// ************************************************************************* //
