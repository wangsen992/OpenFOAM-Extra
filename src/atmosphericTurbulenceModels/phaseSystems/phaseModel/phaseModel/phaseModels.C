/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "rhoThermo.H"
#include "rhoReactionThermo.H"
#include "rhoAtmThermo.H"

#include "combustionModel.H"

#include "phaseModel.H"
#include "ThermoPhaseModel.H"
#include "IsothermalPhaseModel.H"
#include "AnisothermalPhaseModel.H"
#include "PurePhaseModel.H"
#include "MultiComponentPhaseModel.H"
#include "InertPhaseModel.H"
#include "ReactingPhaseModel.H"
#include "AtmMovingPhaseModel.H"
#include "StationaryPhaseModel.H"
#include "PassiveMovingPhaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Native rhoReactionThermo is used for passive phase
    typedef
        AnisothermalPhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    AtmMovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoAtmThermo>
                    >
                >
            >
        >
        pureAtmPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        pureAtmPhaseModel,
        phaseSystem,
        pureAtmPhaseModel
    );

    typedef
        AnisothermalPhaseModel
        <
            MultiComponentPhaseModel
            <
                InertPhaseModel
                <
                    AtmMovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoAtmThermo>
                    >
                >
            >
        >
        multiComponentAtmPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        multiComponentAtmPhaseModel,
        phaseSystem,
        multiComponentAtmPhaseModel
    );

    typedef
        IsothermalPhaseModel
        <
            MultiComponentPhaseModel
            <
                InertPhaseModel
                <
                    AtmMovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoAtmThermo>
                    >
                >
            >
        >
        multiComponentIsothermalAtmPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        multiComponentIsothermalAtmPhaseModel,
        phaseSystem,
        multiComponentIsothermalAtmPhaseModel
    );

    typedef
        AnisothermalPhaseModel
        <
            MultiComponentPhaseModel
            <
                ReactingPhaseModel
                <
                    AtmMovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoAtmThermo>
                    >
                >
            >
        >
        reactingAtmPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        reactingAtmPhaseModel,
        phaseSystem,
        reactingAtmPhaseModel
    );
}

// ************************************************************************* //
