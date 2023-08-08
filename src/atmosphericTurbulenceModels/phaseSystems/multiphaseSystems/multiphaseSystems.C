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

#include "phaseSystem.H"
#include "MomentumTransferPhaseSystem.H"
#include "OneResistanceHeatTransferPhaseSystem.H"
#include "TwoResistanceHeatTransferPhaseSystem.H"
#include "PhaseTransferPhaseSystem.H"
// #include "InterfaceCompositionPhaseChangePhaseSystem.H"
#include "atmInterfaceCompositionPhaseChangePhaseSystem.H"
#include "PopulationBalancePhaseSystem.H"
#include "ThermalPhaseChangePhaseSystem.H"
// #include "DropletNucleationPhaseChangePhaseSystem.H"
#include "atmThermalPhaseChangePhaseSystem.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    typedef
        atmThermalPhaseChangePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                TwoResistanceHeatTransferPhaseSystem
                <
                    MomentumTransferPhaseSystem<phaseSystem>
                >
            >
        >
        atmThermalPhaseChangeMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        atmThermalPhaseChangeMultiphaseSystem,
        dictionary,
        atmThermalPhaseChangeMultiphaseSystem
    );

    typedef
        atmInterfaceCompositionPhaseChangePhaseSystem
        <
          atmThermalPhaseChangePhaseSystem
          <
              PhaseTransferPhaseSystem
              <
                  TwoResistanceHeatTransferPhaseSystem
                  <
                      MomentumTransferPhaseSystem<phaseSystem>
                  >
              >
          >
        >
        interfaceCompositionAtmThermalPhaseChangeMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        interfaceCompositionAtmThermalPhaseChangeMultiphaseSystem,
        dictionary,
        interfaceCompositionAtmThermalPhaseChangeMultiphaseSystem
    );

    typedef
          atmThermalPhaseChangePhaseSystem
          <
            PopulationBalancePhaseSystem
            <
                PhaseTransferPhaseSystem
                <
                    TwoResistanceHeatTransferPhaseSystem
                    <
                        MomentumTransferPhaseSystem<phaseSystem>
                    >
              >
            >
          >
        atmThermalPhaseChangePopulationBalanceMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        atmThermalPhaseChangePopulationBalanceMultiphaseSystem,
        dictionary,
        atmThermalPhaseChangePopulationBalanceMultiphaseSystem
    );

    typedef
        atmInterfaceCompositionPhaseChangePhaseSystem
        <
          atmThermalPhaseChangePhaseSystem
          <
            PopulationBalancePhaseSystem
            <
                PhaseTransferPhaseSystem
                <
                    TwoResistanceHeatTransferPhaseSystem
                    <
                        MomentumTransferPhaseSystem<phaseSystem>
                    >
              >
            >
          >
        >
        interfaceCompositionAtmThermalPhaseChangePopulationBalanceMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        interfaceCompositionAtmThermalPhaseChangePopulationBalanceMultiphaseSystem,
        dictionary,
        interfaceCompositionAtmThermalPhaseChangePopulationBalanceMultiphaseSystem
    );

}


// ************************************************************************* //
