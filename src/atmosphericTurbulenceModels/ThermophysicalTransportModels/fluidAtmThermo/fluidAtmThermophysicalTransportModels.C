/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

#include "fluidAtmThermophysicalTransportModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeThermophysicalTransportModels
(
    ThermophysicalTransportModel,
    dynamicTransportModelCompressibleMomentumTransportModel,
    fluidAtmThermo
);


// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //

#include "Fourier.H"
makeLaminarThermophysicalTransportModel(Fourier);

#include "unityLewisFourier.H"
makeLaminarThermophysicalTransportModel(unityLewisFourier);

#include "FickianFourier.H"
makeLaminarThermophysicalTransportModel(FickianFourier);

#include "MaxwellStefanFourier.H"
makeLaminarThermophysicalTransportModel(MaxwellStefanFourier);


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "eddyDiffusivity.H"
makeRASLESThermophysicalTransportModel(RAS, eddyDiffusivity);

#include "unityLewisEddyDiffusivity.H"
makeRASLESThermophysicalTransportModel(RAS, unityLewisEddyDiffusivity);

#include "nonUnityLewisEddyDiffusivity.H"
makeRASLESThermophysicalTransportModel(RAS, nonUnityLewisEddyDiffusivity);
// 
// #include "FickianEddyDiffusivity.H"
// makeRASLESThermophysicalTransportModel(RAS, FickianEddyDiffusivity);


// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

#include "eddyDiffusivity.H"
makeRASLESThermophysicalTransportModel(LES, eddyDiffusivity);

#include "unityLewisEddyDiffusivity.H"
makeRASLESThermophysicalTransportModel(LES, unityLewisEddyDiffusivity);

#include "nonUnityLewisEddyDiffusivity.H"
makeRASLESThermophysicalTransportModel(LES, nonUnityLewisEddyDiffusivity);

#include "Deardorff1980EddyDiffusivity.H"
makeRASLESThermophysicalTransportModel(LES, Deardorff1980EddyDiffusivity);

// ************************************************************************* //
