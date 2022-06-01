/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "canopyPhysicsModel.H"

namespace Foam
{

canopyPhysicsModel::canopyPhysicsModel
(
    word modelName,
    canopySurfaceModel canopySurface,
    fluidAtmThermophysicalTransportModel& transport,
    radiationModel& radiation
)
:
    modelName_(modelName),
    canopySurface_(canopySurface),
    transport_(transport),
    radiation_(radiation),
    fU_(canopySurface_.canopyCells().size()),
    fT_(canopySurface_.canopyCells().size()),
    fq_(canopySurface_.canopyCells().size()),

    // Initiated with zero size
    fk_(0),
    feps_(0),
    fomega_(0),
    fR_(0)
{
}

void canopyPhysicsModel::correct()
{
    this->correctU();
    this->correctTurb();
    this->correctThermo();
}

}
