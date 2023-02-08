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
Class
    Foam::canopyMomentumTransferModel

Description
    Base class to provide momentumTransport (include velocity, and associated turbulence
    parameters. 
    
SourceFiles
    canopyMomentumTransferModel.C
\*---------------------------------------------------------------------------*/

#ifndef canopyEnergyTransferModel_H
#define canopyEnergyTransferModel_H

#include "treeModel.H"
#include "fluidAtmThermo.H"
#include "cellSetHashTables.H"

namespace Foam
{
namespace urban
{

template<class BaseCanopyModel>
class canopyEnergyTransferModel
{
    
    dimensionedScalarCellSet  Tleaf;
    dimensionedScalarCellSet  Cpleaf;
    

public:
      
    canopyEnergyTransferModel
    (
      const treeModel& tree
    );

    virtual ~canopyEnergyTransferModel();

    void addCanopySupEnergy(fvScalarMatrix&);

};

}
}

#endif

