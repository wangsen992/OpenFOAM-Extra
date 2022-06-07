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
    Foam::canopyPhysicsModel

\*---------------------------------------------------------------------------*/

#include "dragCanopyPhysicsModel.H"
#include "absorptionEmissionModel.H"

namespace Foam
{

dragCanopyPhysicsModel::dragCanopyPhysicsModel
(
    canopySurfaceModel canopySurface,
    const fluidAtmThermophysicalTransportModel& transport,
    const radiationModel& radiation,
    scalar Cd
)
:
canopyPhysicsModel(canopySurface, transport, radiation),
Cd_(Cd),
lad_(this->canopySurface().lad()),
cells_(this->canopySurface().canopyCells().sortedToc())
{
};

void dragCanopyPhysicsModel::correctU()
{
    // Collect references for use in correct
    const vectorField& U(transport().momentumTransport().U().primitiveField());

    // tmp variables for iteration
    vector Ui;
    scalar magUi;
    scalar ladi;
    
    forAll(cells_, i)
    {
        Ui = U[cells_[i]];
        magUi = mag(Ui);
        ladi = lad_[cells_[i]].value();

        fU()[cells_[i]] = Cd_ * ladi * cmptMultiply(cmptSqr(Ui), Ui);
    };
}
void dragCanopyPhysicsModel::correctTurb()
{
}

void dragCanopyPhysicsModel::correctThermo()
{
    // correct heat capacity of tree region
    const radiationModel& radiation(this->radiation());
    const radiationModels::absorptionEmissionModel& absEms(radiation.absorptionEmission());
    const fvMesh& mesh(this->transport().momentumTransport().mesh());
    volScalarField& Cv = mesh.lookupObjectRef<volScalarField>("Cv");
    Info << "Testing access to Cv: " << endl;
    Info << Cv.name() << ": " << Cv.size() << endl;

    // correct thermal diffusivity of tree region
    // correct the radiative flux propreties 
}

void dragCanopyPhysicsModel::correct()
{
    correctU();
    correctTurb();
    correctThermo();
};

}
