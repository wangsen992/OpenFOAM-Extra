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

Application
    testMomentumModel

Description
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "utils.H"
#include "fluidAtmThermo.H"
#include "dynamicMomentumTransportModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList args(argc, argv);
    Time runTime(Time::controlDictName, args);
    fvMesh mesh
    (
      IOobject
      (
        fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
      )
    );

    Info << fvMesh::defaultRegion << endl;

    volVectorField& U(utils::lookupOrConstructVector(mesh, "U"));
    autoPtr<fluidAtmThermo> pthermo (fluidAtmThermo::New(mesh));

    surfaceScalarField rhophi
    (
      linearInterpolate(pthermo->rho()) * fvc::flux(U)
    );
    autoPtr<compressible::momentumTransportModel> pturb
    (
        compressible::momentumTransportModel::New
        (
          pthermo->rho(),
          U,
          rhophi,
          pthermo()
        )
    );

    return 0;
}
