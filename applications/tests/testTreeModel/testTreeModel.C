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
    testThermo

Description

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "fvCFD.H"
#include "compressibleMomentumTransportModel.H"
#include "fluidAtmThermo.H"
#include "fluidAtmThermophysicalTransportModel.H"

#include "canopySurfaceModel.H"
#include "dragCanopyPhysicsModel.H"

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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );
    autoPtr<fluidAtmThermo> tThermo
    (
        fluidAtmThermo::New(mesh)
    );
    volVectorField U
    (
      IOobject
      (
        "U", 
        runTime.timeName(), 
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
      ),
      mesh
    );
    surfaceScalarField phi
    (
      IOobject
      (
        "phi", 
        runTime.timeName(), 
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
      ),
      fvc::flux(U* tThermo->rho()) 
    );
    autoPtr<compressible::momentumTransportModel> tTurbulence
    (
      compressible::momentumTransportModel::New
      (
        tThermo->rho(), U, phi, tThermo
      )
    );
    autoPtr<fluidAtmThermophysicalTransportModel> tTransport
    (
      fluidAtmThermophysicalTransportModel::New
      (
        tTurbulence, tThermo
      )
    );

    autoPtr<radiationModel> tRad
    (
      radiationModel::New(tThermo->T())
    );

    // Initiate the tree canopy model
    triSurfaceMesh canopySurf
    (
      IOobject
      (
        "dense_tree.obj",
        runTime.constant()/"geometry",
        mesh,
        IOobject::MUST_READ
      )
    );
    canopySurfaceModel canopy1(mesh, canopySurf);
    Info << "canopySurfaceModel created." << endl;

    dragCanopyPhysicsModel dragTree(canopy1, tTransport, tRad, 0.5);
    Info << "dragCanopyPhysicsModel created." << endl;

    Info << "Models in treeCanopy: " << endl;
    Info << "type = " << dragTree.transport().momentumTransport().type() << endl;
    Info << "coeffDict = " << dragTree.transport().momentumTransport().coeffDict() << endl;
    
    Info << "Models in global: " << endl;
    Info << "type = " << tTransport->momentumTransport().type() << endl;
    Info << "coeffDict = " << tTransport->momentumTransport().coeffDict() << endl;

    runTime++;
    runTime.write();
    Info << "end" << endl;
    
    return 0;
}
