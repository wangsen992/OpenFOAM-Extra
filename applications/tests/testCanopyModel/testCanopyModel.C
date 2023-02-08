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
#include "treeModel.H"
#include "canopyModel.H"
#include "fluidAtmThermo.H"
#include "dynamicMomentumTransportModel.H"
#include "fvMesh.H"
#include "typeInfo.H"

#include "kEpsilon.H"

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

    #include "readGravitationalAcceleration.H"
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
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      fvc::flux(tThermo->rho() * U)
    );

    autoPtr<compressible::momentumTransportModel> tturbulence
    (
      compressible::momentumTransportModel::New
      (
        tThermo->rho(), 
        U, 
        phi, 
        tThermo()
      )
    );
    compressible::momentumTransportModel& turbulence = tturbulence();

    IOdictionary urbanDict
    (
      IOobject
      (
        "urbanModelProperties",
        mesh.time().constant(),
        mesh.time(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
      )
    );

    dictionary treeDict(urbanDict.subDict("tree"));
    autoPtr<treeModel> tTree
    (
      treeModel::New(mesh, treeDict)
    );

    Info << "Tree creation complete. " << endl;


    autoPtr<canopyModel> tCanopy(tTree->canopy());
    dimensionedVectorCellSet& lad(tCanopy->lad());

    volVectorField ladV
    (
      IOobject
      (
        "lad",
        runTime.constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedVector(dimArea/dimVolume, vector(0,0,0))
    );
    Info << "Loading Fu.." << endl;
    volVectorField Fu
    (
      IOobject
      (
        "Fu",
        runTime.constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedVector(dimVelocity/dimTime, vector(0,0,0))
    );
    Info << "Loading Fk.." << endl;
    volScalarField Fk
    (
      IOobject
      (
        "Fk",
        runTime.constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dimVelocity*dimVelocity/dimTime, 0)
    );

    Info << "Loading Feps.." << endl;
    volScalarField Feps
    (
      IOobject
      (
        "Feps",
        runTime.constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dimVelocity*dimVelocity/dimTime, 0)
    );

    Info << "assigning values to volFields.." << endl;
    forAllConstIter(labelHashSet, tCanopy->canopyCells(), iter)
    {
        ladV[iter.key()] = lad[iter.key()].value();
        Fu[iter.key()] = tCanopy->Fu()[iter.key()].value();
        Fk[iter.key()] = tCanopy->Fturb("k")[iter.key()].value();
        Feps[iter.key()] = tCanopy->Fturb("epsilon")[iter.key()].value();
    }
    
    // Test momentumTransferModel
    // tCanopy->correctMomentumTransfer();
     
    // This is used to validate the type of turbulence model
    Info << isA<RASModels::kEpsilon<compressible::momentumTransportModel>>(turbulence) << endl;

    runTime.writeNow();
    return 0;
}
