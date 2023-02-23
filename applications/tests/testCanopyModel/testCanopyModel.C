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
#include "treeModelList.H"
#include "treeModel.H"
#include "canopyModel.H"
#include "fluidAtmThermo.H"
#include "radiationModel.H"
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


    autoPtr<radiationModel> tradiation
    (
      radiationModel::New(tThermo->T())
    );


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


    treeModelList treelist(mesh, urbanDict.subDict("treeList"));

    canopyModel& canopy(treelist[0].canopy());
    dimensionedScalarCellSet& FheCanopy(canopy.Fhe());
    dimensionedScalarCellSet& FqCanopy(canopy.Fq());
    volScalarField Ru
    (
      IOobject
      (
        "Ru",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      tradiation->Ru()->dimensions()
    );
    volScalarField Fhe
    (
      IOobject
      (
        "Fhe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dimEnergy/dimTime, 0)
    );
    volScalarField Fle
    (
      IOobject
      (
        "Fle",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dimEnergy/dimTime, 0)
    );

    while (runTime.loop())
    {
      Info << "Time: " << runTime.timeName() << endl;

      tThermo->correct();
      treelist.correct();
      tradiation->correct();


      Ru.field() = tradiation->Ru();
      Info << "assigning values to volFields.." << endl;
      forAllConstIter(labelHashSet, canopy.canopyCells(), iter)
      {
          Fhe[iter.key()] = FheCanopy[iter.key()].value();
          Fle[iter.key()] = 1.2 * pow(10,6) * FqCanopy[iter.key()].value();
      }

      runTime.writeNow();
    }

    return 0;
}
