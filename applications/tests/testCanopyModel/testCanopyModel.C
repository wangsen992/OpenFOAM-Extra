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
#include "dynamicMomentumTransportModel.H"
#include "fvMesh.H"
#include "typeInfo.H"

#include "radiationModel.H"
#include "fvDOM.H"
#include "kEpsilon.H"

#include "cyclicFvPatchField.H"

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


    radiationModels::fvDOM dom(tThermo->T());
    Info << dom.dictName() << endl;
    Info << dom.subDict("fvDOMCoeffs") << endl;
    vector d0 = normalised(dom.subDict("fvDOMCoeffs").lookup<vector>("primaryDir"));
    label bandI0 = dom.subDict("fvDOMCoeffs").lookup<label>("primaryBand");
    label rayI0 = dom.subDict("fvDOMCoeffs").lookup<label>("primaryRay");
    
    Info << d0 << endl;
    Info << bandI0 << endl;
    Info << rayI0 << endl;


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
    volScalarField qin_sw
    (
      IOobject
      (
        "qin_sw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dom.qr().dimensions()
    );
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
      dom.Ru()->dimensions()
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

    volScalarField& aTree = mesh.lookupObjectRef<volScalarField>("aTree");
    volScalarField& a = const_cast<volScalarField&>(dom.a());
    PtrList<volScalarField>& aLambda = const_cast<PtrList<volScalarField>&>(dom.aLambda());

    // RunTime Starts

    // delta value for radiation direction
    label nSteps = 48;
    label nCounter = 1;

    while (runTime.loop())
    {
      Info << "Time: " << runTime.timeName() << endl;

      tThermo->correct();
      treelist.correct();

      // Replace dom.correct() with actually correct each ray
      // dom.correct();
      d0.x() = -std::cos(constant::mathematical::pi / nSteps * nCounter);
      d0.y() = 0.3;
      d0.z() = -std::sin(constant::mathematical::pi / nSteps * nCounter);
      nCounter++;

      tensor rotation = rotationTensor(dom.IRay(bandI0).d(), normalised(d0));

      radiationModels::blackBodyEmission& blackbody = const_cast<radiationModels::blackBodyEmission&>(dom.blackBody());
      for (label i = 0; i < dom.nLambda(); i++)
      {
        blackbody.correct(i, dom.absorptionEmission().bands(i));
      }

      for(label i = 0; i < dom.nRay(); i++)
      {
          radiationModels::radiativeIntensityRay& ray = const_cast<radiationModels::radiativeIntensityRay&>(dom.IRay(i));

          // rotate rays to prescribed angle
          vector& rayd = const_cast<vector&>(ray.d());
          vector& raydAve = const_cast<vector&>(ray.dAve());
          rayd = rotation & rayd;
          raydAve = rotation & raydAve;
          Info << "Ray " << i << " d(): " << ray.d() << endl;
          Info << "Ray " << i << " dAve(): " << ray.dAve() << endl;

          // When the two z's are opposite sign, divergence happens
          // This is prone to happen when z is close to zero
          if (rayd.z() * raydAve.z() < 0)
          {
              raydAve.z() = (-1) * raydAve.z();
          }
          // Given periodic horizontal boundary condition, z close to zero will
          // cause divergence
          if (mag(raydAve.z()) < pow(0.1,3))
          {
              Info << "raydAve.z() = " << raydAve.z() << endl;
              raydAve.z() = (-1) * 0.1 / dom.nTheta();
          }

          // Update absorption coefficient
          // [Note] Disabled this part to speed up the test
          // dimensionedScalarCellSet LaCov = treelist[0].canopy().correctLaCov(rayd);
          // forAllConstIter(labelHashSet, treelist[0].canopy().canopyCells(), iter)
          // {
          //   aTree[iter.key()] = LaCov[iter.key()].value();
          // }

          // // Constraining upward shortwave radiation at boundary patch floor
          volScalarField& rayI = const_cast<volScalarField&>(ray.ILambda(0));
          scalarField& Ibfi = rayI.boundaryFieldRef()[3];
          forAll(rayI.boundaryFieldRef(), patchi)
          {
              if (rayI.boundaryFieldRef()[patchi].patch().name() == "floor" && rayd.z() > 0)
              {
                Info << "Constraining floor side upward shortwave radiation to zero" << endl;
                Ibfi = rayI.boundaryFieldRef()[patchi];
                Ibfi = 0.0;
              }
          }
              
          dom.absorptionEmission().correct(a, aLambda);
          ray.correct();
          
          // // Do another round to save results
          // // Constraining upward shortwave radiation at boundary patch floor
          // forAll(rayI.boundaryFieldRef(), patchi)
          // {
          //     if (rayI.boundaryFieldRef()[patchi].patch().name() == "floor" && rayd.z() > 0)
          //     {
          //       Info << "Constraining floor side upward shortwave radiation to zero" << endl;
          //       Ibfi = rayI.boundaryFieldRef()[patchi];
          //       Ibfi = 0.0;
          //     }
          // }
          Info << "Floor shortwave I value mean = " << average(Ibfi) << endl;
      }
      dom.updateG();
      // Update direction

      Ru.field() = dom.Ru();
      Info << "assigning values to volFields.." << endl;
      forAllConstIter(labelHashSet, canopy.canopyCells(), iter)
      {
          Fhe[iter.key()] = FheCanopy[iter.key()].value();
          Fle[iter.key()] = 1.2 * pow(10,6) * FqCanopy[iter.key()].value();
      }

      runTime.writeNow();
      if (nCounter == nSteps){break;}

    }

    return 0;
}
