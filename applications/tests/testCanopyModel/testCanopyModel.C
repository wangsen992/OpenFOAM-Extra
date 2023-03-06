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
    This test is a simple model that does not solve dynamics, but only radiation
    and thermo correction

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
      dimensionedScalar(dom.qr().dimensions(), 0)
    );

    volScalarField& aTree = mesh.lookupObjectRef<volScalarField>("aTree");
    volScalarField& a = const_cast<volScalarField&>(dom.a());
    PtrList<volScalarField>& aLambda = const_cast<PtrList<volScalarField>&>(dom.aLambda());

    // RunTime Starts

    // dom properties
    radiationModels::blackBodyEmission& blackbody = const_cast<radiationModels::blackBodyEmission&>(dom.blackBody());
    radiationModels::absorptionEmissionModel& absorptionEmission = const_cast<radiationModels::absorptionEmissionModel&>(dom.absorptionEmission());

    while (runTime.loop())
    {
      Info << "Time: " << runTime.timeName() << endl;

      tThermo->correct();
      // treelist.correct();

      for (label i = 0; i < dom.nLambda(); i++)
      {
        blackbody.correct(i, dom.absorptionEmission().bands(i));
      }

      // Rotate the solid angles
      for(label i = 0; i < dom.nRay(); i++)
      {
          radiationModels::radiativeIntensityRay& ray = const_cast<radiationModels::radiativeIntensityRay&>(dom.IRay(i));

          dom.absorptionEmission().correct(a, aLambda);
          ray.correct();
      }
      dom.updateG();
      // Update direction
      Info << "Floor average_qin : " 
           << runTime.timeName() << "," 
           << dom.qin().boundaryField()[0].patch().name() << ", "
           << average(dom.qin().boundaryField()[0]) << endl;
      Info << "Floor average_qem : " 
           << runTime.timeName() << "," 
           << dom.qin().boundaryField()[0].patch().name() << ", "
           << average(dom.qem().boundaryField()[0]) << endl;


      runTime.writeNow();
    }

    return 0;
}
