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
    atmTurbSimpleFoam

Description
    SIMPLE Algorithm is borrowed from simpleFoam implementation with key 
differences in the addition of evolution equation of potential temperature, 
evolution equation of water vapor (mixing ratio), and moist air thermodynamics. 

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "fvCFD.H"
#include "simpleControl.H"
#include "pressureReference.H"
#include "kinematicMomentumTransportModel.H"
#include "singlePhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // set root case
    argList::addBoolOption("dryRun", "if specified not solve");
    argList args(argc, argv);
    if (!args.checkRootCase())
    {
        FatalError.exit();
    }

    // create time
    Info << "Creat time" << nl << endl;
    Time runTime(Time::controlDictName, args);

    // create mesh
    Info  << "Create mesh for time = " 
          << runTime.timeName() 
          << nl << Foam::endl;

    fvMesh mesh (
      IOobject
      (
        fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
      )
    );

    // create control
    simpleControl simple(mesh);

    // create fields (prognostic variables)
    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
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

    // create flux (phi)
    Info << "Reading/calculating face flux field phi" << nl << endl;

    surfaceScalarField phi
    (
      IOobject
      (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      fvc::flux(U)
    );

    pressureReference pressureReference(p, simple.dict());

    mesh.setFluxRequired(p.name());



    Info << "Reading Atmospheric Turbulence Properties" << endl;
    IOdictionary atmTurbulenceProperties
    (
      IOobject
      (
        "atmTurbulenceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
      )
    );

    // Issue : if loading dimensionedVector(dict.lookup(name)) there is a
    // reading error. 
    dimensionedVector f("f", atmTurbulenceProperties.lookup("f"));
    dimensionedVector Ug("Ug", atmTurbulenceProperties.lookup("Ug"));
    dimensionedScalar rho_0("rho_0", atmTurbulenceProperties.lookup("rho_0"));
    volVectorField fU_Ug("fU_Ug", (f ^ (U - Ug)));
    
    // initContinuityErrs
    scalar cumulativeContErr = 0;
    
    // Init turbulence model
    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::momentumTransportModel> turbulence
    (
        incompressible::momentumTransportModel::New(U, phi, laminarTransport)
    );

    if (args.optionFound("dryRun"))
    {
      return 0;
    }
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (simple.loop(runTime))
    {
      Info << "Time = " << runTime.timeName() << endl;

      // Solve 
      // Momentum predictor
      tmp<fvVectorMatrix> tUEqn
      (
          fvm::div(phi, U)
        + turbulence->divDevSigma(U)
        + fU_Ug // addition of geostrophic term
      );
      fvVectorMatrix& UEqn = tUEqn.ref();

      UEqn.relax();

      Info << "simple.momentumPredictor() = " 
           << simple.momentumPredictor()
           << endl;
      if (simple.momentumPredictor())
      {
          solve(UEqn == -fvc::grad(p)) ;
      }
      
      // Pressure Corrector
      volScalarField rAU(1.0/UEqn.A());
      volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
      surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
      adjustPhi(phiHbyA, U, p);

      tmp<volScalarField> rAtU(rAU);

      Info  <<  "simple.consistent() = "
            << simple.consistent() 
            << endl;
      if (simple.consistent())
      {
          rAtU = 1.0/(1.0/rAU - UEqn.H1());
          phiHbyA +=
              fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
          HbyA -= (rAU - rAtU())*fvc::grad(p);
      }

      tUEqn.clear();

      // Update the pressure BCs to ensure flux consistency
      // The incompressible form of constrainPressure() used which basically
      // overloaded the compressible form with rho being one field
      
      constrainPressure(p, U, phiHbyA, rAtU());

      // Non-orthogonal pressure corrector loop
      while (simple.correctNonOrthogonal())
      {
          fvScalarMatrix pEqn
          (
              fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
          );

          pEqn.setReference
          (
              pressureReference.refCell(),
              pressureReference.refValue()
          );

          pEqn.solve();

          if (simple.finalNonOrthogonalIter())
          {
              phi = phiHbyA - pEqn.flux();
          }
      }

      // compute cumulativeContErr
      volScalarField contErr(fvc::div(phi));

      scalar sumLocalContErr = runTime.deltaTValue()*
          mag(contErr)().weightedAverage(mesh.V()).value();

      scalar globalContErr = runTime.deltaTValue()*
          contErr.weightedAverage(mesh.V()).value();
      cumulativeContErr += globalContErr;

      Info<< "time step continuity errors : sum local = " << sumLocalContErr
          << ", global = " << globalContErr
          << ", cumulative = " << cumulativeContErr
          << endl;

      // Explicitly relax pressure for momentum corrector
      p.relax();

      // Momentum corrector
      U = HbyA - rAtU()*fvc::grad(p);
      U.correctBoundaryConditions();

      // Write
      runTime.write();
    };

    Info<< "End\n" << endl;

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    

    return 0;
}


// ************************************************************************* //
