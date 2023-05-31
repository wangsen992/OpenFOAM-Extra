/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "hydrostaticInitialisation.H"
#include "atmHydrostaticInitialisation.H"

#include "fluidAtmThermo.H"
#include "fvmLaplacian.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "constrainPressure.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::atmHydrostaticInitialisation
(
    volScalarField& p_rgh,
    volScalarField& rho,
    const volVectorField& U,
    const volScalarField& gh,
    const surfaceScalarField& ghf,
    const uniformDimensionedScalarField& pRef,
    fluidAtmThermo& thermo,
    const dictionary& dict
)
{
    if (dict.lookupOrDefault<bool>("atmHydrostaticInitialisation", false))
    {
        Info << "Starting atmHydrostaticInitialisation.\n" << endl;
        const fvMesh& mesh = p_rgh.mesh();

        volScalarField& ph_rgh = regIOobject::store
        (
            new volScalarField
            (
                IOobject
                (
                    "ph_rgh",
                    "0",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );

        if (!mesh.time().restart())
        {
            word atmHydrostaticInitialisationMode(dict.lookupOrDefault<word>("atmHydrostaticInitialisationMode", word("dryFixedT")));

            if (atmHydrostaticInitialisationMode == "dryFixedT")
            {
              // Original value must be cached for energy setting
              volScalarField T_orig = thermo.T();
              volScalarField& p = thermo.p();
              volScalarField& he = thermo.he();
              p = ph_rgh + rho*gh + pRef;
              thermo.correct();
              rho = thermo.rho();
              label nCorr
              (
                  dict.lookupOrDefault<label>("nHydrostaticCorrectors", 5)
              );

              for (label i=0; i<nCorr; i++)
              {
                  surfaceScalarField rhof("rhof", fvc::interpolate(rho));

                  surfaceScalarField phig
                  (
                      "phig",
                      -rhof*ghf*fvc::snGrad(rho)*mesh.magSf()
                  );

                  // Update the pressure BCs to ensure flux consistency
                  constrainPressure(ph_rgh, rho, U, phig, rhof);
                  // [Temporary fix]
                  // ph_rgh = ph_rgh - max(ph_rgh);
                  // This is to fix the max value of ph_rgh to zero
                  Pout << "gMax(ph_rgh) = " << max(ph_rgh) << endl;
                  Pout << "gMax(ph_rgh) = " << gMax(ph_rgh) << endl;
                  ph_rgh.primitiveFieldRef() = ph_rgh.primitiveField() - gMax(ph_rgh);

                  fvScalarMatrix ph_rghEqn
                  (
                      fvm::laplacian(rhof, ph_rgh) == fvc::div(phig)
                  );

                  Info << "ph_rghEqn eqn symmetric? " << ph_rghEqn.symmetric() << endl;


                  SolverPerformance<scalar> sp = ph_rghEqn.solve();
                  p = ph_rgh + rho*gh + pRef;
                  he = thermo.he(p, T_orig);
                  thermo.correct();
                  rho = thermo.rho();

                  Info<< "Hydrostatic pressure variation "
                      << (max(p) - min(p)).value() << endl;
                  Info << "[hydroInit]" << "max(p)" << max(p) << "; ";
                  Info << "[hydroInit] max poinit : " << mesh.C()[findMax(p)] << endl;
              }

              Info<< "Final Hydrostatic pressure variation "
                  << (max(p) - min(p)).value() << endl;
              p_rgh = ph_rgh;
            }
            else if (atmHydrostaticInitialisationMode == "wetTheta") 
            {
              // Original value must be cached for energy setting
              volScalarField T_orig = thermo.T();
              volScalarField qv_orig = thermo.composition().Y("H2O");
              tmp<volScalarField> trv = qv_orig / (1 - qv_orig);

              // Load water volume fraction and subsequent water mixing ratio
              // This will be updated at the end
              volScalarField& alphaWater(mesh.lookupObjectRef<volScalarField>("alpha.water"));
              volScalarField& alphaAir(mesh.lookupObjectRef<volScalarField>("alpha.air"));
              // Load preset rt and wetTheta from dictionary
              dimensionedScalar rt("rt", dimless, dict.lookupOrDefault<scalar>("rt", 0.02));
              dimensionedScalar wetTheta("wetTheta", dimTemperature, dict.lookupOrDefault<scalar>("wetTheta", 320));
              
              volScalarField& T = thermo.T();
              volScalarField& p = thermo.p();
              volScalarField& qv = thermo.composition().Y("H2O");
              volScalarField& he = thermo.he();
              p = ph_rgh + rho*gh + pRef;
              thermo.correct();
              rho = thermo.rho();
              label nCorr
              (
                  dict.lookupOrDefault<label>("nHydrostaticCorrectors", 5)
              );


              for (label i=0; i<nCorr; i++)
              {
                  for (label j=0; j<nCorr; j++)
                  {
                    surfaceScalarField rhof("rhof", fvc::interpolate(rho));

                    surfaceScalarField phig
                    (
                        "phig",
                        -rhof*ghf*fvc::snGrad(rho)*mesh.magSf()
                    );

                    // Update the pressure BCs to ensure flux consistency
                    constrainPressure(ph_rgh, rho, U, phig, rhof);
                    // [Temporary fix]

                    fvScalarMatrix ph_rghEqn
                    (
                        fvm::laplacian(rhof, ph_rgh) == fvc::div(phig)
                    );

                    SolverPerformance<scalar> sp = ph_rghEqn.solve();
                    p = ph_rgh + rho*gh + pRef;
                    he = thermo.he(p, T_orig);
                    thermo.correct();
                    rho = thermo.rho();

                    Info<< "Hydrostatic pressure variation "
                        << (gMax(p) - gMin(p)) << endl;
                    Info << "[hydroInit]" << "max(p)" << max(p) << "; ";
                    Info << "[hydroInit] max poinit : " << mesh.C()[findMax(p)] << endl;
                  }
              }

              Info<< "Final Hydrostatic pressure variation "
                        << (gMax(p) - gMin(p)) << endl;
              p_rgh = ph_rgh;
                
            }
        }
        else
        {
            Info << "Restart condition of atmHydrostaticInitialisation.." << endl;
            thermo.correct();
            rho = thermo.rho();
        }
    }
    else
    {
        // Force p_rgh to be consistent with p
        // p_rgh = thermo.p() - thermo.rho()*gh - pRef;
    }
}


// ************************************************************************* //
