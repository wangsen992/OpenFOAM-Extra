
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

Description
    Utilities functions to load fields needed. Replacing createFields.H code header file which is not good C++ practice. 
\*---------------------------------------------------------------------------*/

#include "atmTurbModel.H"
#include "CorrectPhi.H"

// Pressure corrector to enforce continuity
// Solution control is ignored here
void Foam::atmTurbModel::pressureCorrect()
{
    Info << "Init pEqn." << endl;
    volScalarField rAU("rAU", 1.0/UEqn_->A());
    volVectorField HbyA
    (
      constrainHbyA(rAU*UEqn_->H(), U_, p_rgh_)
    );
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    if (p_rgh_.needReference())
    {
        fvc::makeRelative(phiHbyA, U_);
        adjustPhi(phiHbyA, U_, p_rgh_);
        fvc::makeAbsolute(phiHbyA, U_);
    }

    tmp<volScalarField> rAtU(rAU);
    // Update the pressure BCs to ensure flux consistency
    // The incompressible form of constrainPressure() used which basically
    // overloaded the compressible form with rho being one field
    if (pimple_.consistent())
    {
      rAtU = 1.0/max(1.0/rAU - UEqn_->H1(), 0.1/rAU);
      phiHbyA += 
          fvc::interpolate(rAtU() - rAU) * fvc::snGrad(p_rgh_)*mesh_.magSf();
      HbyA -= (rAU - rAtU()) * fvc::grad(p_rgh_);
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh_, U_, phiHbyA, rAtU());

    // Non-orthogonal pressure corrector loop ignored here
    while (pimple_.correctNonOrthogonal())
    {
      // RHS multiplied by rho_ to maintain dimension balance
      fvScalarMatrix pEqn
      (
        fvm::laplacian(rAtU(), p_rgh_) == fvc::div(phiHbyA) 
      );
      pEqn.setReference
      (
        pressureReference_.refCell(), 
        pressureReference_.refValue()
      );

      pEqn.solve();

      if (pimple_.finalNonOrthogonalIter())
      {
        Info << "Correcting phi_ in pEqn()." << endl;
        phi_ = phiHbyA - pEqn.flux();
        rhophi_ = rho0f_ * phi_;
        Info << "phi correcting complete." << endl;
      }
    }

    // Explicitly relax pressure for UEqn
    Info << "timeIndex : " << mesh_.time().timeIndex() << endl; 
    if (mesh_.time().timeIndex() > 1)
    {
      Info << "Relaxing p field for U correction." << endl;
      p_rgh_.relax();
    }
    U_ = HbyA - rAtU * fvc::grad(p_rgh_);
    U_.correctBoundaryConditions();
}

void Foam::atmTurbModel::nutCorrect()
{
    thermo_->correct();
    turbulence_->correct();
    transport_->correct();
}

void Foam::atmTurbModel::phiCorrect()
{
    phi_ = fvc::flux(U_);
    Info << "CorrectPhi" << endl;
    CorrectPhi
    (
      phi_, 
      U_, 
      thermo_->p(),
      dimensionedScalar("rAUf", dimTime, 1),
      geometricZeroField(),
      pressureReference_,
      pimple_
    );
    Info << "CorrectPhi completed" << endl;
}
