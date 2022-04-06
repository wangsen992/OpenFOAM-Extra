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

// Momentum Equation Solution setup
tmp<fvVectorMatrix> Foam::atmTurbModel::UEqn()
{
      Info << "Constructing UEqn_ " << endl;
      // Momentum predictor
      tmp<volScalarField> rho_ = thermo_->rho();
      volScalarField& p_ = thermo_->p();
      volScalarField& T_ = thermo_->T();
      tmp<fvVectorMatrix> tUEqn
      (
          fvm::ddt(U_)
        + fvm::div(phi_, U_)
        + fvc::laplacian(turbulence_->nuEff(), U_)
        + fU_Ug() 
        == 
        - g_ * (T_ - T0_) / T0_
      );
      tUEqn->relax();

      if (pimple_.momentumPredictor())
      {
        solve(tUEqn() == -fvc::grad(p_) ); // Geostrohpic Term
      }
      Info << "UEqn Constructed. " << endl;
      UEqn_ = tUEqn;
      return UEqn_;
}


// Other variables
tmp<fvScalarMatrix> Foam::atmTurbModel::TEqn()
{
    Info << "Init TEqn." << endl;

    volScalarField& T_ = thermo_->T();
    tmp<fvScalarMatrix> tTEqn
    (
        fvm::ddt(T_)
      + fvm::div(phi_, T_)
      == 
        // Warning: nut() is used instead of alphaEff T
        fvc::laplacian(transport_->kappaEff()/(thermo_->Cp() * thermo_->rho()), T_)
    );
    tTEqn->relax();
    tTEqn->solve();
    return tTEqn;
}
    

tmp<fvScalarMatrix> Foam::atmTurbModel::qEqn()
{
    Info << "Init qEqn." << endl;
    tmp<fvScalarMatrix> tQEqn
    (
        fvm::ddt(q_)
      + fvm::div(phi_, q_)
      == 
        // Warning: nut() is used instead of alphaEff Q
        fvc::laplacian(transport_->alphaEff(), q_)
    );
    return tQEqn;
}

// Pressure corrector to enforce continuity
// Solution control is ignored here
void Foam::atmTurbModel::pressureCorrect()
{
    Info << "Init pEqn." << endl;
    volScalarField rAU("rAU", 1.0/UEqn_->A());
    volScalarField& p_ = thermo_->p();
    volVectorField HbyA
    (
      constrainHbyA(rAU*UEqn_->H(), U_, p_)
    );
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    if (p_.needReference())
    {
        fvc::makeRelative(phiHbyA, U_);
        adjustPhi(phiHbyA, U_, p_);
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
          fvc::interpolate(rAtU() - rAU) * fvc::snGrad(p_)*mesh_.magSf();
      HbyA -= (rAU - rAtU()) * fvc::grad(p_);
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_, U_, phiHbyA, rAtU());

    // Non-orthogonal pressure corrector loop ignored here
    while (pimple_.correctNonOrthogonal())
    {
      // RHS multiplied by rho_ to maintain dimension balance
      fvScalarMatrix pEqn
      (
        fvm::laplacian(rAtU(), p_) == fvc::div(phiHbyA) 
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
        Info << "phi correcting complete." << endl;
      }
    }

    // Explicitly relax pressure for UEqn
    Info << "timeIndex : " << mesh_.time().timeIndex() << endl; 
    if (mesh_.time().timeIndex() > 1)
    {
      Info << "Relaxing p field for U correction." << endl;
      p_.relax();
    }
    U_ = HbyA - rAtU * fvc::grad(p_);
    U_.correctBoundaryConditions();
}

void Foam::atmTurbModel::nutCorrect()
{
    thermo_->correct();
    turbulence_->correct();
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
