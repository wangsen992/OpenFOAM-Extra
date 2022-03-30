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

#include "atmTurbMesh.H"

// Momentum Equation Solution setup
tmp<fvVectorMatrix> Foam::atmTurbMesh::UEqn()
{
    Info << "Init UEqn() " << endl;
    // Momentum predictor
    Info << "Test Turbulence object" << endl;
    Info << turbulence_->filePath() << endl;
    Info << turbulence_->divDevTau(U_)<< endl;
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U_)
      + fvm::div(phi_, U_)
      // + turbulence_->divDevSigma(U_)
      + fU_Ug() // Geostrohpic Term
      ==
      -fvc::grad(p_) / this->rho0_ 
      + g_ * (theta_ - theta0_) / (theta0_)
    );
    Info << "UEqn Constructed. " << endl;
    return tUEqn;
}

void Foam::atmTurbMesh::solve_U()
{
    // fvVectorMatrix& UEqn(this->UEqn());
    // UEqn.relax();
    // UEqn.solve();
}

// Pressure corrector to enforce continuity
// Solution control is ignored here
tmp<fvScalarMatrix> Foam::atmTurbMesh::pEqn()
{
    fvVectorMatrix& UEqn(this->UEqn().ref());
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA
    (
      constrainHbyA(rAU*UEqn.H(), this->U_, this->p_)
    );
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    adjustPhi(phiHbyA, this->U_, this->p_);

    tmp<volScalarField> rAtU(rAU);
    // Update the pressure BCs to ensure flux consistency
    // The incompressible form of constrainPressure() used which basically
    // overloaded the compressible form with rho being one field
    constrainPressure(this->p_, this->U_, phiHbyA, rAtU());
    tmp<fvScalarMatrix> tpEqn
    (
      fvm::laplacian(rAtU(), this->p_) == fvc::div(phiHbyA)
    );

    return tpEqn; 
}
void Foam::atmTurbMesh::solve_p()
{
  Info << "solve_p() not implemented. " << endl;
  FatalError.exit();
}

// Other variables
tmp<fvScalarMatrix> Foam::atmTurbMesh::thetaEqn()
{
  tmp<fvScalarMatrix> tThetaEqn
  (
      fvm::ddt(this->theta_)
    + fvm::div(this->phi_, this->theta_)
    == 
      // Warning: nut() is used instead of alphaEff T
      fvc::laplacian(this->turbulence_->nut(), this->theta_)
  );
  return tThetaEqn.ref();
}
    
void Foam::atmTurbMesh::solve_theta()
{
  Info << "solve_theta() not implemented. " << endl;
  FatalError.exit();
}

tmp<fvScalarMatrix> Foam::atmTurbMesh::qEqn()
{
  tmp<fvScalarMatrix> tQEqn
  (
      fvm::ddt(this->q_)
    + fvm::div(this->phi_, this->q_)
    == 
      // Warning: nut() is used instead of alphaEff Q
      fvc::laplacian(this->turbulence_->nut(), this->q_)
  );
  return tQEqn;
}

void Foam::atmTurbMesh::solve_q()
{
  Info << "solve_q() not implemented. " << endl;
  FatalError.exit();
}
