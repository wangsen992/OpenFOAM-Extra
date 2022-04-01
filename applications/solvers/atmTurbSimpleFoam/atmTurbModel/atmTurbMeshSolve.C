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
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U_)
      + fvm::div(phi_, U_)
      // + turbulence_->divDevSigma(U_)
      + fU_Ug() // Geostrohpic Term
      ==
      -fvc::grad(p_) / rho_
    );
    Info << "UEqn Constructed. " << endl;
    tUEqn->relax();
    tUEqn->solve();
    return tUEqn;
}


// Pressure corrector to enforce continuity
// Solution control is ignored here
tmp<fvScalarMatrix> Foam::atmTurbMesh::pEqn()
{
    volScalarField rAU(1.0/UEqn()->A());
    volVectorField HbyA
    (
      constrainHbyA(rAU*UEqn()->H(), this->U_, this->p_)
    );
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    adjustPhi(phiHbyA, U_, p_);

    tmp<volScalarField> rAtU(rAU);
    // Update the pressure BCs to ensure flux consistency
    // The incompressible form of constrainPressure() used which basically
    // overloaded the compressible form with rho being one field
    constrainPressure(p_, U_, phiHbyA, rAtU());
    if (pimple_.consistent())
    {
      rAtU = 1.0/max(1.0/rAU - tUEqn->H1(), 0.1/rAU);
      phiHbyA += 
          fvc::interpolate(rAtU() - rAU) * fvc::snGrad(p)*this->magSf();
    tmp<fvScalarMatrix> tpEqn
    (
      fvm::laplacian(rAtU(), this->p_) == fvc::div(phiHbyA)
    );

    return tpEqn; 
}
tmp<fvScalarMatrix> Foam::atmTurbMesh::solve_p()
{
  tmp<fvScalarMatrix> tpEqn(pEqn());
  Info << "solve_p() not implemented. " << endl;
  return tpEqn;
}

// Other variables
tmp<fvScalarMatrix> Foam::atmTurbMesh::TEqn()
{
  tmp<fvScalarMatrix> tTEqn
  (
      fvm::ddt(T_)
    + fvm::div(phi_, T_)
    == 
      // Warning: nut() is used instead of alphaEff T
      fvc::laplacian(this->turbulence_->nut(), T_)
  );
  return tTEqn;
}
    
void Foam::atmTurbMesh::solve_T()
{
  Info << "solve_theta() not implemented. " << endl;
  FatalError.exit();
}

tmp<fvScalarMatrix> Foam::atmTurbMesh::qEqn()
{
  tmp<fvScalarMatrix> tQEqn
  (
      fvm::ddt(q_)
    + fvm::div(phi_, q_)
    == 
      // Warning: nut() is used instead of alphaEff Q
      fvc::laplacian(this->turbulence_->nut(), q_)
  );
  return tQEqn;
}

void Foam::atmTurbMesh::solve_q()
{
  Info << "solve_q() not implemented. " << endl;
  FatalError.exit();
}
