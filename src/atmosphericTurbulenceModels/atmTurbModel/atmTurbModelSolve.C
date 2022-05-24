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
      tmp<fvVectorMatrix> tUEqn
      (
          fvm::ddt(U_)
        + fvm::div(phi_, U_)
        - fvm::laplacian(turbulence_->nuEff(), U_)
        == 
        fvModels_.source(U_)
      );
      tUEqn->relax();
      fvConstraints_.constrain(tUEqn.ref());

      if (pimple_.momentumPredictor())
      {
        solve(tUEqn() == -fvc::grad(p_rgh_) ); 
      }
      Info << "UEqn Constructed. " << endl;
      fvConstraints_.constrain(U_);
      UEqn_ = tUEqn;
      return UEqn_;
}

// Other variables
tmp<fvScalarMatrix> Foam::atmTurbModel::thetaEqn()
{
    Info << "Init TEqn." << endl;

    // Explicit coding of radiation term 
    tmp<volScalarField> tExner
    (
        thermo_->exner(thermo_->p(), thermo_->p0(), thermo_->gamma())
    );
    volScalarField& Exner = tExner.ref();
        
    dimensionedScalar rhoCp = average(thermo_->rho() * thermo_->Cp());
    tmp<volScalarField> tRadSourceT
    (
      - (radiation_->Rp() * pow4(theta_/Exner) / rhoCp)
    );

    volScalarField& radSourceT = tRadSourceT.ref();
    radSourceT.primitiveFieldRef() += radiation_->Ru() / rhoCp;

    tmp<volScalarField> radSourceTheta(radSourceT / Exner);

    // Solve theta equation
    tmp<fvScalarMatrix> tThetaEqn
    (
        fvm::ddt(theta_)
      + fvm::div(phi_, theta_)
      == 
        fvm::laplacian(transport_->alphaEff()/thermo_->rho(), theta_)
      + fvModels_.source(theta_)
      + fvc::Su(radSourceTheta, theta_)
    );
    tThetaEqn->relax();
    fvConstraints_.constrain(tThetaEqn.ref());
    tThetaEqn->solve();
    fvConstraints_.constrain(theta_);
    thetaEqn_ = tThetaEqn;
    return thetaEqn_;
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
        fvc::laplacian(transport_->alphaEff() / thermo_->rho(), q_)
      + fvModels_.source(q_)
    );
    fvConstraints_.constrain(tQEqn.ref());
    tQEqn->relax();
    tQEqn->solve();
    fvConstraints_.constrain(q_);
    qEqn_ = tQEqn;
    return qEqn_;
}

tmp<fvScalarMatrix> Foam::atmTurbModel::lwcEqn()
{
    Info << "Init lwcEqn." << endl;
    tmp<fvScalarMatrix> tlwcEqn
    (
        fvm::ddt(lwc_)
      + fvm::div(phi_, lwc_)
      == 
        // Warning: nut() is used instead of alphaEff Q
        fvc::laplacian(transport_->alphaEff() / thermo_->rho(), lwc_)
      + fvModels_.source(lwc_)
    );
    tlwcEqn->relax();
    fvConstraints_.constrain(tlwcEqn.ref());
    tlwcEqn->solve();
    fvConstraints_.constrain(lwc_);
    lwcEqn_ = tlwcEqn;
    return lwcEqn_;
}
