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
      volScalarField& p_ = thermo_->p();
      volScalarField& T_ = thermo_->T();
      tmp<fvVectorMatrix> tUEqn
      (
          fvm::ddt(U_)
        + fvm::div(phi_, U_)
        - fvm::laplacian(turbulence_->nuEff(), U_)
        == 
        fvModels_.source(U_)
        - fU_Ug() 
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
        fvm::laplacian(transport_->alphaEff()/thermo_->rho0(), T_)
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
        fvc::laplacian(transport_->alphaEff() / thermo_->rho0(), q_)
    );
    return tQEqn;
}

