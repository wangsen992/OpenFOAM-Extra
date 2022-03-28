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

template<class TurbulenceModel, class ThermoModel>
void Foam::atmTurbMesh<TurbulenceModel, ThermoModel>::createFields(fvMesh& mesh)
{
    
    const Time& runTime(mesh.time());
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
}

template<class TurbulenceModel, class ThermoModel>
void Foam::atmTurbMesh<TurbulenceModel, ThermoModel>::solveU()
{
    // Momentum predictor
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(this->phi_, this->U_)
      + this->turbulence_->divDevSigma(this->U_)
      + this->fU_Ug() // Geostrohpic Term
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    solve(UEqn == -fvc::grad(this->p_)) ;
}

template<class TurbulenceModel, class ThermoModel>
void Foam::atmTurbMesh<TurbulenceModel, ThermoModel>::solveP()
{
}

