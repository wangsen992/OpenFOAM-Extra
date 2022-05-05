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

#include "IOobject.H"
#include "atmTurbModel.H"

// Protected Member Function

// Constructor
Foam::atmTurbModel::atmTurbModel(IOobject io)
:
  mesh_(io),
  pimple_(mesh_),
  atmTurbDict_
  (
    IOobject
    ( "atmTurbulenceProperties",
      mesh_.time().constant(),
      mesh_,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  ),
  thermo_(fluidAtmThermo::New(mesh_)),
  rho0f_("rho0f", linearInterpolate(thermo_->rho0())),
  U_
  (
    IOobject
    (
        "U",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_
  ),
  p_rgh_
  (
    IOobject
    (
        "p_rgh",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_
  ),
  theta_
  (
    thermo_->theta()
  ),
  phi_
  (
    IOobject
    (
      "phi",
      mesh_.time().timeName(),
      mesh_,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
    ),
    fvc::flux(U_)
  ),
  rhophi_
  (
    IOobject
    (
      "rhophi",
      mesh_.time().timeName(),
      mesh_,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
    ),
    rho0f_ * phi_
  ),
  q_
  (
    thermo_->q()
  ),
  lwc_
  (
    thermo_->lwc()
  ),
  pressureReference_(thermo_->p(), pimple_.dict()),

  f_
  (
      "f", 
      atmTurbDict_.lookup("f")
  ),
  g_
  (
      "g", 
      atmTurbDict_.lookup("g")
  ),
  Ug_
  (
      "Ug", 
      atmTurbDict_.lookup("Ug")
  ),
  theta0_
  (
      "theta0", 
      atmTurbDict_.lookup("theta0")
  ),
  turbulence_
  (
    compressible::momentumTransportModel::New
    (
      thermo_->rho(),
      U_,
      rhophi_,
      thermo_()
    )
  ),
  transport_
  (
    fluidThermophysicalTransportModel::New(turbulence_, thermo_)
  ),
  fvModels_(fvModels::New(mesh_)),
  fvConstraints_(fvConstraints::New(mesh_)),
  UEqn_(),
  thetaEqn_(),
  qEqn_(),
  lwcEqn_()
{
  // creating fields
  mesh_.setFluxRequired(p_rgh_.name());
};
// Access Functions
