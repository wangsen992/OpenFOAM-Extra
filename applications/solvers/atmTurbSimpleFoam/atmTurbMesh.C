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

// Constructor
Foam::atmTurbMesh::atmTurbMesh(IOobject io)
:
  fvMesh(io),
  atmTurbDict_
  (
    IOobject
    ( "atmTurbulenceProperties",
      this->time().constant(),
      *this,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  ),
  U_
  (
    IOobject
    (
        "U",
        this->time().timeName(),
        *this,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    *this
  ),
  p_
  (
    IOobject
    (
        "p",
        this->time().timeName(),
        *this,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    *this
  ),
  theta_
  (
    IOobject
    (
        "theta",
        this->time().timeName(),
        *this,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    *this
  ),
  q_
  (
    IOobject
    (
        "q",
        this->time().timeName(),
        *this,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    *this
  ),
  f_
  (
      "f", 
      this->atmTurbDict_.lookup("f")
  ),
  g_
  (
      "g", 
      this->atmTurbDict_.lookup("g")
  ),
  Ug_
  (
      "Ug", 
      this->atmTurbDict_.lookup("Ug")
  ),
  theta0_
  (
      "theta0", 
      this->atmTurbDict_.lookup("theta0")
  ),
  p0_
  (
      "p0", 
      this->atmTurbDict_.lookup("p0")
  ),
  rho0_
  (
      "rho0", 
      this->atmTurbDict_.lookup("rho0")
  ),
  phi_
  (
    IOobject
    (
      "phi",
      this->time().timeName(),
      *this,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
    ),
    fvc::flux(this->U_)
  ),
  thermo_(fluidThermo::New(*this)()),
  turbulence_
  (
    incompressible::momentumTransportModel::New
    (
      this->U_, this->phi_, this->thermo_
    )()
  )
{};
// Access Functions


