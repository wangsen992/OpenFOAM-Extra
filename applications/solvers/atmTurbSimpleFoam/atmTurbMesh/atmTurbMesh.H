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
Class
    Foam::atmTurbMesh

Description
    
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "dynamicMomentumTransportModel.H"
#include "kinematicMomentumTransportModel.H"


namespace Foam
{

class atmTurbMesh
  : public fvMesh
{

protected:
    // Protected Data
    // Note: Similar to fvMesh that contains all relavant information of the
    // model setup (numerical schemes, solution control, variables etc), this
    // is the starting model that handles the fluid mechanics solver behaviour

    IOdictionary atmTurbDict_;
    autoPtr<fluidThermo>     thermo_;
    volScalarField rho_; // density
    volScalarField p_; // Perturbation Pressure to continuity

    // Fields (fields can be loaded from mesh directly)
    volVectorField U_; // Velocity Field
    surfaceScalarField phi_; // face flux
    volScalarField T_; // Potential Temperature
    volScalarField q_; // specific humidity 

    dimensionedVector f_; // coriolis parameter
    dimensionedVector g_; // gravititational acceleration
    dimensionedVector Ug_; // geostrohpic wind
    dimensionedScalar p0_; // reference density

    autoPtr<incompressible::momentumTransportModel> turbulence_;

  // Protected Static Function
  // Protected Member Function


public:

    // constructor
    atmTurbMesh
    (
        IOobject io
    );
    
    // Access Functions
    volVectorField U(){return this->U_;}

    // construct prognostic equations
    tmp<fvVectorMatrix> UEqn();
    tmp<fvScalarMatrix> pEqn();
    tmp<fvScalarMatrix> thetaEqn();
    tmp<fvScalarMatrix> qEqn();

    // member functions
    void solve_U();
    void solve_p();
    void solve_T();
    void solve_q();

    tmp<volVectorField> fU_Ug()
    {
      return f_ ^ (U_ - Ug_);
    }
    
};
}; // end of namespace Foam


