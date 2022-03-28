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
#include "kinematicMomentumTransportModel.H"


namespace Foam
{

template<class TurbulenceModel, class ThermoModel>
class atmTurbMesh
  : public fvMesh
{

protected:
  // Protected Data
    // Note: Similar to fvMesh that contains all relavant information of the
    // model setup (numerical schemes, solution control, variables etc), this
    // is the starting model that handles the fluid mechanics solver behaviour

    fvMesh& mesh_;

    // Fields
    volVectorField U_; // Velocity Field
    volScalarField p_; // Perturbation Pressure to continuity
    volScalarField theta_; // Potential Temperature
    volScalarField q_; // specific humidity 
    uniformDimensionedVectorField f_; // coriolis parameter
    uniformDimensionedVectorField Ug_; // geostrohpic wind
    uniformDimensionedScalarField rho_0_; // reference density

    surfaceScalarField phi_; // face flux
    TurbulenceModel& turbulence_;
    ThermoModel&     thermo_;

  // Protected Member Function
    static void createFields(fvMesh&);

public:

  // constructor
    atmTurbMesh
    (
      fvMesh& mesh,
      volVectorField& U,
      volScalarField& p,
      volScalarField& theta,
      volScalarField& q
    );
  
  // member functions

  void solveU();
  void solveP();

  tmp<volVectorField> fU_Ug()
  {
    return f_ ^ (U_ - Ug_);
  }
    
};

} // end of namespace Foam

