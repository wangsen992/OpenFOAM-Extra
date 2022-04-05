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
#include "pimpleControl.H"
#include "pressureReference.H"


namespace Foam
{

class atmTurbModel
{

protected:
    // Protected Data
    // Note: Similar to fvMesh that contains all relavant information of the
    // model setup (numerical schemes, solution control, variables etc), this
    // is the starting model that handles the fluid mechanics solver behaviour

    
    fvMesh  mesh_;
    pimpleControl pimple_;
    IOdictionary atmTurbDict_;
    autoPtr<fluidThermo>     thermo_;

    // Fields (fields can be loaded from mesh directly)
    volVectorField U_; // Velocity Field
    surfaceScalarField phi_; // face flux
    volScalarField q_; // specific humidity 

    pressureReference pressureReference_;

    dimensionedVector f_; // coriolis parameter
    dimensionedVector g_; // gravititational acceleration
    dimensionedVector Ug_; // geostrohpic wind
    dimensionedScalar p0_; // reference density
    dimensionedScalar T0_; // reference density

    autoPtr<incompressible::momentumTransportModel> turbulence_;

    tmp<fvVectorMatrix> UEqn_;
    tmp<fvScalarMatrix> TEqn_;
    tmp<fvScalarMatrix> qEqn_;

    // Protected Static Function
    // Protected Member Function
      // exner function
      volScalarField exner (){ return pow((thermo_->p()/p0_), thermo_->gamma()())();}


public:

    // constructor
    atmTurbModel
    (
        IOobject io
    );
    
    // Access Functions
    fvMesh& mesh(){return mesh_;}
    const Time& runTime(){return mesh_.time();}
    pimpleControl& pimple(){return pimple_;}
    volVectorField& U(){return U_;}
    volScalarField& p(){return thermo_->p();}
    volScalarField& T(){return thermo_->T();}
    volScalarField& q(){return q_;}
    autoPtr<fluidThermo> thermo(){return thermo_;}
    autoPtr<incompressible::momentumTransportModel> turbulence(){return turbulence_;}

    // construct prognostic equations
    // Construct momentum equation
    tmp<fvVectorMatrix> UEqn();
    tmp<fvScalarMatrix> TEqn();
    tmp<fvScalarMatrix> qEqn();

    void pressureCorrect();
    void nutCorrect();
    void phiCorrect();


    tmp<volVectorField> fU_Ug() {return f_ ^ (U_ - Ug_);}
};
}; // end of namespace Foam
