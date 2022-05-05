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
#include "fluidAtmThermo.H"
#include "dynamicMomentumTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "fluidThermophysicalTransportModel.H"
#include "fvModels.H"
#include "fvConstraints.H"


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

    // Instantiate for thermodynamics modelling
    autoPtr<fluidAtmThermo>     thermo_;
    surfaceScalarField rho0f_;

    // Fields (fields can be loaded from mesh directly)
    volVectorField U_; // Velocity Field
    volScalarField p_rgh_; // pressure perturbation
    volScalarField& theta_; // potential temperature
    surfaceScalarField phi_; // face velocity flux
    surfaceScalarField rhophi_; // face mass flux
    volScalarField& q_; // specific humidity 
    volScalarField& lwc_; // specific humidity 

    pressureReference pressureReference_;

    dimensionedVector f_; // coriolis parameter
    dimensionedVector g_; // gravititational acceleration
    dimensionedVector Ug_; // geostrohpic wind
    dimensionedScalar theta0_; // reference density

    // Instantiate for turbulence modeling (eddy viscosity)
    autoPtr<compressible::momentumTransportModel> turbulence_;
    // Instantiate for thermal turbulent modeling & scalar modeling (eddy
    // diffusivities)
    autoPtr<fluidThermophysicalTransportModel> transport_;
    fvModels& fvModels_;
    fvConstraints& fvConstraints_;

    tmp<fvVectorMatrix> UEqn_;
    tmp<fvScalarMatrix> thetaEqn_;
    tmp<fvScalarMatrix> qEqn_;
    tmp<fvScalarMatrix> lwcEqn_;

    // Protected Static Function
    static volScalarField& lookupOrConstructScalar
    (
        const fvMesh& mesh,
        const char* name
    );

    static volVectorField& lookupOrConstructVector
    (
        const fvMesh& mesh,
        const char* name
    );

    // Protected Member Function
      // exner function


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
    volScalarField& q(){return q_;}
    volScalarField& lwc(){return lwc_;}
    volScalarField& theta(){return theta_;}
    volScalarField& p_rgh(){return p_rgh_;}

    autoPtr<fluidAtmThermo> thermo(){return thermo_;}
    volScalarField& p(){return thermo_->p();}
    volScalarField& T(){return thermo_->T();}

    autoPtr<compressible::momentumTransportModel> turbulence(){return turbulence_;}

    // construct prognostic equations
    // Construct momentum equation
    tmp<fvVectorMatrix> UEqn();
    tmp<fvScalarMatrix> thetaEqn();
    tmp<fvScalarMatrix> qEqn();

    // Note: The implementation of liquid water content may be abstracted since
    // there are multiple approaches to this variable
    tmp<fvScalarMatrix> lwcEqn();

    void pressureCorrect();
    void nutCorrect();
    void phiCorrect();


    tmp<volVectorField> fU_Ug() {return f_ ^ (U_ - Ug_);}
};
}; // end of namespace Foam
