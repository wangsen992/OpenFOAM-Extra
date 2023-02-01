/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "Ostream.H"
#include "incompressibleKEqn.H"
#include "bound.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include <bits/types/locale_t.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicMomentumTransportModel>
void incompressibleKEqn<BasicMomentumTransportModel>::correctNut()
{
    scalar Ck = 0.10;
    this->nut_ = Ck * l_ *sqrt(k_);
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}

template<class BasicMomentumTransportModel>
void incompressibleKEqn<BasicMomentumTransportModel>::correctCeps()
{
    Ceps_ = 0.19 + (0.51 * l_ / delta_);
}

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> incompressibleKEqn<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
incompressibleKEqn<BasicMomentumTransportModel>::incompressibleKEqn
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    LESeddyViscosity<BasicMomentumTransportModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),

    Ceps_
    (
        IOobject
        (
            IOobject::groupName("Ceps", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimless
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    l_
    (
        IOobject
        (
            IOobject::groupName("l", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->delta()
    ), 

    delta_
    (
        IOobject
        (
            "delta",
            this->runTime_.constant(), // for non-adaptive mesh
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->delta()
    ),
    kTransportAdv_
    (
        IOobject
        (
            IOobject::groupName("kTransportAdv", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(), // for non-adaptive mesh
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimVelocity * dimVelocity * dimDensity / dimTime
        
    ),
    kTransportDiff_
    (
        IOobject
        (
            IOobject::groupName("kTransportDiff", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(), // for non-adaptive mesh
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimVelocity * dimVelocity * dimDensity / dimTime
    ),
    kProdShear_
    (
        IOobject
        (
            IOobject::groupName("kProdShear", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(), // for non-adaptive mesh
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimVelocity * dimVelocity * dimDensity / dimTime
    ),
    kProdDilate_
    (
        IOobject
        (
            IOobject::groupName("kProdDilate", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(), // for non-adaptive mesh
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimVelocity * dimVelocity * dimDensity / dimTime
    ),
    kDissipation_
    (
        IOobject
        (
            IOobject::groupName("kDissipation", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(), // for non-adaptive mesh
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimVelocity * dimVelocity * dimDensity / dimTime
    ),
    kSource_
    (
        IOobject
        (
            IOobject::groupName("kSource", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(), // for non-adaptive mesh
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimVelocity * dimVelocity * dimDensity / dimTime
    )
{
    bound(k_, this->kMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
    correctCeps();
    correctNut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool incompressibleKEqn<BasicMomentumTransportModel>::read()
{
    if (LESeddyViscosity<BasicMomentumTransportModel>::read())
    {
        //Ck_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> incompressibleKEqn<BasicMomentumTransportModel>::epsilon() const
{
    return volScalarField::New
    (
        IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
        this->Ce_*k()*sqrt(k())/this->delta()
    );
}


template<class BasicMomentumTransportModel>
void incompressibleKEqn<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    LESeddyViscosity<BasicMomentumTransportModel>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));
    volScalarField divalphaRhoPhi(fvc::div(fvc::absolute(alphaRhoPhi, U)));

    tmp<volTensorField> tgradU(fvc::grad(U));
    volScalarField G(this->GName(), nut*(tgradU() && dev(twoSymm(tgradU()))));
    tgradU.clear();
    
    // IO K equation contributors
    kTransportAdv_ = -fvc::div(alphaRhoPhi, k_);
    kTransportDiff_ = fvc::laplacian(alpha*rho*DkEff(), k_);
    
    kProdShear_ = alpha*rho*G;
    kProdDilate_ = - fvc::SuSp((2.0/3.0)*alpha*rho*divU, k_);
    kDissipation_ = - Ceps_ * alpha * rho * pow(k_, 1.5) / l_;
    kDissipation_.boundaryFieldRef() = 0;

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)  // convective transport of k
      - fvm::laplacian(alpha*rho*DkEff(), k_) // turbulent diffusion of k
     ==
        alpha*rho*G - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_) // production term
      - fvm::Sp(Ceps_*alpha*rho*sqrt(k_)/l_, k_) // dissipation
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);

    correctCeps();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
