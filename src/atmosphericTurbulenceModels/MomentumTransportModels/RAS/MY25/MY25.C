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

#include "utils.H"
#include "MY25.H"
#include "bound.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void MY25<BasicMomentumTransportModel>::correctNut()
{
    NotImplemented;
    // this->nut_ = this->Cmu_*sqr(k_)/epsilon_;
    // this->nut_.correctBoundaryConditions();
    // fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> MY25<BasicMomentumTransportModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
MY25<BasicMomentumTransportModel>::MY25
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
    ReynoldsStress<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),

    A1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A1",
            this->coeffDict_,
            0.09
        )
    ),
    A2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A2",
            this->coeffDict_,
            0.09
        )
    ),
    B1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "B1",
            this->coeffDict_,
            0.09
        )
    ),
    B2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "B2",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.8
        )
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.5*tr(this->R_)
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    theta_(utils::lookupOrConstructScalar(this->mesh_, "theta"))
{
    if (type == typeName)
    {
        this->printCoeffs(type);

        this->boundNormalStress(this->R_);
        bound(epsilon_, this->epsilonMin_);
        k_ = 0.5*tr(this->R_);

        // Test theta loading (successful)
        // Info << "theta test: " << endl;
        // Info << average(theta_) << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool MY25<BasicMomentumTransportModel>::read()
{
    if (ReynoldsStress<RASModel<BasicMomentumTransportModel>>::read())
    {
        A1_.readIfPresent(this->coeffDict());
        A2_.readIfPresent(this->coeffDict());
        B1_.readIfPresent(this->coeffDict());
        B2_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
tmp<volSymmTensorField> MY25<BasicMomentumTransportModel>::DREff() const
{
    NotImplemented;
    // return volSymmTensorField::New
    // (
    //     "DREff",
    //     (Cs_*(this->k_/this->epsilon_))*this->R_ + I*this->nu()
    // );
}


template<class BasicMomentumTransportModel>
tmp<volSymmTensorField> MY25<BasicMomentumTransportModel>::DepsilonEff() const
{
    NotImplemented;
    // return volSymmTensorField::New
    // (
    //     "DepsilonEff",
    //     (Ceps_*(this->k_/this->epsilon_))*this->R_ + I*this->nu()
    // );
}


template<class BasicMomentumTransportModel>
void MY25<BasicMomentumTransportModel>::correct()
{
    NotImplemented;
    // if (!this->turbulence_)
    // {
    //     return;
    // }

    // // Local references
    // const alphaField& alpha = this->alpha_;
    // const rhoField& rho = this->rho_;
    // const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    // const volVectorField& U = this->U_;
    // volSymmTensorField& R = this->R_;
    // const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    // const Foam::fvConstraints& fvConstraints
    // (
    //     Foam::fvConstraints::New(this->mesh_)
    // );

    // ReynoldsStress<RASModel<BasicMomentumTransportModel>>::correct();

    // tmp<volTensorField> tgradU(fvc::grad(U));
    // const volTensorField& gradU = tgradU();

    // volSymmTensorField P(-twoSymm(R & gradU));
    // volScalarField G(this->GName(), 0.5*mag(tr(P)));

    // // Update epsilon and G at the wall
    // epsilon_.boundaryFieldRef().updateCoeffs();

    // // Dissipation equation
    // tmp<fvScalarMatrix> epsEqn
    // (
    //     fvm::ddt(alpha, rho, epsilon_)
    //   + fvm::div(alphaRhoPhi, epsilon_)
    //   - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
    //  ==
    //     Ceps1_*alpha*rho*G*epsilon_/k_
    //   - fvm::Sp(Ceps2_*alpha*rho*epsilon_/k_, epsilon_)
    //   + epsilonSource()
    //   + fvModels.source(alpha, rho, epsilon_)
    // );

    // epsEqn.ref().relax();
    // fvConstraints.constrain(epsEqn.ref());
    // epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    // solve(epsEqn);
    // fvConstraints.constrain(epsilon_);
    // bound(epsilon_, this->epsilonMin_);


    // // Correct the trace of the tensorial production to be consistent
    // // with the near-wall generation from the wall-functions
    // const fvPatchList& patches = this->mesh_.boundary();

    // forAll(patches, patchi)
    // {
    //     const fvPatch& curPatch = patches[patchi];

    //     if (isA<wallFvPatch>(curPatch))
    //     {
    //         forAll(curPatch, facei)
    //         {
    //             label celli = curPatch.faceCells()[facei];
    //             P[celli] *= min
    //             (
    //                 G[celli]/(0.5*mag(tr(P[celli])) + small),
    //                 1.0
    //             );
    //         }
    //     }
    // }

    // // Reynolds stress equation
    // tmp<fvSymmTensorMatrix> REqn
    // (
    //     fvm::ddt(alpha, rho, R)
    //   + fvm::div(alphaRhoPhi, R)
    //   - fvm::laplacian(alpha*rho*DREff(), R)
    //   + fvm::Sp(C1_*alpha*rho*epsilon_/k_, R)
    //   ==
    //     alpha*rho*P
    //   - (2.0/3.0*(1 - C1_)*I)*alpha*rho*epsilon_
    //   - C2_*alpha*rho*dev(P)
    //   + this->RSource()
    //   + fvModels.source(alpha, rho, R)
    // );

    // // Optionally add wall-refection term
    // if (wallReflection_)
    // {
    //     const volVectorField& n_(wallDist::New(this->mesh_).n());
    //     const volScalarField& y_(wallDist::New(this->mesh_).y());

    //     const volSymmTensorField reflect
    //     (
    //         Cref1_*R - ((Cref2_*C2_)*(k_/epsilon_))*dev(P)
    //     );

    //     REqn.ref() +=
    //         ((3*pow(Cmu_, 0.75)/kappa_)*(alpha*rho*sqrt(k_)/y_))
    //        *dev(symm((n_ & reflect)*n_));
    // }

    // REqn.ref().relax();
    // fvConstraints.constrain(REqn.ref());
    // solve(REqn);
    // fvConstraints.constrain(R);

    // this->boundNormalStress(R);

    // k_ = 0.5*tr(R);

    // correctNut();

    // // Correct wall shear-stresses when applying wall-functions
    // this->correctWallShearStress(R);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
