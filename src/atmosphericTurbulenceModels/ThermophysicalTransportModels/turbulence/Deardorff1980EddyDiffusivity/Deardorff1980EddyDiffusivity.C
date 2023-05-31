/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

// #include "Deardorff1980EddyDiffusivity.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class TurbulenceThermophysicalTransportModel>
void Deardorff1980EddyDiffusivity<TurbulenceThermophysicalTransportModel>::
correctl()
{
    // Load variable data into space
    const volScalarField& k(this->momentumTransport().k()());

    volVectorField& b
    (
      mesh_.lookupObjectRef<volVectorField>
      (
        IOobject::groupName
        (
          "b",
          this->momentumTransport().alphaRhoPhi().group()
        )
      )
    );
    // Test code
    b.correctBoundaryConditions();

    tmp<volScalarField> tdbdz
    (
      fvc::grad(b)->component(8)
    );
    const volScalarField& dbdz(tdbdz.ref());

    
    // Testing Stage
    forAll(dbdz, celli)
    {
      if (dbdz[celli] <= 0)
      {
         l_[celli] = delta_[celli];
      }
      else 
      {
          l_[celli] = 0.76 
            * sqrt(k[celli] / dbdz[celli]);

          if (l_[celli] > delta_[celli])
          {
              l_[celli] = delta_[celli];
          }
      }
    }
    l_.correctBoundaryConditions();
    Info << "[Deardorff1980EddyDiffusivity] max(l) = " << gMax(l_) << endl;
}

template<class TurbulenceThermophysicalTransportModel>
void Deardorff1980EddyDiffusivity<TurbulenceThermophysicalTransportModel>::
correctAlphat()
{
    // As boundary condition matter in this case, the original
    // code doesn't work in parallel for multi-phase case
    forAll(alphat_, i)
    {
        alphat_[i] = (1 + 2 * l_[i]/delta_[i])
                    * this->momentumTransport().rho()[i]
                    * this->momentumTransport().nut()()[i];
    }
    // alphat_ =
    //     (1 + 2 * l_ / delta_)
    //    *this->momentumTransport().rho()
    //    *this->momentumTransport().nut();
    alphat_.correctBoundaryConditions();
    Info << "[Deardorff1980EddyDiffusivity] min(delta) = " << gMin(delta_) << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
Deardorff1980EddyDiffusivity<TurbulenceThermophysicalTransportModel>::
Deardorff1980EddyDiffusivity
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    Deardorff1980EddyDiffusivity
    (
        typeName,
        momentumTransport,
        thermo,
        false
    )
{}


template<class TurbulenceThermophysicalTransportModel>
Deardorff1980EddyDiffusivity<TurbulenceThermophysicalTransportModel>::
Deardorff1980EddyDiffusivity
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo,
    const bool allowDefaultPrt
)
:
    TurbulenceThermophysicalTransportModel
    (
        type,
        momentumTransport,
        thermo
    ),
    mesh_(momentumTransport.mesh()),
    alphat_
    (
        IOobject
        (
            IOobject::groupName
            (
                "alphat",
                this->momentumTransport().alphaRhoPhi().group()
            ),
            momentumTransport.time().timeName(),
            momentumTransport.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        momentumTransport.mesh()
    ),
    l_
    (
      mesh_.lookupObjectRef<volScalarField>
      (
        IOobject::groupName
        (
            "l",
            this->momentumTransport().alphaRhoPhi().group()
        )
      )
    ),
    delta_
    (
      mesh_.lookupObjectRef<volScalarField>("delta")
    )
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
bool Deardorff1980EddyDiffusivity<TurbulenceThermophysicalTransportModel>::read()
{
    if (TurbulenceThermophysicalTransportModel::read())
    {
        // Prt_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class TurbulenceThermophysicalTransportModel>
tmp<surfaceScalarField>
Deardorff1980EddyDiffusivity<TurbulenceThermophysicalTransportModel>::q() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName
        (
            "q",
            this->momentumTransport().alphaRhoPhi().group()
        ),
       -fvc::interpolate(this->alphaEff()*this->alpha())
       *fvc::snGrad(this->thermo().he())
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
Deardorff1980EddyDiffusivity<TurbulenceThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
    return -fvm::laplacian(this->alpha()*this->alphaEff(), he);
}


template<class TurbulenceThermophysicalTransportModel>
tmp<surfaceScalarField>
Deardorff1980EddyDiffusivity<TurbulenceThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
    return surfaceScalarField::New
    (
        IOobject::groupName
        (
            "j(" + Yi.name() + ')',
            this->momentumTransport().alphaRhoPhi().group()
        ),
       -fvc::interpolate(this->DEff(Yi)*this->alpha())*fvc::snGrad(Yi)
    );
}


template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
Deardorff1980EddyDiffusivity<TurbulenceThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    return -fvm::laplacian(this->alpha()*this->DEff(Yi), Yi);
}


template<class TurbulenceThermophysicalTransportModel>
void Deardorff1980EddyDiffusivity<TurbulenceThermophysicalTransportModel>::
correct()
{
    TurbulenceThermophysicalTransportModel::correct();
    correctl();
    correctAlphat();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
