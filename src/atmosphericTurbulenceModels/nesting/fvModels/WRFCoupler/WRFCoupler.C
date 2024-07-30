/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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

#include "WRFCoupler.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include <algorithm>
#include <vector>

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(WRFCoupler, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        WRFCoupler,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::WRFCoupler::readCoeffs()
{
    phaseName_ = coeffs().lookupOrDefault<word>("phase", word::null);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::WRFCoupler::WRFCoupler
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    mesh_(mesh),
    wrfTime_
    (
      Time::controlDictName, 
      dict.lookup<string>("wrf_case_root"),
      dict.lookup<string>("wrf_case_name")
    ),
    pwrfMesh_(),
    thermo_
    (
      mesh.lookupObjectRef<fluidAtmThermo>
      (
        IOobject::groupName
        (
          "thermophysicalProperties",
          "air"
        )
      )
    ),
    varList_(dict.subDict("variables")),
    nestingCells_(),
    nestingCellTbl_(),
    cellWeights_
    ( 
      IOobject
      (
        "cellWeights",
        mesh.time().constant(),
        mesh.time(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar(dimless/dimTime, dict.lookupOrDefault<scalar>("nudgingCoeff", 0.0003)),
      "zeroGradient"
    ),
    nestingDist_(dict.lookupOrDefault<scalar>("nestingDist", 500)),
    nestingDistTop_(dict.lookupOrDefault<scalar>("nestingDistTop", 100)),
    relaxationFactor_(dict.lookupOrDefault<scalar>("relaxationFactor", 0.5)),
    currTimeInd_(-1),
    phaseName_(word::null),
    wrfi_(-1),
    projUold_
    (
      IOobject
      (
        IOobject::groupName("U.air" , "projold"),
        mesh.time().timeName(),
        mesh.time(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimVelocity,
      "zeroGradient"
    ),
    projUnew_
    (
      IOobject
      (
        IOobject::groupName("U.air" , "projnew"),
        mesh.time().timeName(),
        mesh.time(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimVelocity,
      "zeroGradient"
    )
{
    Info << "WRF Loading starts" << endl;
    readCoeffs();
    // Set up the nesting cells
    for(word pn: std::vector<word>{"east", "west", "south", "north"})
    {
      // nestingCells_.append(getPatchCloseCells(mesh, pn, nestingDist_).first);
      combineCloseCellTables(nestingCellTbl_ , getPatchCloseCells(mesh, pn, nestingDist_));
    }
    nestingCells_ = nestingCellTbl_.sortedToc();
    Info << "nestingCells size: " << nestingCells_.size() << endl;

    nestingCells_.resize(nestingCells_.size());
    nestingCellCentres_.resize(nestingCells_.size());
    std::transform
    (
      nestingCells_.cbegin(),
      nestingCells_.cend(),
      nestingCellCentres_.begin(),
      [&](label i){return mesh.cellCentres()[i];}
    );
    forAll(nestingCells_, i)
    {
      label celli = nestingCells_[i];
      cellWeights_[celli] = 1-nestingCellTbl_[celli]/nestingDist_;
      cellWeights_[celli] *= relaxationFactor_;
    }

    nestingCells_.clear();
    nestingCellTbl_.clear();
    combineCloseCellTables(nestingCellTbl_, getPatchCloseCells(mesh, "top", nestingDistTop_));
    forAll(nestingCells_, i)
    {
      label celli = nestingCells_[i];
      cellWeights_[celli] = 1 - nestingCellTbl_[celli]/nestingDistTop_;
      cellWeights_[celli] *= relaxationFactor_;
    }

    // Retrieve the pwrfMesh_
    pwrfMesh_.set
    (
      new fvMesh
      (
        IOobject
        (
          fvMesh::defaultRegion,
          wrfTime_.timeName(),
          wrfTime_,
          IOobject::MUST_READ
        )
      )
    );
    tsearcher_.set
    (
      new meshSearch(pwrfMesh_())
    );


    // Create scalar hashtables for variables and cells
    // Both new and old variables are selected, based on the relative position
    // of the current time in the wrf data timeseries
    Info << "creating old varialbes" << endl;
    projVolScalarFieldPtrTableOld_.set
    (
      "T.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            IOobject::groupName("T.air", "projold"),
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimTemperature,
          "zeroGradient"
        )
      )
    );
    projVolScalarFieldPtrTableOld_.set
    (
      "e.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            IOobject::groupName("e.air", "projold"),
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimEnergy / dimMass,
          "zeroGradient"
        )
      )
    );
    for(auto k : varList_.keys())
    {
      Info << k << endl;
      projVolScalarFieldPtrTableOld_.set
      (
        k,
        autoPtr<volScalarField>
        (
          new volScalarField
          (
            IOobject
            (
              IOobject::groupName(k, "projold"),
              mesh.time().timeName(),
              mesh.time(),
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
            ),
            mesh,
            varList_.subDict(k).lookup<dimensionSet>("dim"),
            "zeroGradient"
          )
        )
      );
    }

    Info << "creating new varialbes" << endl;
    projVolScalarFieldPtrTableNew_.set
    (
      "T.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            IOobject::groupName("T.air", "projnew"),
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimTemperature,
          "zeroGradient"
        )
      )
    );
    projVolScalarFieldPtrTableNew_.set
    (
      "e.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            IOobject::groupName("e.air", "projnew"),
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimEnergy / dimMass,
          "zeroGradient"
        )
      )
    );
    for(auto k : varList_.keys())
    {
      projVolScalarFieldPtrTableNew_.set
      (
        k,
        autoPtr<volScalarField>
        (
          new volScalarField
          (
            IOobject
            (
              IOobject::groupName(k, "projnew"),
              mesh.time().timeName(),
              mesh.time(),
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
            ),
            mesh,
            varList_.subDict(k).lookup<dimensionSet>("dim"),
            "zeroGradient"
          )
        )
      );
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::WRFCoupler::addSupFields() const
{
    return wordList{"U.air", "e.air", "H2O.air"};
    // return wordList{"U.air"};
}

void Foam::fv::WRFCoupler::correct()
{
    Info << "Loading wrf time" << endl;
    Time& wrfTime ( wrfTime_ );
    label curi
    (
      pwrfMesh_->time().findClosestTimeIndex
      (
        pwrfMesh_->time().times(),
        mesh().time().value()
      )
    );

    // Update old and new data for interpolation
    if (curi > wrfi_ )
    {
        Info << "Correcting WRF variables at t = " << curi << endl;
        wrfi_ = curi;
        wrfTime.setTime(wrfTime.times()[wrfi_], wrfi_);

        // load the wrfMesh varialbes
      {
        projUold_.primitiveFieldRef() = interpolate(mesh().points(), volVectorField( IOobject ( "U.air", wrfTime.timeName(), wrfTime, IOobject::MUST_READ), pwrfMesh_()));
        projVolScalarFieldPtrTableOld_["T.air"]->primitiveFieldRef() = interpolate(mesh().points(), volScalarField ( IOobject ( "T.air", wrfTime.timeName(), wrfTime, IOobject::MUST_READ), pwrfMesh_()));
        projVolScalarFieldPtrTableOld_["e.air"]->primitiveFieldRef() 
          = thermo_.he(thermo_.p(), projVolScalarFieldPtrTableOld_["T.air"]);
        for(auto k : projVolScalarFieldPtrTableOld_.toc())
        {
          if( k != "T.air" && k != "e.air" )
          {
            projVolScalarFieldPtrTableOld_[k]->primitiveFieldRef() = interpolate(mesh().points(), volScalarField ( IOobject ( k, wrfTime.timeName(), wrfTime, IOobject::MUST_READ), pwrfMesh_()));
          }
        }
      }
          
        // load the new ones
      {
        wrfTime.setTime(wrfTime.times()[wrfi_+1], wrfi_+1);
        projUnew_.primitiveFieldRef() = interpolate(mesh().points(), volVectorField( IOobject ( "U.air", wrfTime.timeName(), wrfTime, IOobject::MUST_READ), pwrfMesh_()));
        projVolScalarFieldPtrTableNew_["T.air"]->primitiveFieldRef() = interpolate(mesh().points(), volScalarField ( IOobject ( "T.air", wrfTime.timeName(), wrfTime, IOobject::MUST_READ), pwrfMesh_()));
        projVolScalarFieldPtrTableNew_["e.air"]->primitiveFieldRef() 
          = thermo_.he(thermo_.p(), projVolScalarFieldPtrTableNew_["T.air"]);
        for(auto k : projVolScalarFieldPtrTableNew_.toc())
        {
          if( k != "T.air" && k != "e.air" )
          {
            projVolScalarFieldPtrTableNew_[k]->primitiveFieldRef() = interpolate(mesh().points(), volScalarField ( IOobject ( k, wrfTime.timeName(), wrfTime, IOobject::MUST_READ), pwrfMesh_()));
          }
        }
      }
    }
}


bool Foam::fv::WRFCoupler::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}

void Foam::fv::WRFCoupler::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<Foam::vector>& eqn,
    const word& fieldName
) const
{
  Info << "[fvModel] adding wrf field " << fieldName << endl;
  typedef GeometricField<Foam::vector, fvPatchField, volMesh> psiType;
  auto psi_foam = mesh().lookupObjectRef<psiType>(fieldName);
  auto V = mesh().V();

  auto t = mesh().time().value();
  auto dt = mesh().time().deltaTValue();
  auto told = wrfTime_.times()[wrfi_].value();
  auto tnew = wrfTime_.times()[wrfi_+1].value();
  auto tprojPsi = (tnew - t)/(tnew-told) * projUold_ 
              +  (t - told)/(tnew-told) * projUnew_;
  volVectorField& projPsi(tprojPsi.ref());
  tmp<volVectorField> tdeltaPsi = projPsi - eqn.psi();
  volVectorField& deltaPsi(tdeltaPsi.ref());

  // // Remove vertical velocity addition
  // std::for_each
  // (
  //   deltaPsi.begin(), 
  //   deltaPsi.end(), 
  //   [](vector& v){v.z() = 0;}
  // );

  // // Set patchField values
  // forAll(deltaPsi.boundaryFieldRef(), i)
  // {
  //   deltaPsi.boundaryFieldRef()[i] = deltaPsi.boundaryFieldRef()[i].patchInternalField();
  // }

  // Computed smoothed field
  // tmp<volVectorField> deltaPsiSmoothed = smooth<vector>(deltaPsi, 5);
  // Info << "[fvModel] averageDeltaPsi = " << average(mag(deltaPsi)) << endl;
  // eqn.source() += 0.1 * (alpha * rho * cellWeights_ * deltaPsi* relaxationFactor_)->field()
  //                     * V.field();
  eqn += (alpha * rho * cellWeights_ * projPsi - fvm::Sp(alpha * rho * cellWeights_, eqn.psi()))/dt;
  // forAll(eqn.source(), i)
  // {
  //   eqn.source()[i] -= (alpha[i] * rho[i] * cellWeights_[i] * (deltaPsi[i] - 0.0 * deltaPsiSmoothed.ref()[i]) ) * V[i];
  // }
}

void Foam::fv::WRFCoupler::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<Foam::scalar>& eqn,
    const word& fieldName
) const
{
  Info << "[fvModel] addSup for var " << fieldName << endl;
  typedef GeometricField<Foam::scalar, fvPatchField, volMesh> psiType;
  auto psi_foam = mesh().lookupObjectRef<psiType>(fieldName);
  auto V = mesh().V();
  auto t = mesh().time().value();
  auto dt = mesh().time().deltaTValue();
  auto told = wrfTime_.times()[wrfi_].value();
  auto tnew = wrfTime_.times()[wrfi_+1].value();
  auto tprojPsi = (tnew - t)/(tnew-told) * projVolScalarFieldPtrTableOld_[fieldName]() 
              +  (t - told)/(tnew-told) * projVolScalarFieldPtrTableNew_[fieldName]();
  Info << told << ";" << tnew << endl;
  volScalarField& projPsi(tprojPsi.ref());
  tmp<volScalarField> tdeltaPsi = projPsi - eqn.psi();
  Info << told << ";" << tnew << endl;
  volScalarField& deltaPsi(tdeltaPsi.ref());
  Info << told << ";" << tnew << endl;
   
  // Set patchField values
  // forAll(deltaPsi.boundaryFieldRef(), i)
  // {
  //   deltaPsi.boundaryFieldRef()[i] = deltaPsi.boundaryFieldRef()[i].patchInternalField();
  // }
  // tmp<volScalarField> deltaPsiSmoothed = smooth<scalar>(deltaPsi, 5);
  Info << "[fvModel] t= " << mesh().time().value() << ", " << "averageDeltaPsi = " << average(mag(deltaPsi)) << endl;

  eqn += (alpha * rho * cellWeights_  * projPsi - fvm::Sp(alpha * rho * cellWeights_, eqn.psi()))/dt;

  // forAll(eqn.source(), i)
  // {
  //   eqn.source()[i] -= (alpha[i] * rho[i] * cellWeights_[i] * (deltaPsi[i] - 0.0 * deltaPsiSmoothed.ref()[i])) * V[i];
  // }
  Info << "[fvModel] source added" << endl;
}

// ************************************************************************* //
