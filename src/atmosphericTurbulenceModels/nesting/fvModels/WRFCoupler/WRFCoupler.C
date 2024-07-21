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
    wrf_
    (
      mesh.time().lookupObjectRef<WRF>("WRF")
    ),
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
      dimensionedScalar(dimless, 0),
      "zeroGradient"
    ),
    nestingDist_(dict.lookupOrDefault<scalar>("nestingDist", 500)),
    nestingDistTop_(dict.lookupOrDefault<scalar>("nestingDistTop", 100)),
    relaxationFactor_(dict.lookupOrDefault<scalar>("relaxationFactor", 0.5)),
    currTimeInd_(-1),
    phaseName_(word::null)
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
    }

    nestingCells_.clear();
    nestingCellTbl_.clear();
    combineCloseCellTables(nestingCellTbl_, getPatchCloseCells(mesh, "top", nestingDistTop_));
    forAll(nestingCells_, i)
    {
      label celli = nestingCells_[i];
      cellWeights_[celli] = 1 - nestingCellTbl_[celli]/100.0;
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::WRFCoupler::addSupFields() const
{
    return wordList{"U.air", "e.air", "H2O.air", "dryAir.air","thermo:rho.air"};
    // return wordList{"U.air", "e.air", "H2O.air"};
}

void Foam::fv::WRFCoupler::correct()
{
    label newTimeInd
    (
      std::floor
      (
        mesh_.time().value() / wrf_.dt()
      )
    );

    if (newTimeInd > currTimeInd_)
    {
      Info << "Updating wrf variables" << endl;
      wrf_.updateVars(newTimeInd);
      currTimeInd_ = newTimeInd;
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
    fvMatrix<Foam::scalar>& eqn,
    const word& fieldName
) const
{

  Info << "[fvModel] addSup for var " << fieldName << endl;
  typedef GeometricField<Foam::scalar, fvPatchField, volMesh> psiType;
  auto V = mesh().V();
  auto psi_foam = mesh().lookupObjectRef<psiType>(fieldName);
  tmp<volScalarField> deltaPsi = wrf_.var(fieldName) - psi_foam;

  // eqn.source() -= 0.1 * (alpha * cellWeights_ * deltaPsi* relaxationFactor_)->field() 
  //                     * V.field();
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
  tmp<volVectorField> tdeltaPsi = wrf_.U() - eqn.psi();
  volVectorField& deltaPsi(tdeltaPsi.ref());

  // Remove vertical velocity addition
  // std::for_each
  // (
  //   deltaPsi.begin(), 
  //   deltaPsi.end(), 
  //   [](vector& v){v.z() = 0;}
  // );

  // Set patchField values
  forAll(deltaPsi.boundaryFieldRef(), i)
  {
    deltaPsi.boundaryFieldRef()[i] = deltaPsi.boundaryFieldRef()[i].patchInternalField();
  }

  // Computed smoothed field
  tmp<volVectorField> deltaPsiSmoothed = smooth<vector>(deltaPsi, 5);
  Info << "[fvModel] averageDeltaPsi = " << average(mag(deltaPsi)) << endl;
  // eqn.source() += 0.1 * (alpha * rho * cellWeights_ * deltaPsi* relaxationFactor_)->field()
  //                     * V.field();
  forAll(eqn.source(), i)
  {
    eqn.source()[i] -= 0.1 * (alpha[i] * rho[i] * cellWeights_[i] * (deltaPsi[i] - 0.2 * deltaPsiSmoothed.ref()[i]) * relaxationFactor_) * V[i];
  }
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
  const psiType& psi(wrf_.var(fieldName));
  tmp<volScalarField> tdeltaPsi = psi - eqn.psi();
  volScalarField& deltaPsi(tdeltaPsi.ref());
   
  // Set patchField values
  forAll(deltaPsi.boundaryFieldRef(), i)
  {
    deltaPsi.boundaryFieldRef()[i] = deltaPsi.boundaryFieldRef()[i].patchInternalField();
  }
  tmp<volScalarField> deltaPsiSmoothed = smooth<scalar>(deltaPsi, 5);
  Info << "[fvModel] averageDeltaPsi = " << average(mag(deltaPsi)) << endl;

  // eqn.source() += 0.1 * (alpha * rho * cellWeights_ * deltaPsi * relaxationFactor_)->field()
  //                     * V.field();
  forAll(eqn.source(), i)
  {
    eqn.source()[i] -= 0.1 * (alpha[i] * rho[i] * cellWeights_[i] * (deltaPsi[i] - 0.2 * deltaPsiSmoothed.ref()[i]) * relaxationFactor_) * V[i];
  }
  Info << "[fvModel] source added" << endl;
}

// ************************************************************************* //
