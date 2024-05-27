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

void Foam::fv::WRFCoupler::updateVars(label it)
{
    wrf_.time().setTime(wrf_.dt() * it, it);
    U_ = wrf_.U(it);
    volScalarFieldPtrTable_["p"]() = wrf_.var("P", dimPressure, it) + wrf_.var("PB", dimPressure, it);

    volScalarFieldPtrTable_["T.air"]() 
      = (wrf_.var("T", dimTemperature, it) + wrf_.T0()) 
        * pow
          (
            volScalarFieldPtrTable_["p"]()/wrf_.P0(), 
            0.286
          );
    volScalarFieldPtrTable_["H2O.air"]() = wrf_.var("QVAPOR", dimless, it);
    volScalarFieldPtrTable_["thermo:rho.air"]() = volScalarFieldPtrTable_["p"]() / (volScalarFieldPtrTable_["T.air"]() * dimensionedScalar(dimEnergy/(dimMass*dimTemperature), 287.05));

    // Interpolate to the fields
    Info << "[WRFCoupler] Interpolating cell var values" << endl;
    projU_.primitiveFieldRef() = wrf_.interpolate(mesh_.cellCentres(), U_);
    
    projVolScalarFieldPtrTable_["T.air"]().primitiveFieldRef() = wrf_.interpolate(mesh_.cellCentres(), volScalarFieldPtrTable_["T.air"]());
    projVolScalarFieldPtrTable_["p"]().primitiveFieldRef() = wrf_.interpolate(mesh_.cellCentres(), volScalarFieldPtrTable_["p"]());
    projVolScalarFieldPtrTable_["e.air"]().primitiveFieldRef() = thermo_.he
      (
       projVolScalarFieldPtrTable_["p"](),
       projVolScalarFieldPtrTable_["T.air"]()
      );
    projVolScalarFieldPtrTable_["H2O.air"]().primitiveFieldRef() = wrf_.interpolate(mesh_.cellCentres(), volScalarFieldPtrTable_["H2O.air"]());
    projVolScalarFieldPtrTable_["thermo:rho.air"]().primitiveFieldRef() = wrf_.interpolate(mesh_.cellCentres(), volScalarFieldPtrTable_["thermo:rho.air"]());
    Info << "[WRFCoupler] Interpolation complete" << endl;
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
    U_
    (
      IOobject
      (
        "U.air",
        wrf_.time().timeName(),
        wrf_.time(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      wrf_.mesh(),
      dimVelocity
    ),
    projU_
    (
      IOobject
      (
        IOobject::groupName("U.air" , "proj"),
        mesh.time().timeName(),
        mesh.time(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      mesh,
      dimVelocity
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
      dimensionedScalar(dimless, 0)
    ),
    nestingDist_(dict.lookupOrDefault<scalar>("nestingDist", 500)),
    relaxationFactor_(dict.lookupOrDefault<scalar>("relaxationFactor", 20)),
    currTimeInd_(-1),
    phaseName_(word::null)
{
    Info << "WRF Loading starts" << endl;
    readCoeffs();

    // Set up the nesting cells
    for(word pn: std::vector<word>{"east", "west", "south", "north", "top"})
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
      cellWeights_[celli] = (nestingDist_ - nestingCellTbl_[celli])/nestingDist_;
    }

    // Create scalar hashtables for variables and cells
    volScalarFieldPtrTable_.set
    (
      "T.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            "T.air",
            wrf_.time().timeName(),
            wrf_.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          wrf_.mesh(),
          dimTemperature
        )
      )
    );
    volScalarFieldPtrTable_.set
    (
      "H2O.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            "H2O.air",
            wrf_.time().timeName(),
            wrf_.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          wrf_.mesh(),
          dimless
        )
      )
    );
    volScalarFieldPtrTable_.set
    (
      "thermo:rho.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            "thermo.rho.air",
            wrf_.time().timeName(),
            wrf_.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          wrf_.mesh(),
          dimDensity
        )
      )
    );
    volScalarFieldPtrTable_.set
    (
      "p",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            "p",
            wrf_.mesh().time().timeName(),
            wrf_.mesh().time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          wrf_.mesh(),
          dimPressure
        )
      )
    );
    
    // Create scalar hashtables for variables and cells
    projVolScalarFieldPtrTable_.set
    (
      "T.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            IOobject::groupName("T.air", "proj"),
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimTemperature
        )
      )
    );
    projVolScalarFieldPtrTable_.set
    (
      "H2O.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            IOobject::groupName("H2O.air", "proj"),
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimless
        )
      )
    );
    projVolScalarFieldPtrTable_.set
    (
      "thermo:rho.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            IOobject::groupName("thermo:rho.air", "proj"),
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimDensity
        )
      )
    );
    projVolScalarFieldPtrTable_.set
    (
      "p",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            IOobject::groupName("p", "proj"),
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimPressure
        )
      )
    );
    projVolScalarFieldPtrTable_.set
    (
      "e.air",
      autoPtr<volScalarField>
      (
        new volScalarField
        (
          IOobject
          (
            IOobject::groupName("e.air", "proj"),
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimEnergy / dimMass
        )
      )
    );

        
    // Update vars
    updateVars(0);

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::WRFCoupler::addSupFields() const
{
    return wordList{"U.air", "e.air", "H2O.air", "thermo:rho.air"};
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
      updateVars(newTimeInd);
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
  auto C = mesh().cellCentres();
  auto psi_foam = mesh().lookupObjectRef<psiType>(fieldName);
  tmp<volScalarField> deltaPsi = projVolScalarFieldPtrTable_[fieldName]() - psi_foam;

    eqn -= 0.1 * alpha * cellWeights_ * deltaPsi / mesh_.time().deltaT() * relaxationFactor_;
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
  tmp<volVectorField> deltaPsi = projU_ - psi_foam;

  // Remove vertical velocity addition
  std::for_each
  (
    deltaPsi.ref().begin(), 
    deltaPsi.ref().end(), 
    [](vector& v){v.z() = 0;}
  );
  eqn -= 0.1 * alpha * rho * cellWeights_ * deltaPsi / mesh_.time().deltaT() * relaxationFactor_ ;

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
  Info << psi_foam.dimensions() << endl;
  Info << projVolScalarFieldPtrTable_[fieldName]().dimensions();
  tmp<volScalarField> deltaPsi = projVolScalarFieldPtrTable_[fieldName]() - psi_foam;

    eqn -= 0.1 * alpha * rho * cellWeights_ * deltaPsi / mesh_.time().deltaT() * relaxationFactor_ ;
}

// ************************************************************************* //
