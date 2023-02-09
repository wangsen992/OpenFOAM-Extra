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

#include "treeModelSource.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(treeModelSource, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        treeModelSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
dimensionedScalarCellSet Foam::fv::treeModelSource::getScalarSource(const word& fieldName) const
{
    dimensionedScalarCellSet Fi(tree_->canopy().canopyCells().size());
    if (fieldName == "k")
    {
        Fi.clear();
        Fi = tree_->canopy().Fturb("k");
    }
    else if (fieldName == "epsilon")
    {
        Fi.clear();
        Fi = tree_->canopy().Fturb("k");
    }

    return Fi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::treeModelSource::treeModelSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    mesh_(mesh),
    urbanDict_
    (
      IOobject
      (
        "urbanModelProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
      )
    ),
    tree_
    (
      treeModel::New
      (
        mesh, 
        urbanDict_.subDict("tree")
      )
    )
{
    Info << "[debug] treeModelSource Init: " << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::treeModelSource::addSupFields() const
{
    Info << "treeModelSource addSupFields: " << endl;
    wordList supFields = tree_->canopy().addSupFields();
    return supFields;
}

void Foam::fv::treeModelSource::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    Info << "treeModelSource addSup: " << endl;
    vectorField& Usource =  eqn.source();
    const dimensionedVectorCellSet& Fu = tree_->canopy().Fu();
    forAllConstIter(labelHashSet, tree_->canopy().canopyCells(), iter)
    {
        Usource[iter.key()] -= Fu[iter.key()].value() * mesh_.V()[iter.key()];
    }
}


void Foam::fv::treeModelSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    Info << "treeModelSource addRhoSup: " << endl;
    vectorField& Usource =  eqn.source();
    const dimensionedVectorCellSet& Fu = tree_->canopy().Fu();
    forAllConstIter(labelHashSet, tree_->canopy().canopyCells(), iter)
    {
        Usource[iter.key()] -= rho[iter.key()] * Fu[iter.key()].value() * mesh_.V()[iter.key()];
    }
}


void Foam::fv::treeModelSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    Info << "treeModelSource addAlphaRhoSup: " << endl;
    vectorField& Usource =  eqn.source();
    const dimensionedVectorCellSet& Fu = tree_->canopy().Fu();
    forAllConstIter(labelHashSet, tree_->canopy().canopyCells(), iter)
    {
        Usource[iter.key()] -= alpha[iter.key()] * rho[iter.key()] * Fu[iter.key()].value() * mesh_.V()[iter.key()];
    }
}

// Correcting all scalar equations
void Foam::fv::treeModelSource::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    Info << "treeModelSource addSupS: " << endl;
    scalarField& source =  eqn.source();
    dimensionedScalarCellSet Fi = getScalarSource(fieldName);
    forAllConstIter(labelHashSet, tree_->canopy().canopyCells(), iter)
    {
        source[iter.key()] -= Fi[iter.key()].value() * mesh_.V()[iter.key()];
    }
        
}


void Foam::fv::treeModelSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    Info << "treeModelSource addRhoSupC: " << endl;
    scalarField& source =  eqn.source();
    dimensionedScalarCellSet Fi = getScalarSource(fieldName);
    forAllConstIter(labelHashSet, tree_->canopy().canopyCells(), iter)
    {
        source[iter.key()] -= rho[iter.key()] * Fi[iter.key()].value() * mesh_.V()[iter.key()];
    }
}


void Foam::fv::treeModelSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    Info << "treeModelSource addAlphaRhoSupC: " << endl;
    scalarField& source =  eqn.source();
    dimensionedScalarCellSet Fi = getScalarSource(fieldName);
    forAllConstIter(labelHashSet, tree_->canopy().canopyCells(), iter)
    {
        source[iter.key()] -= alpha[iter.key()] * rho[iter.key()] * Fi[iter.key()].value() * mesh_.V()[iter.key()];
    }
}

void Foam::fv::treeModelSource::correct()
{
    Info << "TreeModelSource Correcting.." << endl;
    tree_->canopy().correctMomentumTransfer();
}

bool Foam::fv::treeModelSource::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
