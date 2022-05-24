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
\*---------------------------------------------------------------------------*/

#include "treeModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //
namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(treeModel, 0);
    addToRunTimeSelectionTable(fvModel, treeModel, dictionary);
}
}

Foam::wordList Foam::fv::treeModel::addSupFields() const
{
    return wordList({UName_, TName_, qName_});
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::treeModel::readCoeffs()
{
    UName_ = coeffs().lookupOrDefault<word>("UName_", "U");
    TName_ = coeffs().lookupOrDefault<word>("TName_", "T");
    qName_ = coeffs().lookupOrDefault<word>("qName_", "q");
}

void Foam::fv::treeModel::update()
{}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::fv::treeModel::treeModel
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    UName_(word::null),
    TName_(word::null),
    qName_(word::null)
{
    readCoeffs();
}
    

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// ************************************************************************* //
