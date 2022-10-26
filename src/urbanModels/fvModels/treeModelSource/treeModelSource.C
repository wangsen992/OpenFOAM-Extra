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

void Foam::fv::treeModelSource::readCoeffs()
{
    phaseName_ = coeffs().lookupOrDefault<word>("phase", word::null);

    UName_ =
        coeffs().lookupOrDefault<word>
        (
            "U",
            IOobject::groupName("U", phaseName_)
        );
}

bool Foam::fv::treeModelSource::validate_turb()
{
  return false;    
}

void Foam::fv::treeModelSource::apply
(
  const volScalarField& rho,
  fvMatrix<scalar>& eqn
)
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (eqn.psi().dimensions() == dimTemperature)
    {
        // eqn -= L/Cp*(fvc::ddt(rho, alpha1_));
    }
    else
    {
        // eqn -= L*(fvc::ddt(rho, alpha1_));
    }
}

Foam::treeModel Foam::fv::treeModelSource::constructTreeModel
(dictionary& dict)
{
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
    phaseName_(word::null),
    UName_(word::null),
    treeModelList_(1)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::treeModelSource::addSupFields() const
{
    
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>(basicThermo::dictName);

    return wordList({UName_, thermo.he().name()});
}


void Foam::fv::treeModelSource::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
}


void Foam::fv::treeModelSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
}


void Foam::fv::treeModelSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
}


bool Foam::fv::treeModelSource::read(const dictionary& dict)
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


// ************************************************************************* //
