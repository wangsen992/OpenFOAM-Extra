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


#include "IOdictionary.H"
#include "IOobject.H"
#include "basicCanopy.H"
#include "fvCellSet.H"

const Foam::word Foam::basicCanopy::dictName("canopyProperties");

namespace Foam
{
    defineTypeNameAndDebug(basicCanopy, 1);
    defineRunTimeSelectionTable(basicCanopy, fvMesh);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicCanopy::implementation::implementation
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
      IOobject
      (
        dictName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
      )
    ),
    mesh_(mesh)
{
    if (basicCanopy::debug)
    {
        Info << properties() << endl;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::basicCanopy::~basicCanopy() {};

Foam::basicCanopy::implementation::~implementation() {};
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
