/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

#include "basicCanopy.H"
#include "wordIOList.H"
#include "compileTemplate.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
template<class Canopy>
Foam::autoPtr<Canopy> Foam::basicCanopy::New
(
    const fvMesh& mesh
)
{
    IOdictionary canopyDict
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    const word canopyType = canopyDict.lookup("canopyType");

    Info << "Selecting canopy model type " << canopyType << endl;

    typename Canopy::fvMeshConstructorTable::iterator cstrIter = 
        fvMeshConstructorTablePtr_->find(canopyType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown canopy model type "
            << canopyType << nl << nl 
            << "Valid canopy model types: " << endl
            << fvMeshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<basicCanopy>
    (
        cstrIter()(mesh)
    );
}

