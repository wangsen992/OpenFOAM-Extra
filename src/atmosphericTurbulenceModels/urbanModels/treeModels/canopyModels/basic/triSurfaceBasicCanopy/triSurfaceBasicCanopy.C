/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::triSurfaceBasicCanopy

Description
    A test implementation that does not include a drag model yet

SourceFiles
    triSurfaceBasicCanopy.C

\*---------------------------------------------------------------------------*/

#include "triSurfaceBasicCanopy.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicCanopy>
void Foam::triSurfaceBasicCanopy<BasicCanopy>::calculate()
{
    if (BasicCanopy::debug)
    {
        Info << "Running calculate..." << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicCanopy>
Foam::triSurfaceBasicCanopy<BasicCanopy>::triSurfaceBasicCanopy
(
    const fvMesh& mesh
)
:
    triSurfaceBasicCanopy<BasicCanopy>(mesh)
{
    calculate();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicCanopy>
void Foam::triSurfaceBasicCanopy<BasicCanopy>::correct()
{
    
    if (BasicCanopy::debug)
    {
        Info << "Running correct..." << endl;
    }

    calculate();

}
