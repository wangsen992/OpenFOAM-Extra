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
    Foam::triSurfaceDragCanopy

Description
    A test implementation that does not include a drag model yet

SourceFiles
    triSurfaceDragCanopy.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceBasicCanopy_H
#define triSurfaceBasicCanopy_H

#include "basicCanopy.H"
#include "triSurfaceCanopy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class triSurfaceBasicCanopy Declaration
\*---------------------------------------------------------------------------*/

template<class BasicDragCanopy>
class triSurfaceDragCanopy
:
    public triSurfaceCanopy<BasicDragCanopy>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

public:

    //- Runtime type information
    TypeName("triSurfaceDragCanopy");


    // Constructors

        //- Construct from mesh and phase name
        triSurfaceDragCanopy
        (
            const fvMesh&
        );

        //- Disallow default bitwise copy construction
        triSurfaceDragCanopy(const triSurfaceDragCanopy<BasicDragCanopy>&) = delete;


    //- Destructor
    virtual ~triSurfaceDragCanopy(){};
    
    // Member Functions

        //- Update properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const triSurfaceDragCanopy<BasicDragCanopy>&
        ) = delete;

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "triSurfaceDragCanopy.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
