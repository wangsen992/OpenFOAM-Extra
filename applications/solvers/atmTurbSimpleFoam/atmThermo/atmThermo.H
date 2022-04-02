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

Class
    Foam::atmThermo

Description

SourceFiles
    atmThermoI.H
    atmThermo.C
    atmThermoIO.C

\*---------------------------------------------------------------------------*/

#ifndef atmThermo_H
#define atmThermo_H

#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class atmThermo;
Istream& operator>>(Istream&, atmThermo&);
Ostream& operator<<(Ostream&, const atmThermo&);


/*---------------------------------------------------------------------------*\
                          Class atmThermo Declaration
\*---------------------------------------------------------------------------*/

class atmThermo
:
    public fluidThermo
{
    // Private Data

        //- Description of data_
        dataType data_;


    // Private Member Functions


public:

    // Static Data Members

        //- Static data staticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        atmThermo();

        //- Construct from components
        atmThermo(const dataType& data);

        //- Construct from Istream
        atmThermo(Istream&);

        //- Copy constructor
        atmThermo(const atmThermo&);

        //- Disallow default bitwise copy construction
        atmThermo(const atmThermo&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<atmThermo> New();


    //- Destructor
    ~atmThermo();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const atmThermo&);

        //- Disallow default bitwise assignment
        void operator=(const atmThermo&) = delete;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, atmThermo&);
        friend Ostream& operator<<(Ostream&, const atmThermo&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "atmThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
