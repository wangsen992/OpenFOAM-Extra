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
    Foam::basicET

Description
    A basic evapotranspiration model to help set up the class templating 
    structure.

    Evapotranspiration dictates the transport of moisture from the foliage 
    and intake of heat from the evaporation of transpired water within trees. 
    Evapotranspiation class compute the source of moisture and heat at each 
    cell. Convective transfer is controled by Transport and the modulation of 
    ET rate is driven also by radiation.

\*---------------------------------------------------------------------------*/

#ifndef basicET_H
#define basicET_H

#include "word.H"
#include "scalar.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class heThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Transport, class Radiation>
class basicET
{

public:
    
    basicET(){};

    ~basicET(){};

    //Member Functions

        // Return the instantiated type name
        static word typeName()
        {
            return "basicET";
        }

        // Moisture release rate at cell [g m^-3 s^-1]
        scalar fq
        (
            const scalar T,
            const scalar q
        );

        // Heat absorption rate at cell due to ET [K s^-1]
        scalar fT
        (
            const scalar T,
            const scalar q
        );

};
}

#endif
