/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
    Foam::potentialTemperature

Description
    Thermodynamics mapping class to expose the sensible enthalpy functions.

\*---------------------------------------------------------------------------*/

#ifndef potentialTemperature_H
#define potentialTemperature_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class sensibleEnthalpy Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class potentialTemperature
{

public:

    // Constructors

        //- Construct
        potentialTemperature()
        {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "potentialTemperature";
        }

        // Fundamental properties

            static word energyName()
            {
                return "h";
            }

            // Heat capacity at constant pressure [J/kg/K]
            scalar Cpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.Cp(p, T);
            }

            // Sensible enthalpy [J/kg]
            scalar THETA
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.theta(p, T);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
