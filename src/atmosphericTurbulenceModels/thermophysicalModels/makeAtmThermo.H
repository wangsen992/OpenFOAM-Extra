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

#ifndef makeAtmThermo_H
#define makeAtmThermo_H

#include "makeThermo.H"
#include "SpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineAtmThermo(BaseAtmThermo, CThermo, Mixture, ThermoPhys) \
                                                                               \
    typedef CThermo                                                            \
    <                                                                          \
        BaseAtmThermo::composite,                                         \
        SpecieMixture                                                          \
        <                                                                      \
            Mixture                                                            \
            <                                                                  \
                ThermoPhys                                                     \
            >                                                                  \
        >                                                                      \
    > CThermo##Mixture##ThermoPhys;                                            \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CThermo##Mixture##ThermoPhys,                                          \
        (                                                                      \
            word(CThermo##Mixture##ThermoPhys::typeName_()) + "<"              \
          + Mixture<ThermoPhys>::typeName() + ">"                              \
        ).c_str(),                                                             \
        0                                                                      \
    )


#define addAtmThermo(BaseAtmThermo, CThermo, Mixture, ThermoPhys)    \
                                                                               \
    addThermo(BaseAtmThermo, CThermo, Mixture, ThermoPhys)


#define makeAtmThermo(BaseAtmThermo, CThermo, Mixture, ThermoPhys)   \
                                                                               \
    defineAtmThermo(BaseAtmThermo, CThermo, Mixture, ThermoPhys);    \
                                                                               \
    addAtmThermo(fluidAtmThermo, CThermo, Mixture, ThermoPhys);      \
    addAtmThermo(BaseAtmThermo, CThermo, Mixture, ThermoPhys)


#define makeAtmThermos(                                                   \
    BaseThermo, BaseAtmThermo, CThermo, Mixture, ThermoPhys)              \
                                                                               \
    defineAtmThermo(BaseAtmThermo, CThermo, Mixture, ThermoPhys);    \
                                                                               \
    addAtmThermo(basicThermo, CThermo, Mixture, ThermoPhys);              \
    addAtmThermo(fluidThermo, CThermo, Mixture, ThermoPhys);              \
    addAtmThermo(BaseThermo, CThermo, Mixture, ThermoPhys);               \
    addAtmThermo(fluidAtmThermo, CThermo, Mixture, ThermoPhys);      \
    addAtmThermo(BaseAtmThermo, CThermo, Mixture, ThermoPhys)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
