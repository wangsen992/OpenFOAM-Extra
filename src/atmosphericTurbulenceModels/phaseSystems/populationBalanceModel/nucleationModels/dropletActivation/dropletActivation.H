/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::diameterModels::nucleationModels::dropletActivation

Description
    Nucleation source due to droplet activation. Applicable for irreversible
    reactions only.


SourceFiles
    dropletActivation.C

\*---------------------------------------------------------------------------*/

#ifndef dropletActivation_H
#define dropletActivation_H

#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace nucleationModels
{

/*---------------------------------------------------------------------------*\
                       Class dropletActivation Declaration
\*---------------------------------------------------------------------------*/

class dropletActivation
:
    public nucleationModel
{
    // Private Data

        //- Sphere-equivalent nucleation diameter
        const dimensionedScalar dNuc_;

        //- Velocity group in which the nucleation occurs
        const velocityGroup& velGroup_;

        //- Phase in which the nucleation occurs
        const phaseModel& nucPhase_;

        //- PhasePair between which the nucleation occurs
        const phasePair pair_;

        //- Name of the specie which nucleates
        const word specieName_;

        //- vapor pressure
        volScalarField  e_;

        //- saturation vapor pressure
        volScalarField  es_;

        //- Supersaturation
        volScalarField  S_;


public:

    //- Runtime type information
    TypeName("dropletActivation");

    // Constructor

        dropletActivation
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~dropletActivation()
    {}


    // Member Functions

        //- Add to nucleationRate
        virtual void addToNucleationRate
        (
            volScalarField& nucleationRate,
            const label i
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
