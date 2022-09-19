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
    Foam::heRhoAtmThermo

Description
    An analytical engine for the thermodynamics of a mixture in atm context

SourceFiles
    heRhoAtmThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heRhoAtmThermo_H
#define heRhoAtmThermo_H

#include "fvCFD.H"
#include "volFields.H"
#include "uniformGeometricFields.H"
#include "rhoAtmThermo.H"
#include "SpecieMixture.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heRhoAtmThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicRhoThermo, class MixtureType>
class heRhoAtmThermo
:
    public heThermo<BasicRhoThermo, MixtureType>
{
    // Private Member Functions

        // Validate basic mixture species
        // Species should include at least dryAir and H2O (moisture)
        void validate_mixture();

        //- Calculate the thermo variables
        void calculate();


public:

    //- Runtime type information
    TypeName("heRhoAtmThermo");


    // Constructors

        //- Construct from mesh and phase name
        heRhoAtmThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        heRhoAtmThermo(const heRhoAtmThermo<BasicRhoThermo, MixtureType>&) = delete;


    //- Destructor
    virtual ~heRhoAtmThermo();


    // Member Functions
        //- Exner Function
        //scalar exner(const scalar, const scalar);



            //- Update properties
            virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const heRhoAtmThermo<BasicRhoThermo, MixtureType>&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "heRhoAtmThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
