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
    Foam::thetaRhoThermo

Description
    A basic thermo class that includes potential temperature as a thermodynamic 
state variable (instead of a derived variable like "T" as defined in basicThermo)

SourceFiles
    thetaRhoThermoI.H
    thetaRhoThermo.C
    thetaRhoThermoIO.C

\*---------------------------------------------------------------------------*/

#ifndef thetaRhoThermo_H
#define thetaRhoThermo_H

#include "rhoThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class thetaRhoThermo;
Istream& operator>>(Istream&, thetaRhoThermo&);
Ostream& operator<<(Ostream&, const thetaRhoThermo&);


/*---------------------------------------------------------------------------*\
                          Class thetaRhoThermo Declaration
\*---------------------------------------------------------------------------*/

class thetaRhoThermo
:
    virtual public rhoThermo
{

    // Public Classes
    
        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;


    //- Runtime type information
    TypeName("thetaRhoThermo");


    //- Declare run-time constructor selection table
    declareRunTimeNewSelectionTable
    (
        autoPtr, 
        thetaRhoThermo, 
        fvMesh, 
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors
        
        //- Standard selection based on fvMesh
        static autoPtr<thetaRhoThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );

    //- Destructor
    virtual ~thetaRhoThermo();
        
    // Member Functions

        // Access to thermodynamic state variables

            //- Potential Temperature [K]
            //  Non-const access allowed for transport equations
            virtual volScalarField& theta() = 0;

            //- Potential Temperature [K]
            virtual const volScalarField& theta() const = 0;

            //- Potential Temperature for patch [K]
            virtual tmp<volScalarField> theta(const label patchi) const = 0;

            //- Potential Temperature for cell-set [K]
            virtual tmp<scalarField> theta
            (
                const labelList& cells
            ) const = 0;


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
        thetaRhoThermo();

        //- Construct from components
        thetaRhoThermo(const dataType& data);

        //- Construct from Istream
        thetaRhoThermo(Istream&);

        //- Copy constructor
        thetaRhoThermo(const thetaRhoThermo&);

        //- Disallow default bitwise copy construction
        thetaRhoThermo(const thetaRhoThermo&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<thetaRhoThermo> New();


    //- Destructor
    ~thetaRhoThermo();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const thetaRhoThermo&);

        //- Disallow default bitwise assignment
        void operator=(const thetaRhoThermo&) = delete;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, thetaRhoThermo&);
        friend Ostream& operator<<(Ostream&, const thetaRhoThermo&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thetaRhoThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
