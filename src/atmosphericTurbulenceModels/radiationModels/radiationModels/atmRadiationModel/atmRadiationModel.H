/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::atmRadiationModel

Description
    Top level model for radiation modelling, adapted from official radiationModel
    to account for directional variances in various coefficients.

    This allows different coefficients for each radiativeIntensityRay to adjust 
    the coefficients according to the configurations in the absorption emission 
    models.

SourceFiles
    atmRadiationModel.C
    radiationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModel_H
#define radiationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fvMatricesFwd.H"
#include "Switch.H"
#include "fluidAtmThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

namespace atmRadiationModels
{
class atmAbsorptionEmissionModel;
class scatterModel;
}

/*---------------------------------------------------------------------------*\
                       Class atmRadiationModel Declaration
\*---------------------------------------------------------------------------*/

class atmRadiationModel
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- Reference to the temperature field
        const volScalarField& T_;

        //- Radiation model dictionary
        dictionary coeffs_;

        //- Radiation solver frequency - number flow solver iterations per
        //  radiation solver iteration
        label solverFreq_;

        //- Flag to enable radiation model to be evaluated on first iteration
        bool firstIter_;


        // References to the radiation sub-models

            //- Absorption/emission model
            autoPtr<atmRadiationModels::atmAbsorptionEmissionModel>
                absorptionEmission_;

            //- Scatter model
            autoPtr<atmRadiationModels::scatterModel> scatter_;

private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Initialise
        void initialise();


public:

    //- Runtime type information
    TypeName("atmRadiationModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            atmRadiationModel,
            T,
            (
                const volScalarField& T
            ),
            (T)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            atmRadiationModel,
            dictionary,
            (
                const dictionary& dict,
                const volScalarField& T
            ),
            (dict, T)
        );


    // Constructors

        //- Null constructor
        atmRadiationModel(const volScalarField& T);

        //- Construct from components
        atmRadiationModel(const word& type, const volScalarField& T);

        //- Construct from components
        atmRadiationModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& T
        );

        //- Disallow default bitwise copy construction
        atmRadiationModel(const atmRadiationModel&) = delete;


    // Selectors

        //- Return a reference to the selected radiation model
        static autoPtr<atmRadiationModel> New(const volScalarField& T);

        //- Return a reference to the selected radiation model
        static autoPtr<atmRadiationModel> New
        (
            const dictionary& dict,
            const volScalarField& T
        );


    //- Destructor
    virtual ~atmRadiationModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();

            //- Solve radiation equation(s)
            //- Source term component (for power of T^4)
            virtual void calculate() = 0;

            //- Read radiationProperties dictionary
            virtual bool read() = 0;


        // Access

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const = 0;
            
            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp(label LambdaI) const = 0;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const = 0;
            
            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru(label LambdaI) const = 0;

            //- Energy source term
            virtual tmp<fvScalarMatrix> Sh
            (
                const basicThermo& thermo,
                const volScalarField& he
            ) const;

            //- Temperature source term
            virtual tmp<fvScalarMatrix> ST
            (
                const dimensionedScalar& rhoCp,
                volScalarField& T
            ) const;

            //- Access to absorptionEmission model
            const atmRadiationModels::atmAbsorptionEmissionModel&
                absorptionEmission() const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const atmRadiationModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToRadiationRunTimeSelectionTables(model)                            \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        atmRadiationModel,                                                        \
        model,                                                                 \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        atmRadiationModel,                                                        \
        model,                                                                 \
        T                                                                      \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
