/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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
    Foam::fv::treeModel

Description
    Interface to include tree models into the transport equations

Usage
    Example usage:
    \verbatim
    tree1
    {
        type        treeModel;

        U           U;  // Name of the velocity field
    }
    \endverbatim

SourceFiles
    buoyancyForce.C

\*---------------------------------------------------------------------------*/

#ifndef treeModel_H
#define treeModel_H

#include "fvModel.H"
#include "treeModel.H"

namespace Foam
{
namespace fv
{
/*---------------------------------------------------------------------------*\
                        Class treeModelSource Declaration
\*---------------------------------------------------------------------------*/

class treeModelSource
:
    public fvModel
{
    // Private Data
    
        //- Optional phase name
        word phaseName_;

        //- Name of the velocity field
        word UName_;

        //- Reference to embedded models
        thermophysicalTransportModel& transport_;
        radiationModel& radiation_;

        //- Pointer list of tree models to ingest into transport equation
        PtrList<treeModel> treeModelList_;

    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- construct tree model 
        treeModel constructTreeModel(dictionary& dict);

public:

    //- Runtime type information
    TypeName("treeModelSource");


    // Constructors

        //- Construct from explicit source name and mesh
        treeModelSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        treeModelSource(const treeModelSource&) = delete;
        
    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluate

            //- Add explicit contribution to incompressible momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const treeModelSource&) = delete;
};


}
}

#endif
