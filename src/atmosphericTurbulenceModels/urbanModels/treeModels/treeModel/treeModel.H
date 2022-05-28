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
    Foam::treeModel

Description
    A base class for vegetation related models based on fvModel, such that 
    IO and interface with the equations can be handled through fvModels in 
    a centralized manner
    
SourceFiles
    treeModel.C
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvCellSet.H"
#include "canopyModel.H"
#include "trunkModel.H"

namespace Foam
{


class treeModel
{
    
    // Private Member Data
    
        //- Reference to components of trees
        canopyModel& canopy_;
        trunkModel& trunk_;

    // Private Member Functions
    
        //- Update the model
        virtual void update();


public: 
    
    //- Runtime type information
    TypeName("treeModel");

    // Constructors
    
        //- Constructor from components
        treeModel
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~treeModel(){}

    // Member Functions

        // Mesh motion

            //- Update for mesh changes
            virtual void updateMesh(const mapPolyMesh&);

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);

};

}
