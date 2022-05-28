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
    Foam::canopyPhysicsModel

Description
    Base class for the handling the evapotranspiration of leaves & canopies, 
    which provides source terms to the treeModel that ultimately interacts with 
    the solver. Only to be called by the treeModel. 

    The key is on how leaves spatial domain are defined. With OpenFOAM utilities,
    cellSet can be used (or should ultimately be used), but the given the 
    potential large number of instances, there could be many surface files, or 
    many dictionaries. So this selection should provide multiple methods and 
    potentially easier interface to set up the vegetation coverage. 

Description
    Base class for tree/plant canopy model based on bulk parametrization of the tree canopy with aerodynamic drag force. 
    
SourceFiles
    dragCanopyPhysicsModel.C
\*---------------------------------------------------------------------------*/

#include "canopyPhysicsModel.H"

namespace Foam
{

class dragCanopyPhysicsModel
  : virtual public canopyPhysicsModel
{
    // Private Member Data
        //- Spatially variable drag coefficient
        Field<scalar>   Cd_;

        //- Leaf area density (projected)
        Field<vector>   LAD_;

    // Private Member Functions
    
public:
  
    //- Runtime type information
    TypeName("dragCanopyPhysicsModel");

    // Constructors
        //- Construct with model name
        dragCanopyPhysicsModel
        (
            word    modelName,
            fvCellSet   set,
            fluidAtmThermophysicalTransportModel& transport,
            radiationModel& radiation
        );

    // Destructor
        virtual ~dragCanopyPhysicsModel(){};

    // Member Functions

        // Source terms
            //- Source term to the momentum equation
            //model) 
            virtual Field<vector> fU() const;

            //- Source term to the TKE equation
            virtual Field<scalar> fk() const;

            //- Source term to the TKE Dissipation equation
            virtual Field<scalar> feps() const;

            //- Source term to the Reynolds Stress equation
            virtual Field<tensor> fR() const;

            //- Source term to the T/theta equation
            virtual Field<scalar> fT() const;

            //- Source term to the specific humidity equation 
            virtual Field<scalar> fq() const;
        
        //- Correct underlying model
        virtual void correct();


};

}
