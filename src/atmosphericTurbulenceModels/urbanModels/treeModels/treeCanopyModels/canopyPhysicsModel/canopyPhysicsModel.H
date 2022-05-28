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
    Base class for the handling the physics & evapotranspiration of leaves & canopies, 
    which provides source terms to the treeModel that ultimately interacts with 
    the solver. Only to be called by the treeModel. 

    The key is on how leaves spatial domain are defined. With OpenFOAM utilities,
    cellSet can be used (or should ultimately be used), but the given the 
    potential large number of instances, there could be many surface files, or 
    many dictionaries. So this selection should provide multiple methods and 
    potentially easier interface to set up the vegetation coverage. 
    
SourceFiles
    canopyModel.C
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "canopySurfaceModel.H"
#include "triSurfaceMesh.H"
#include "runTimeSelectionTables.H"
#include "fluidAtmThermophysicalTransportModel.H"
#include "radiationModel.H"

namespace Foam
{
class canopyPhysicsModel
{
    
    // Private Member Data

        //- Model name
        word      modelName_;

        //- Cells where source terms are added to
        canopySurfaceModel   canopySurface_;

        //- Reference that provide access to transport
        fluidAtmThermophysicalTransportModel&   transport_;

        radiationModel& radiation_;

    // Private Member Function

public:
     
    //- Runtime type information
    TypeName("canopyModel");
    
    // Constructors
        //- Construct with components
        canopyPhysicsModel
        (
            word    modelName,
            fvCellSet   set,
            fluidAtmThermophysicalTransportModel& transport,
            radiationModel& radiation
        );

    // Destructor
        virtual ~canopyModel(){};

    // Member Functions
        
        // Access
            //- Model Name
            inline const word modelName(){return modelName_;}

            //- fvCellSet
            inline const fvCellSet& set(){return canopyCells_;}

            //- thermophysical model 
            inline const fluidAtmThermophysicalTransportModel& transport(){return transport_;}

            //- radiation model
            inline const radiationModel& radiation(){return radiation_;}

        // Source terms
            //- Source term to the momentum equation
            //- Aerodynamic Drag of Tree Canopy (should be implemented as a
            //model) 
            virtual Field<vector> fU() const = 0;

            //- Source term to the TKE equation
            virtual Field<scalar> fk() const = 0;

            //- Source term to the TKE Dissipation equation
            virtual Field<scalar> feps() const = 0;

            //- Source term to the specific TKE dissipation equation
            virtual Field<scalar> fomega() const = 0;
            
            //- Source term to the Reynolds Stress equation
            virtual Field<tensor> fR() const = 0;

            //- Source term to the T/theta equation
            virtual Field<scalar> fT() const = 0;

            //- Source term to the specific humidity equation 
            virtual Field<scalar> fq() const = 0;
        
        //- Correct underlying model
        virtual void correct();
};

}
    
