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

#ifndef canopyPhysicsModel_H
#define canopyPhysicsModel_H

#include "fvCFD.H"
#include "canopySurfaceModel.H"
#include "triSurfaceMesh.H"
#include "runTimeSelectionTables.H"
#include "fluidAtmThermophysicalTransportModel.H"
#include "radiationModel.H"

#include "canopySurfaceModel.H"

namespace Foam
{
class canopyPhysicsModel
{
    
    // Private Member Data

        //- Model name
        word      modelName_;

        //- canopy surface model
        canopySurfaceModel   canopySurface_;
        
        //- Reference that provide access to transport
        fluidAtmThermophysicalTransportModel&   transport_;

        radiationModel& radiation_;

        //- internal source terms 
        vectorField fU_;
        scalarField fT_;
        scalarField fq_;

        scalarField fk_;
        scalarField feps_;
        scalarField fomega_;
        tensorField fR_;

    // Private Member Function
    
        //- Correct underlying model
        virtual void correctU() = 0;
        virtual void correctTurb() = 0;
        virtual void correctThermo() = 0;

public:
     
    //- Runtime type information
    // TypeName("canopyModel");
    
    // Constructors
        //- Construct with components
        canopyPhysicsModel
        (
            word    modelName,
            canopySurfaceModel canopySurface,
            fluidAtmThermophysicalTransportModel& transport,
            radiationModel& radiation
        );

    // Destructor
        virtual ~canopyPhysicsModel(){};

    // Member Functions
        
        // Access
            //- Model Name
            inline const word modelName() const{return modelName_;}

            //- canopy surface model acess
            inline const canopySurfaceModel& canopySurface() const 
            {return canopySurface_;}

            //- cellLabels
            inline labelList cells()
            {return canopySurface_.canopyCells().toc();};

            //- thermophysical model 
            inline const fluidAtmThermophysicalTransportModel& transport() const {return transport_;}

            //- radiation model
            inline const radiationModel& radiation() const {return radiation_;}

        // Source terms
            //- Source term to the momentum equation
            //- Aerodynamic Drag of Tree Canopy (should be implemented as a
            //model) 
            vectorField& fU(){return fU_;};

            //- Source term to the TKE equation
            scalarField& fk(){return fk_;};

            //- Source term to the TKE Dissipation equation
            scalarField& feps(){return feps_;};

            //- Source term to the specific TKE dissipation equation
            scalarField& fomega(){return fomega_;};
            
            //- Source term to the Reynolds Stress equation
            tensorField& fR(){return fR_;} ;

            //- Source term to the T/theta equation
            scalarField& fT(){return fT_;};

            //- Source term to the specific humidity equation 
            scalarField& fq(){return fq_;};
        
        
        virtual void correct();
};

}
    
#endif
