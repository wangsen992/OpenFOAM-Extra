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
    Foam::canopyModel

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

#ifndef canopyModel_H
#define canopyModel_H

#include "fvCFD.H"
#include "runTimeSelectionTables.H"
#include "canopySurfaceModel.H"
#include "triSurfaceMesh.H"
#include "runTimeSelectionTables.H"
#include "fluidAtmThermophysicalTransportModel.H"
#include "radiationModel.H"

#include "canopySurfaceModel.H"
#include "cellSetHashTables.H"

namespace Foam
{
class canopyModel
{
    
    // Private Member Data

        //- canopy surface model
        canopySurfaceModel   canopySurface_;
        
        //- Reference that provide access to transport
        const fluidAtmThermophysicalTransportModel&   transport_;

        const radiationModel& radiation_;

        //- internal source terms 
        dimensionedVectorCellSet fU_;
        dimensionedScalarCellSet fT_;
        dimensionedScalarCellSet fq_;

        dimensionedScalarCellSet fk_;
        dimensionedScalarCellSet feps_;
        dimensionedScalarCellSet fomega_;
        dimensionedTensorCellSet fR_;

    // Private Member Function
    
        //- Correct underlying model
        virtual void correctU() = 0;
        virtual void correctTurb() = 0;
        virtual void correctThermo() = 0;

public:
     
    // Declare run-time construtor selection table
    declareRunTimeNewSelectionTable
    (
      autoPtr, 
      canopyModel, 
      dictionary, 
      (
            canopySurfaceModel canopySurface,
            fluidAtmThermophysicalTransportModel& transport,
            radiationModel& radiation
      ), 
      (canopySurface, transport, radiation)
    );
    
    // Constructors
        //- Construct with components
        canopyModel
        (
            canopySurfaceModel canopySurface,
            const fluidAtmThermophysicalTransportModel& transport,
            const radiationModel& radiation
        );

    // Selectors

        //- Return areference to the selected canopy physics model
        static autoPtr<canopyModel> New
        (
            canopySurfaceModel canopySurface,
            const fluidAtmThermophysicalTransportModel& transport,
            const radiationModel& radiation
        );
            

    // Destructor
        virtual ~canopyModel(){};

    // Member Functions
        
        // Access
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
            dimensionedVectorCellSet& fU(){return fU_;};

            //- Source term to the TKE equation
            dimensionedScalarCellSet& fk(){return fk_;};

            //- Source term to the TKE Dissipation equation
            dimensionedScalarCellSet& feps(){return feps_;};

            //- Source term to the specific TKE dissipation equation
            dimensionedScalarCellSet& fomega(){return fomega_;};
            
            //- Source term to the Reynolds Stress equation
            dimensionedTensorCellSet& fR(){return fR_;} ;

            //- Source term to the T/theta equation
            dimensionedScalarCellSet& fT(){return fT_;};

            //- Source term to the specific humidity equation 
            dimensionedScalarCellSet& fq(){return fq_;};
        
        
        virtual void correct();
};

}
    
#endif
