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
    Base-class for canopy type sources with both fluid and solid properties. 

    Canopy type sources are similar to porous medium (actually is derived from
    porous medium) and dispersed phase flow (but different in the sense of 
    stationary/quasi-stationary solid phase). 

    This canopyModels package is built following the design of 
    `thermophysicalModels`. `canopy` classes (similar to `thermo` classes 
    such as `basicThermo`, `fluidThermo` etc) set up the general interface
    of the canopy properties itself, with increasing complexity (e.g. 
    basic->fluid through the addition of pressure), in canopy case, canopies
    can be interpreted as lad-based, or with greater grangularity and fully
    resovled leaves (in more advanced canopy model, inherited from basicCanopy)
    
SourceFiles
    canopyModel.C
\*---------------------------------------------------------------------------*/

#ifndef basicCanopy_H
#define basicCanopy_H


#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"
#include "cellSet.H"

// Include user defined HashTables
#include "cellSetHashTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class basicCanopy Declaration
\*---------------------------------------------------------------------------*/
class basicCanopy
{
protected:
    // Protected Static Member Functions
    

    // Protected Member Functions


public:
    // Public Classes
    class implementation;

    // Static Member data
      
        //- Name of the canopy properties dictionary
        static const word dictName;

    //- Runtime type information
    TypeName("basicCanopy");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicCanopy,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );
    
    // Static Member Functions
    
    // Selectors
    
        //- Generic New for each of the related canopy packages with
        // Uniform LAD and Cd
        template<class Canopy>
        static autoPtr<Canopy> New
        (
            const fvMesh&
        );
        
    //- Destructor
    virtual ~basicCanopy();
    
    // Member Functions

        //- Return the dictionary
        virtual const IOdictionary& properties() const = 0;

        //- Access the dictionary
        virtual IOdictionary& properties() = 0;

        //- Update properties
        virtual void correct() = 0;

        //- Access to canopy form-related variables

            //- Cells
            virtual const labelHashSet& canopyCells() const = 0;

            //- Leaf Area Density [m^-1]
            virtual const dimensionedScalarCellSet& lad() const = 0;

            //- Drag Coefficients
            virtual const scalarCellSet& Cd() const = 0;



};


/*---------------------------------------------------------------------------*\
                 Class basicThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class basicCanopy::implementation
:
    virtual public basicCanopy,
    public IOdictionary
{
protected:

    // Protected data


public:
            
    // Static Member data

        //- Name of the thermophysical properties dictionary
        using basicCanopy::dictName;


    // Constructors

        //- Construct from mesh
        implementation(const fvMesh&);
        

    //- Destructor
    virtual ~implementation();

    // Member Functions

        //- Return the dictionary
        virtual const IOdictionary& properties() const
        {
            return *this;
        }

        //- Access the dictionary
        virtual IOdictionary& properties()
        {
            return *this;
        }
        
};
}

#ifdef NoRepository
    #include "basicCanopyTemplates.C"
#endif

#endif
