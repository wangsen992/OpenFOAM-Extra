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
        (const fvMesh& mesh, const cellSet& cells),
        (mesh, cells)
    );
    
    // Static Member Functions
    
    // Selectors
    
        //- Generic New for each of the related canopy packages with
        // Uniform LAD and Cd
        template<class Canopy>
        static autoPtr<Canopy> New
        (
            const fvMesh&,
            const cellSet& cells
        );
        
        //- Generic New for each of the related thermodynamics packages with
        // dictionary (within which cellSet is specified, and more control 
        // for LAD and Cd)
        template<class Canopy>
        static autoPtr<Canopy> New
        (
            const fvMesh&,
            const dictionary&
        );

        //- To-Do: Other specializations of New
        
    
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

            //- Leaf Area Density [m^-1]
            virtual const dimensionedScalarCellSet& lad() = 0;

            //- Drag Coefficients
            virtual const scalarCellSet& Cd() = 0;

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

        //- Cell Set of Canopy Influence
        const cellSet& cells_;

        // Cell Set Fields

            //- Leaf Area Density [m^-1]
            dimensionedScalarCellSet lad_;

            //- Drag Coefficients
            scalarCellSet   Cd_;

public:
            
    // Static Member data

        //- Name of the thermophysical properties dictionary
        using basicCanopy::dictName;


    // Constructors

        //- Construct from mesh and cellSet
        implementation(const fvMesh&, const cellSet& cells);

        //- Construct from mesh and arbitrary dictionary
        implementation(const fvMesh&, const dictionary&);

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
        
        //- Access to canopy form-related variables

            //- Leaf Area Density [m^-1]
            virtual const dimensionedScalarCellSet& lad()
            {
                return lad_;
            }

            //- Drag Coefficients
            virtual const scalarCellSet& Cd()
            {
                return Cd_;
            }
};
}

#endif
