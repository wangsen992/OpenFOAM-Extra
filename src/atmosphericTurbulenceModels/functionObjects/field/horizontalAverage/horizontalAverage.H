/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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
    Foam::functionObjects::horizontalAverage

Description
    Performance horizontal averaging on a strcutured mesh and 
    cache the averaged variables to memory and IO. 

    Must be performed on a rectangular mesh. 

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    horizontalAverage.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_horizontalAverage_H
#define functionObjects_horizontalAverage_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "HashSet.H"
//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class horizontalAverage Declaration
\*---------------------------------------------------------------------------*/
class horizontalAverage
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{

    // Private Member Functions
      
        // Private Data

          // List of variable names to perform horizontal averaing
          wordHashSet fieldSet_;

        //- File header information
        virtual void writeFileHeader(const label i);

public:

    //- Runtime type information
    TypeName("horizontalAverage");


    // Constructors

        //- Construct from Time and dictionary
        horizontalAverage
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        horizontalAverage(const horizontalAverage&) = delete;


    //- Destructor
    virtual ~horizontalAverage();
    
    // Member Functions

        virtual bool read(const dictionary&);

        //- Calculate the horizontalAverage field
        virtual bool execute();

        //- Write the horizontalAverage field
        virtual bool write();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const horizontalAverage&) = delete;
};
}
}

#endif
