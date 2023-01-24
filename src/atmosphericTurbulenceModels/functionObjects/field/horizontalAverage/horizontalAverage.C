/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "horizontalAverage.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "objectRegistryFuncs.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(horizontalAverage, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        horizontalAverage,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::horizontalAverage::writeFileHeader(const label i)
{
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::horizontalAverage::horizontalAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    fieldSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::horizontalAverage::~horizontalAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::horizontalAverage::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    if (dict.found("field"))
    {
        fieldSet_.insert(word(dict.lookup("field")));
    }
    else
    {
        fieldSet_.insert(wordList(dict.lookup("fields")));
    }

    if (fieldSet_.size())
    {
        Info<< "Horizontally averaging fields:" << nl;
        forAllConstIter(wordHashSet, fieldSet_, iter)
        {
            Info<< "    "
                << IOobject::groupName(iter.key(), "mean") << nl;
        }
        Info<< endl;
    }
    else
    {
        Info<< "no fields requested to be stored" << nl << endl;
    }

    return true;
}


bool Foam::functionObjects::horizontalAverage::execute()
{
    forAllConstIter(wordHashSet, fieldSet_, iter)
    {
        Info << "Testing on variable: " << iter.key() << endl;
        //if (mesh_.objectRegistry::foundObject<volScalarField>(iter.key().c_str()))
        //{
            Info << "Performing averaging on field: " << iter.key() << endl;

            // Load field ref from objectRegistry
            volScalarField& rawField
            (
              Foam::lookupOrConstruct
              (
                mesh_, 
                iter.key().c_str()
              )
            );
            volScalarField meanField
            (
              IOobject
              (
                IOobject::groupName(iter.key(),"mean").c_str(),
                obr_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
              ),
              mesh_,
              rawField.dimensions()
            );
                

            scalarField spanavgField
            (
                meshIndexer_.collapse(rawField)
            );

            Info << "No. count meanField: " << meanField.size() << endl;
            Info << "No. count spanavgField: " << spanavgField.size() << endl;

            meanField.field() = spanavgField;
            meanField.write();
        //}
        //else
        //{
         //   Info << "Didn't find variable: " << iter.key() << endl;
        //}
    }
    return true;
}


bool Foam::functionObjects::horizontalAverage::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    Log << endl;

    return true;
}


// ************************************************************************* //
