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

Application
    testThermo

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidAtmThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList args(argc, argv);
    Time runTime(Time::controlDictName, args);
    fvMesh mesh
    (
      IOobject
      (
        fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
      )
    );

    Info << fvMesh::defaultRegion << endl;

    autoPtr<fluidAtmThermo> pthermo (fluidAtmThermo::New(mesh));
    Info << "Pointer generated." << endl;
    fluidAtmThermo& thermo = pthermo();
    Info << "Reference generated." << endl;
    Info << thermo.type() << endl;

    Info << pthermo->thermoName() << endl;
    Info << pthermo->phaseName() << endl;
    Info << pthermo->typeName_() << endl;

    pthermo->correct();
    Info << pthermo->type() << endl;
    // Info << pthermo->theta() << endl;
    // Info << pthermo->theta() * pthermo->rho() * pthermo->Cp() << endl;
    // Info << pthermo->he() << endl;
    Info << pthermo->Wi(0) << endl;
    Info << pthermo->species() << endl;
    Info << pthermo->Y(0) << endl;
    Info << pthermo->Y("dryAir") << endl;
    Info << pthermo->Y("H2O") / pthermo->Y("dryAir") << endl;


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
