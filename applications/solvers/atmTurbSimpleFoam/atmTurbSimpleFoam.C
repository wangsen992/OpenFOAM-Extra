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
    atmTurbSimpleFoam

Description
    This is not a application. This is a testing script for the libraries
being developed right now. 

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "fvCFD.H"
#include "simpleControl.H"
#include "pressureReference.H"
#include "kinematicMomentumTransportModel.H"
#include "singlePhaseTransportModel.H"

#include "atmTurbModel/atmTurbModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // set root case
    argList::addBoolOption("dryRun", "if specified not solve");
    argList args(argc, argv);
    if (!args.checkRootCase())
    {
        FatalError.exit();
    }

    // create time
    Info << "Creat time" << nl << endl;
    Time runTime(Time::controlDictName, args);

    // create mesh
    Info  << "Create mesh for time = " 
          << runTime.timeName() 
          << nl << Foam::endl;

    atmTurbModel atm (
      IOobject
      (
        fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
      )
    );
    Info << "atmTurbMesh created" << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    if (args.optionFound("dryRun"))
    {
      Info << "Dry Run: " << nl << endl;
      return 0;
    }
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    while(atm.pimple().run(runTime))
    {
      runTime++;
      Info << "Time = " << runTime.timeName() << endl;
      while(atm.pimple().loop())
      {
        fvVectorMatrix UEqn = atm.UEqn()();
        while (atm.pimple().correct())
        {
          atm.pressureCorrect();
        }

        if (atm.pimple().turbCorr())
        {
          atm.nutCorrect();
        }
        atm.TEqn();


        dimensionedScalar rangeT
        (
          "rangeT",
          max(atm.T()) - min(atm.T())
        );
        Info << "rangeT = " << rangeT.value() << endl;

        dimensionedScalar contErr
        (
          "contErr", 
          runTime.deltaT0Value()*mag(fvc::div(atm.U()))().weightedAverage(atm.mesh().V())
        );
        Info << "contErr: " << contErr.value() << endl;
      }
      Info << "mean(U)= " << average(atm.U().primitiveField()) << endl;
      runTime.write();
    }

    Info<< "End\n" << endl;

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    
    return 0;
}

// ************************************************************************* //
