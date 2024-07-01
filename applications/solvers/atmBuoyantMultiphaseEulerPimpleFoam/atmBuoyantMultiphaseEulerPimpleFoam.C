/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    multiphaseEulerFoam

Description
    Solver for a system of any number of compressible fluid phases with a
    common pressure, but otherwise separate properties. The type of phase model
    is run time selectable and can optionally represent multiple species and
    in-phase reactions. The phase system is also run time selectable and can
    optionally represent different types of momentum, heat and mass transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "phaseSystem.H"
#include "phaseDynamicMomentumTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "fluidAtmThermo.H"
#include "atmHydrostaticInitialisation.H"
#include "referenceStateInitialisation.H"
#include "WRF.H"

#include "IOmanip.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    // #include "createFields.H"
    #include "createFieldsTest.H"
    #include "createFieldRefs.H"
    // #include "createDebugFields.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }


    Switch faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );
    Switch partialElimination
    (
        pimple.dict().lookupOrDefault<Switch>("partialElimination", false)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    runTime++;
    phaseModel& phase = fluid.movingPhases()[0];

    Info << "Initializing variables with WRF data" << endl;
    
    // Variables used from WRF: U, H2O.air, T.air, p, e.air, rho
    Info << "On U...." << endl;
    phase.URef() = wrf.U();
    phase.URef().correctBoundaryConditions();
    phase.phiRef() = fvc::flux(phase.URef());
    phase.alphaPhiRef() = fvc::flux(phase.URef());

    Info << "On T...." << endl;
    phase.thermoRef().T() = wrf.var("T.air");
    // phase.thermoRef().T().correctBoundaryConditions();

    phase.YRef()[0] = dimensionedScalar(dimless, 1) - wrf.var("H2O.air");
    phase.YRef()[1] = wrf.var("H2O.air");
    phase.YRef()[0].correctBoundaryConditions();
    phase.YRef()[1].correctBoundaryConditions();

    phase.thermoRef().p() = wrf.var("p");
    phase.thermoRef().p().correctBoundaryConditions();

    volScalarField& wrf_e(const_cast<volScalarField&>(wrf.var("e.air")));
    wrf_e = phase.thermoRef().he(phase.thermoRef().p(), phase.thermoRef().T());
    phase.thermoRef().he().primitiveFieldRef() = wrf.var("e.air");
    phase.thermoRef().he().correctBoundaryConditions();
    he = phase.thermoRef().he();

    phase.thermoRef().rho().primitiveFieldRef() = wrf.var("p") / (wrf.var("T.air") * dimensionedScalar(dimEnergy/(dimMass*dimTemperature), 287.05));
    forAll(phase.thermoRef().rho().boundaryFieldRef(), i)
    {
      phase.thermoRef().rho().boundaryFieldRef()[i] = phase.thermoRef().rho().boundaryFieldRef()[i].patchInternalField();
    }

    p_rgh = phase.thermoRef().p() - phase.thermoRef().rho() * gh - pRef;
    forAll(p_rgh.boundaryFieldRef(), i)
    {
      p_rgh.boundaryFieldRef()[i] = p_rgh.boundaryFieldRef()[i].patchInternalField();
    }
    
    // phase.thermoRef().correct();

    Info << "[Debug] average(U)= " << average(phase.URef()) << endl;
    Info << "[Debug] average(T)= " << average(phase.thermoRef().T()) << endl;
    Info << "[Debug] average(alphaPhi)= " << average(phase.alphaPhiRef()) << endl;
    Info << "[Debug] average(rho) = " << average(rho) << endl;
    Info << "[Debug] average(thermo.rho) = " << average(phase.thermoRef().rho()) << endl;
    
    Info << "Writetime after setting variables: " << runTime.value() << endl;
    // runTime.write();
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "createRDeltaTf.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"

        int nEnergyCorrectors
        (
            pimple.dict().lookupOrDefault<int>("nEnergyCorrectors", 1)
        );

        if (LTS)
        {
            #include "setRDeltaT.H"
            if (faceMomentum)
            {
                #include "setRDeltaTf.H"
            }
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // [Debug]
        phaseModel& phase = fluid.phases()[0];
        volScalarField qv_old("H2O", phase.Y("H2O"));
        volScalarField rho_v_old(qv_old * phase.thermo().rho());
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (!pimple.flow())
            {
                if (pimple.models())
                {
                    fvModels.correct();
                }

                if (pimple.thermophysics())
                {
                    fluid.solve(rAUs, rAUfs);
                    fluid.correct();
                    fluid.correctContinuityError();

                    #include "YEqns.H"

                    #include "EEqns.H"
                    #include "pEqnComps.H"

                    forAll(phases, phasei)
                    {
                        phases[phasei].divU(-pEqnComps[phasei] & p_rgh);
                    }
                }
            }
            else
            {
                if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
                {
                    // Store divU from the previous mesh so that it can be
                    // mapped and used in correctPhi to ensure the corrected phi
                    // has the same divergence
                    tmp<volScalarField> divU;

                    if
                    (
                        correctPhi
                    )
                    {
                        // Construct and register divU for mapping
                        divU = new volScalarField
                        (
                            "divU0",
                            fvc::div
                            (
                                fvc::absolute(phi, fluid.movingPhases()[0].U())
                            )
                        );
                    }

                    fvModels.preUpdateMesh();

                    mesh.update();

                    if (mesh.changing())
                    {
                        gh = (g & mesh.C()) - ghRef;
                        ghf = (g & mesh.Cf()) - ghRef;

                        fluid.meshUpdate();

                        if (correctPhi)
                        {
                            fluid.correctPhi
                            (
                                p_rgh,
                                divU,
                                pressureReference,
                                pimple
                            );
                        }

                        if (checkMeshCourantNo)
                        {
                            #include "meshCourantNo.H"
                        }
                    }
                }

                if (pimple.models())
                {
                    fvModels.correct();
                }

                fluid.solve(rAUs, rAUfs);
                fluid.correct();
                fluid.correctContinuityError();
                Info << "average(alphaPhi)= " << average(phase.alphaPhiRef()) << endl;

                // [Debug]
                qv_old = phase.Y("H2O");
                rho_v_old = qv_old * phase.thermo().rho();
                if (pimple.thermophysics())
                {
                    #include "YEqns.H"
                }

                if (faceMomentum)
                {
                    #include "pUf/UEqns.H"

                    if (pimple.thermophysics())
                    {
                        #include "EEqns.H"
                    }

                    #include "pUf/pEqn.H"
                }
                else
                {
                    Info << "[Debug] Header ----------------------" << endl;
                    #include "pU/UEqns.H"

                    if (pimple.thermophysics())
                    {
                        #include "EEqns.H"
                        Trad = phases[0].thermo().T();
                        Trad.correctBoundaryConditions();
                        
                    }

                    #include "pU/pEqnTest.H"
                }

                fluid.correctKinematics();

                if (pimple.turbCorr())
                {
                    fluid.correctTurbulence();
                }
                // fluid.correctEnergyTransport();

            }
        }    

        runTime.write();

          // // [Debug]
          // Info << setprecision(15);
          // phaseModel& otherPhase = fluid.phases()[1];
          // volScalarField qv("H2O", phase.Y("H2O"));
          // volScalarField rho_v(qv * phase * phase.thermo().rho());
          // volScalarField qt
          // (
          //   qv * phase * phase.thermo().rho() + otherPhase * otherPhase.thermo().rho()
          // );
          // Info << "[Debug] Phase1.pure : " << phase.pure() << " ; " 
          //      << "phase2.pure : " << otherPhase.pure() << endl;
          // 
          // 
          // Info  << "[totalWater] Time ,"     
          //       << runTime.timeName() << ", "                                       // << "qt = "                  
          //       << gAverage(phase) << ", "
          //       << gAverage(otherPhase) << ", "
          //       << gAverage(qt) << ", "                                                         // << "rho_v = "
          //       << gAverage(rho_v_old) << ", "                          // << ", rho_l = "               
          //       << gAverage(rho_v) << ", "                          // << ", rho_l = "               
          //       << gAverage((phase * phase.thermo().rho()).ref()) << ", "
          //       << gAverage((otherPhase * otherPhase.thermo().rho()).ref()) << ", "             // << ", dm.air = "              
          //       << gAverage((fluid.dmdts()[0] * phase.mesh().time().deltaT()).ref()) << ", "          // << ", dm.pos = "              
          //       << gAverage((fluid.dmdts()[1] * phase.mesh().time().deltaT()).ref()) << ", "          // << ", dm.pos = "              
          //       << gSum((neg(fluid.dmdts()[0])).ref()) << ", " // << ", dm.neg = "              
          //       << gAverage((posPart(fluid.dmdts()[0]) * phase.mesh().time().deltaT()).ref()) << ", " // << ", dm.neg = "              
          //       << gAverage((negPart(fluid.dmdts()[0]) * phase.mesh().time().deltaT()).ref()) << ", "// << "; sign = "                          << sign << " : " << pair 
          //       << gAverage((posPart(fluid.dmdts()[1]) * phase.mesh().time().deltaT()).ref()) << ", " // << ", dm.neg = "              
          //       << gAverage((negPart(fluid.dmdts()[1]) * phase.mesh().time().deltaT()).ref()) // << "; sign = "                          << sign << " : " << pair 
          //       << endl;

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
