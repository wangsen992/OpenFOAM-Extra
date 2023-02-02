#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating phaseSystem\n" << endl;

autoPtr<phaseSystem> fluidPtr
(
    phaseSystem::New(mesh)
);
phaseSystem& fluid = fluidPtr();
phaseSystem::phaseModelList& phases = fluid.phases();

#include "gh.H"

volScalarField& p = phases[0].thermoRef().p();

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

if (fluid.found("pMin"))
{
    IOWarningInFunction(fluid)
        << "Pressure limits, pMin and pMax, are now read from "
        << pimple.dict().name() << endl;
}

#include "readpRef.H"

pressureReference pressureReference
(
    p,
    p_rgh,
    pimple.dict(),
    fluid.incompressible()
);

if (fluid.incompressible())
{
    Info << "Running p init 1" << endl;
    p = p_rgh + fluid.rho()*gh;
}

if (p_rgh.needReference() && fluid.incompressible())
{
    Info << "Running p init 2" << endl;
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
}

mesh.setFluxRequired(p_rgh.name());

// thermo initialisation
forAll(fluid.movingPhases(), movingPhasei)
{
    phaseModel& phase = fluid.movingPhases()[movingPhasei];

    if (phase.name() == "air")
    {

        volVectorField U_air
        (
            IOobject
            (
                "U_air",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phase.U()
        );

        // Cast to custom fluidThermo type for init
        Info << phase.type() << endl;
        fluidAtmThermo & thermo
        (
          dynamic_cast<fluidAtmThermo&>(phase.thermoRef())
        );

        volScalarField rho_air
        (
            IOobject
            (
                "rhoairTest",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            phase.thermo().rho()
        );


        Info << "Reference state initialisation on phase " 
             << phase.name() 
             << "\n" << endl;

        referenceStateInitialisation
        (
            rho_air,
            phase.U(),
            gh,
            ghf,
            pRef,
            thermo,
            pimple.dict()
        );

        Info << "Hydrostatic initialisation on phase " 
             << phase.name() 
             << "\n" << endl;

        atmHydrostaticInitialisation
        (
            p_rgh,
            rho_air,
            phase.U(),
            gh,
            ghf,
            pRef,
            thermo,
            pimple.dict()
        );

      }    
}

PtrList<volScalarField> rAUs;
PtrList<surfaceScalarField> rAUfs;

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.rho()
);
