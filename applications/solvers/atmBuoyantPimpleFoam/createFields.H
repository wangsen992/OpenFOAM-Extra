#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidAtmThermo> pThermo(fluidAtmThermo::New(mesh));
fluidAtmThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

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
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



#include "readpRef.H"


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

Info << "Reading pressure reference.\n" << endl;
pressureReference pressureReference
(
    p,
    p_rgh,
    pimple.dict(),
    thermo.incompressible()
);

mesh.setFluxRequired(p_rgh.name());

Info << "Reference state initialisation.\n" << endl;
referenceStateInitialisation
(
    rho,
    U,
    gh,
    ghf,
    pRef,
    thermo,
    pimple.dict()
);

Info << "Hydrostatic initialisation.\n" << endl;
atmHydrostaticInitialisation
(
    p_rgh,
    rho,
    U,
    gh,
    ghf,
    pRef,
    thermo,
    pimple.dict()
);
Info << "Update Reference p and rho.\n" << endl;
// Update reference values for p and rho
thermo.pRef() = p-p_rgh;
thermo.rhoRef() = thermo.rho();

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

// Phi creation is moved after hydrostatic initialization
#include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidAtmThermophysicalTransportModel> thermophysicalTransport
(
    fluidAtmThermophysicalTransportModel::New(turbulence(), thermo)
);

Info << "Correct nut and alphat based on prescribed U and T values.\n" << endl;
turbulence->correct();
thermophysicalTransport->correct();

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
#include "checkRadiationModel.H"

