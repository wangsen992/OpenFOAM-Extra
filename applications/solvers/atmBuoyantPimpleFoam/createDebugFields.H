volVectorField phiAirDebugc
(
    IOobject
    (
        "phiAirDebugc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh, 
    dimVelocity * dimArea
);
volVectorField alphaRhoPhiAirDebugc
(
    IOobject
    (
        "alphaRhoPhiAirDebugc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh, 
    dimVelocity * dimArea * dimDensity 
);
volScalarField divphiAirDebug
(
    IOobject
    (
        "divphiAirDebug",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh, 
    dimVelocity * dimArea / dimVolume
);
volScalarField divalphaRhoPhiAirDebug
(
    IOobject
    (
        "divalphaRhoPhiAirDebug",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh, 
    dimVelocity * dimArea * dimDensity / dimVolume
);
