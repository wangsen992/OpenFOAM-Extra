/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "Ostream.H"
#include "referenceStateInitialisation.H"

#include "fluidAtmThermo.H"
#include "fvmLaplacian.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "constrainPressure.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::referenceStateInitialisation
(
    volScalarField& rho,
    const volVectorField& U,
    const volScalarField& gh,
    const surfaceScalarField& ghf,
    const uniformDimensionedScalarField& pRef,
    fluidAtmThermo& thermo,
    const dictionary& dict
)
{
    if (dict.lookupOrDefault<bool>("referenceStateInitialisation", false))
    {
        Info << "Starting referenceStateInitialisation.\n" << endl;
        const fvMesh& mesh = U.mesh();

        volScalarField& ph_rgh = regIOobject::store
        (
            new volScalarField
            (
                IOobject
                (
                    "ph_rgh",
                    "0",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );

        if (!mesh.time().restart())
        {
            volScalarField& p = thermo.p();
            p = ph_rgh + rho*gh + pRef;
            // Enforce temperature in thermo to compute hydrostatic reference
            // state
            volScalarField T_orig = thermo.T();
            volScalarField he_orig = thermo.he();
            dimensionedScalar Tb = average(T_orig);
            //- Enforce temperature and energy profile
            // Note energy must be enforced as temperature is 
            // derived from temperature. 
            volScalarField& T = thermo.T();
            volScalarField& he = thermo.he();
            T = Tb + gh / thermo.Cp();
            he = thermo.he(p, T);

            // Update steps
            thermo.correct();
            rho = thermo.rho();
            label nCorr
            (
                dict.lookupOrDefault<label>("nHydrostaticCorrectors", 5)
            );

            for (label i=0; i<nCorr; i++)
            {
                surfaceScalarField rhof("rhof", fvc::interpolate(rho));

                surfaceScalarField phig
                (
                    "phig",
                    -rhof*ghf*fvc::snGrad(rho)*mesh.magSf()
                );

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(ph_rgh, rho, U, phig, rhof);

                fvScalarMatrix ph_rghEqn
                (
                    fvm::laplacian(rhof, ph_rgh) == fvc::div(phig)
                );

                ph_rghEqn.solve();
                p = ph_rgh + rho*gh + pRef;
                thermo.correct();
                rho = thermo.rho();

                Info<< "Reference hydrostatic pressure variation "
                    << (max(ph_rgh) - min(ph_rgh)).value() << endl;
            }

            ph_rgh.write();

            thermo.pRef() = p - rho*gh;
            thermo.rhoRef() = rho;
            T = T_orig;
            he = he_orig;
        }
        else
        {
            Info << "Restart condition of atmHydrostaticInitialisation.." << endl;
            thermo.correct();
            rho = thermo.rho();
        }
    }
    else
    {
        // Force p_rgh to be consistent with p
        // p_rgh = thermo.p() - thermo.rho()*gh - pRef;
    }
}


// ************************************************************************* //
