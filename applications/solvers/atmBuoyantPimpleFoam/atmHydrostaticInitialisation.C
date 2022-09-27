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

#include "atmHydrostaticInitialisation.H"

#include "fluidAtmThermo.H"
#include "fvmLaplacian.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "constrainPressure.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::atmHydrostaticInitialisation
(
    volScalarField& p_rgh,
    volScalarField& rho,
    volScalarField& theta,
    const volVectorField& U,
    const volScalarField& gh,
    const surfaceScalarField& ghf,
    const uniformDimensionedScalarField& pRef,
    fluidAtmThermo& thermo,
    const dictionary& dict
)
{
    if (dict.lookupOrDefault<bool>("atmHydrostaticInitialisation", false))
    {
        const fvMesh& mesh = p_rgh.mesh();

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
            volScalarField& T = thermo.T();
            volScalarField& he = thermo.he();
            // theta can be substituted with other potential temperature
            // variant. To differentia different theta types and allow
            // automatic conversion in the poisson equation, name of theta can
            // be used as an identifier to conditional statements. 
            volScalarField theta_orig = theta;

            thermo.correct();
            rho = thermo.rho();
            p = ph_rgh + rho*gh + pRef;
            T = theta_orig * thermo.exner(p, thermo.p0(), 0.2854);
            he = thermo.he(p, T);

            // Report range of states 
            Info << "ph_rgh: " << max(ph_rgh)<< " " <<  min(ph_rgh) << endl;
            Info << "rho: " << max(rho) << " " << min(rho) << endl;
            Info << "gh: " << max(gh) << " " << min(gh)<< endl;

            Info << "p: " << max(p) << " " << min(p) << endl;
            Info << "p0: " << thermo.p0() << endl;
            Info << "T: " <<  max(T) << " " << min(T) << endl;
            Info << "theta: " <<  max(theta) << " " << min(theta) << endl;
            Info << "theta_orig: " <<  max(theta_orig) << " " << min(theta_orig) << endl;

            // thermo.correct();
            // rho = thermo.rho();

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
                T = theta_orig * thermo.exner(p, thermo.p0(), 0.2854);
                he = thermo.he(p, T);
                thermo.correct();
                rho = thermo.rho();
                Info << "ph_rgh: " << max(ph_rgh)<< " " <<  min(ph_rgh) << endl;
                Info << "rho: " << max(rho) << " " << min(rho) << endl;
                Info << "gh: " << max(gh) << " " << min(gh)<< endl;
                Info << "p: " << max(p) << " " << min(p) << endl;
                Info << "T: " <<  max(T) << " " << min(T) << endl;
                Info << "theta: " <<  max(thermo.theta()) << " " << min(thermo.theta()) << endl;
                Info << "theta_orig: " <<  max(theta_orig) << " " << min(theta_orig) << endl;

                Info<< "Hydrostatic pressure variation "
                    << (max(ph_rgh) - min(ph_rgh)).value() << endl;
            }

            ph_rgh.write();

            p_rgh = ph_rgh;
        }
        else
        {
            Info << "Non-Restart condition of atmHydrostaticInitialisation.." << endl;
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
