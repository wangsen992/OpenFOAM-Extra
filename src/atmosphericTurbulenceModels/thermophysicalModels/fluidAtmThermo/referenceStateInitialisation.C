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
#include "meshSearch.H"

#include "specie.H"
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
                    "ph_rgh0",
                    "0",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );

        //if (!mesh.time().restart())
        //{
            volScalarField& p = thermo.p();
            p = ph_rgh + rho*gh + pRef;

            // Enforce temperature in thermo to compute hydrostatic reference
            // state, which is T = T_0 - g/cp * z
            // For the perfect gas case, cp is constant and set in dict
            volScalarField T_orig = thermo.T();
            volScalarField he_orig = thermo.he();
            volScalarField p_orig = thermo.p();

            // Use average temeprature, Tb
            dimensionedScalar Tb = dimensionedScalar(dimTemperature, 290);
            volScalarField T_neutral(Tb + gh / thermo.Cp());

            // Enforce empty species concentration for reference state
            // calculation
            PtrList<volScalarField> Y_orig
            (
              thermo.composition().species().size()
            );

            forAll(thermo.composition().species(), speciei)
            { 
                
                Y_orig.set
                (
                  speciei,
                  new volScalarField
                  (
                    thermo.composition().Y()[speciei]
                  )
                );

                if (thermo.composition().defaultSpecie() != speciei)
                {
                    Info << "Set non default specie " << thermo.composition().Y()[speciei].name() << "to 0." << endl;
                    thermo.composition().Y()[speciei] = 0;
                }
                else
                {
                    Info << "Set default specie " << thermo.composition().Y()[speciei].name() << "to 1." << endl;
                    thermo.composition().Y()[speciei] = 1;
                }
                    
            }

            thermo.correct();

            //- Enforce temperature and energy profile
            // Note energy must be enforced as temperature is 
            // derived from temperature. 
            volScalarField& T = thermo.T();
            volScalarField& he = thermo.he();
            T = T_neutral;
            he = thermo.he(p, T_neutral);

            // Update steps
            // After thermo correction, get the updated density
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

                // [Temporary fix]
                ph_rgh = ph_rgh - max(ph_rgh);
                fvScalarMatrix ph_rghEqn
                (
                    fvm::laplacian(rhof, ph_rgh) == fvc::div(phig)
                );

                ph_rghEqn.solve();

                // Compute pressure using updated values of p_rgh, hydrostatic
                // pressure and reference base pressure value
                p = ph_rgh + rho*gh + pRef;

                // Energy value is updated with the updated pressure 
                he = thermo.he(p, T_neutral);
                thermo.correct();
                rho = thermo.rho();

                Info<< "Reference hydrostatic pressure variation"
                    << (max(ph_rgh) - min(ph_rgh)).value() << endl;
            }

            thermo.pRef() = p;
            thermo.rhoRef() = rho;

            // Reset the temperature and energy field as before
            T = T_orig;
            he = he_orig;
            thermo.p() = p_orig;

            forAll(thermo.composition().species(), speciei)
            { 
              thermo.composition().Y()[speciei] = Y_orig[speciei];
            }
            // thermo.correct();

        // }
        // else
        // {
        //     Info << "Restart condition of atmHydrostaticInitialisation.." << endl;
        //     thermo.correct();
        //     rho = thermo.rho();
        // }
    }
    else
    {
        // Force p_rgh to be consistent with p
        // p_rgh = thermo.p() - thermo.rho()*gh - pRef;
    }
}


// ************************************************************************* //
