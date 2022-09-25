/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
#include "autoPtr.H"
#include "dimensionSet.H"
#include "thetaRhoAtmThermo.H"
#include "speciesTable.H"
#include "thermodynamicConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
template<class BasicRhoThermo, class MixtureType>
void Foam::thetaRhoAtmThermo<BasicRhoThermo, MixtureType>::validate_mixture()
{
    wordList basicSpecie{"dryAir", "H2O"};
    Info << "Validating mixture species elements" << endl;
    const speciesTable& species(this->MixtureType::species());
    Info << species << endl;
    forAll(basicSpecie, i)
    {
      Info << basicSpecie[i] << endl;
      if (! species.found(basicSpecie[i]))
      {
          Foam::FatalError << "Basic species requirement" 
                           << basicSpecie 
                           << " not met." << endl
                           << basicSpecie[i]
                           << " is missing."
                           << Foam::exit(Foam::FatalError);
      }
    }
}


template<class BasicRhoThermo, class MixtureType>
void Foam::thetaRhoAtmThermo<BasicRhoThermo, MixtureType>::calculate()
{

    Info << "Running calculate() ..." << endl;
    // [To-Do] This can be updated for more general method, including all
    // the gamma below
    // [To-Do] Gamma, Cp, Cv and T are iterative. Solve this
    scalar gamma = 0.2854; // Currently using gamma value of dry air
    scalar p0 = this->p0_.value();
    // Initiate Initial value loading (those two should be updated) 
    const scalarField& pCells = this->p_;
    const scalarField& thetaCells = this->theta_;

    scalarField& hCells = this->he();
    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    // Info << "Data loading completed" << endl;
    // Update internal cell values
    forAll(TCells, celli)
    {
        // Info << "Working on " << celli <<  endl;
        const typename MixtureType::thermoMixtureType& thermoMixture =
            this->cellThermoMixture(celli);

        const typename MixtureType::transportMixtureType& transportMixture =
            this->cellTransportMixture(celli, thermoMixture);


        TCells[celli] = thetaCells[celli] 
                  * this->BasicRhoThermo::exner(pCells[celli], p0, gamma);
        hCells[celli] = thermoMixture.HE
        (
            pCells[celli],
            TCells[celli]
        );

        CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
        CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
        psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = thermoMixture.rho(pCells[celli], TCells[celli]);

        muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] =
            transportMixture.kappa(pCells[celli], TCells[celli])
           /thermoMixture.Cp(pCells[celli], TCells[celli]);
        // Info << "Update on " << celli << " completed" << endl;
    }
    // Info << "Internal Cell update completed" << endl;

    // Init reference of boundaryField
    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& thetaBf =
        this->theta_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& CpBf =
        this->Cp_.boundaryFieldRef();

    volScalarField::Boundary& CvBf =
        this->Cv_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    // Update boundaryField values
    forAll(this->theta_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& ptheta = thetaBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        // Update if T is fixes boundary
        if (ptheta.fixesValue())
        {
            forAll(ptheta, facei)
            {
                const typename MixtureType::thermoMixtureType& thermoMixture =
                    this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                pT[facei] = ptheta[facei] 
                  / this->BasicRhoThermo::exner(pp[facei], p0, gamma);
                phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);

                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                prho[facei] = thermoMixture.rho(pp[facei], pT[facei]);

                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                palpha[facei] =
                    transportMixture.kappa(pp[facei], pT[facei])
                   /thermoMixture.Cp(pp[facei], pT[facei]);
            }
        }

        // Update if theta is not fixing value
        else
        {
            forAll(ptheta, facei)
            {
                const typename MixtureType::thermoMixtureType& thermoMixture =
                    this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                // gamma = pCp[facei] / pCv[facei];
                gamma = 0.2854;
                pT[facei] = ptheta[facei] * this->BasicRhoThermo::exner(pp[facei], p0, gamma);
                phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);

                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                prho[facei] = thermoMixture.rho(pp[facei], pT[facei]);

                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                palpha[facei] =
                    transportMixture.kappa(pp[facei], pT[facei])
                   /thermoMixture.Cp(pp[facei], pT[facei]);
            }
        }
    }
    Info << "Boundary cell update completed." << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
Foam::thetaRhoAtmThermo<BasicRhoThermo, MixtureType>::thetaRhoAtmThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicRhoThermo, MixtureType>(mesh, phaseName)
{
    validate_mixture();
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
Foam::thetaRhoAtmThermo<BasicRhoThermo, MixtureType>::~thetaRhoAtmThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
void Foam::thetaRhoAtmThermo<BasicRhoThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    Info << "Running correct in thetaRhoAtmThermo.." << endl;

    // Now run calculate
    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


// ************************************************************************* //
