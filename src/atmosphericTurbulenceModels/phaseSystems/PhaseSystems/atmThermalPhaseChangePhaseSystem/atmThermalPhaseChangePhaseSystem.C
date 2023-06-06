/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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

#include "atmThermalPhaseChangePhaseSystem.H"
#include "alphatPhaseChangeWallFunctionFvPatchScalarField.H"
#include "fvcVolumeIntegrate.H"
#include "fvmSup.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::addDmdts
(
    PtrList<volScalarField>& dmdts
) const
{
    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs_, dmdtfIter)
    {
        const phasePair& pair = this->phasePairs_[dmdtfIter.key()];
        const volScalarField& dmdtf = *dmdtfIter();

        addField(pair.phase1(), "dmdt", dmdtf, dmdts);
        addField(pair.phase2(), "dmdt", - dmdtf, dmdts);
    }

    forAllConstIter(phaseSystem::dmdtfTable, nDmdtfs_, nDmdtfIter)
    {
        const phasePair& pair = this->phasePairs_[nDmdtfIter.key()];
        const volScalarField& nDmdtf = *nDmdtfIter();

        addField(pair.phase1(), "dmdt", nDmdtf, dmdts);
        addField(pair.phase2(), "dmdt", - nDmdtf, dmdts);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::
atmThermalPhaseChangePhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),
    volatile_(this->template lookupOrDefault<word>("volatile", "none")),
    condensePhase_(this->template lookupOrDefault<word>("condensePhase", "air")),
    dmdt0s_(this->phases().size())
{
    this->generatePairsAndSubModels
    (
        "saturation",
        saturationModels_
    );

    // Check that models have been specified in the correct combinations
    forAllConstIter
    (
        saturationModelTable,
        saturationModels_,
        saturationModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[saturationModelIter.key()];

        this->template validateMassTransfer<saturationModel>(pair);

        if
        (
            !this->heatTransferModels_.found(pair)
         || !this->heatTransferModels_[pair].first().valid()
         || !this->heatTransferModels_[pair].second().valid()
        )
        {
             FatalErrorInFunction
                 << "A heat transfer model for both sides of the " << pair
                 << "pair is not specified. This is required by the "
                 << "corresponding saturation model"
                 << exit(FatalError);
        }
    }

    // Generate interfacial mass transfer fields, initially assumed to be zero
    forAllConstIter
    (
        saturationModelTable,
        saturationModels_,
        saturationModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[saturationModelIter.key()];

        dmdtfs_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "thermalPhaseChange:dmdtf",
                        pair.name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );

        d2mdtdpfs_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "thermalPhaseChange:d2mdtdpf",
                        pair.name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar((dimDensity/dimTime)/dimPressure, 0)
            )
        );

        Tfs_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "thermalPhaseChange:Tf",
                        pair.name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                (pair.phase1().thermo().T() + pair.phase2().thermo().T())/2
            )
        );

        nDmdtfs_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        IOobject::groupName
                        (
                          "thermalPhaseChange:nucleation:dmdtf",
                          "H2O"
                        ),
                        pair.name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::
~atmThermalPhaseChangePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
const Foam::saturationModel&
Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::saturation
(
    const phasePairKey& key
) const
{
    return saturationModels_[key];
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::dmdtf
(
    const phasePairKey& key
) const
{
    tmp<volScalarField> tDmdtf = BasePhaseSystem::dmdtf(key);

    const scalar dmdtfSign =
        Pair<word>::compare(this->phasePairs_[key], key);

    if (dmdtfs_.found(key))
    {
        tDmdtf.ref() += dmdtfSign**dmdtfs_[key];
    }

    if (nDmdtfs_.found(key))
    {
        tDmdtf.ref() += dmdtfSign**nDmdtfs_[key];
    }

    return tDmdtf;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    addDmdts(dmdts);

    return dmdts;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::d2mdtdps() const
{
    PtrList<volScalarField> d2mdtdps(BasePhaseSystem::d2mdtdps());

    forAllConstIter(phaseSystem::dmdtfTable, d2mdtdpfs_, d2mdtdpfIter)
    {
        const phasePair& pair = this->phasePairs_[d2mdtdpfIter.key()];
        const volScalarField& d2mdtdpf = *d2mdtdpfIter();

        addField(pair.phase1(), "d2mdtdp", d2mdtdpf, d2mdtdps);
        addField(pair.phase2(), "d2mdtdp", - d2mdtdpf, d2mdtdps);
    }

    return d2mdtdps;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);
    this->addDmdtUfs(nDmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);
    this->addDmdtUfs(nDmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    // Create temperatures at which to evaluate nucleation mass transfers
    phaseSystem::dmdtfTable Tns;
    forAllConstIter(phaseSystem::dmdtfTable, nDmdtfs_, nDmdtfIter)
    {
        const phasePair& pair = this->phasePairs_[nDmdtfIter.key()];
        const saturationModel& satModel = this->saturation(nDmdtfIter.key());

        // const volScalarField Tsat(saturationModelIter()->Tsat(thermo1.p()));
        const volScalarField psat(satModel.pSat(pair.phase1().thermo().T()));
        // [Note] This is an ad hoc implementation. Should take R directly
        // from thermo info
        tmp<volScalarField> Tsat
        (
          pair.phase1().thermo().T()
        );

        Tns.insert(pair, Tsat.ptr());
    }

    // Mass transfer terms
    if (volatile_ != "none")
    {
        {
            phaseSystem::dmidtfTable dmidtfs;

            forAllConstIter(phaseSystem::dmdtfTable, dmdtfs_, dmdtfIter)
            {
                const phasePair& pair = this->phasePairs_[dmdtfIter.key()];
                const volScalarField& dmdtf = *dmdtfIter();

                dmidtfs.insert(pair, new HashPtrTable<volScalarField>());
                dmidtfs[pair]->insert(volatile_, new volScalarField(dmdtf));
            }

            this->addDmidtHefs
            (
                dmidtfs,
                Tfs_,
                BasePhaseSystem::latentHeatScheme::upwind,
                BasePhaseSystem::latentHeatTransfer::heat,
                eqns
            );
        }
        {
            phaseSystem::dmidtfTable nDmidtfs;

            forAllConstIter(phaseSystem::dmdtfTable, nDmdtfs_, nDmdtfIter)
            {
                const phasePair& pair = this->phasePairs_[nDmdtfIter.key()];
                const volScalarField& nDmdtf = *nDmdtfIter();

                nDmidtfs.insert(pair, new HashPtrTable<volScalarField>());
                nDmidtfs[pair]->insert(volatile_, new volScalarField(nDmdtf));
            }

            this->addDmidtHefs
            (
                nDmidtfs,
                Tns,
                BasePhaseSystem::latentHeatScheme::upwind,
                BasePhaseSystem::latentHeatTransfer::heat,
                eqns
            );
        }
    }
    else
    {
        this->addDmdtHefs
        (
            dmdtfs_,
            Tfs_,
            BasePhaseSystem::latentHeatScheme::upwind,
            BasePhaseSystem::latentHeatTransfer::heat,
            eqns
        );
        this->addDmdtHefs
        (
            nDmdtfs_,
            Tns,
            0,
            BasePhaseSystem::latentHeatScheme::upwind,
            eqns
        );
    }

    // Lagging
    {
        PtrList<volScalarField> dmdts(this->phases().size());

        addDmdts(dmdts);

        forAllConstIter(phaseSystem::phaseModelList, this->phases(), phaseIter)
        {
            const phaseModel& phase = phaseIter();

            if (dmdt0s_.set(phase.index()))
            {
                *eqns[phase.name()] +=
                    fvm::Sp
                    (
                        dmdt0s_[phase.index()] - dmdts[phase.index()],
                        phase.thermo().he()
                    );
            }
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::specieTransfer() const
{
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    phaseSystem::specieTransferTable& eqns = eqnsPtr();

    if (volatile_ != "none")
    {
        {
            phaseSystem::dmidtfTable dmidtfs;

            forAllConstIter(phaseSystem::dmdtfTable, dmdtfs_, dmdtfIter)
            {
                const phasePair& pair = this->phasePairs_[dmdtfIter.key()];
                const volScalarField& dmdtf = *dmdtfIter();

                dmidtfs.insert(pair, new HashPtrTable<volScalarField>());
                dmidtfs[pair]->insert(volatile_, new volScalarField(dmdtf));
            }

            this->addDmidtYf(dmidtfs, eqns);
        }

        {
            phaseSystem::dmidtfTable nDmidtfs;

            forAllConstIter(phaseSystem::dmdtfTable, nDmdtfs_, nDmdtfIter)
            {
                const phasePair& pair = this->phasePairs_[nDmdtfIter.key()];
                const volScalarField& nDmdtf = *nDmdtfIter();

                nDmidtfs.insert(pair, new HashPtrTable<volScalarField>());
                nDmidtfs[pair]->insert(volatile_, new volScalarField(nDmdtf));
            }

            this->addDmidtYf(nDmidtfs, eqns);
        }
    }
    else
    {
        this->addDmdtYfs(dmdtfs_, eqns);
        this->addDmdtYfs(nDmdtfs_, eqns);
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::
correctContinuityError()
{
    dmdt0s_ = PtrList<volScalarField>(this->phases().size());

    addDmdts(dmdt0s_);

    BasePhaseSystem::correctContinuityError();
}


template<class BasePhaseSystem>
void
Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::correctInterfaceThermo()
{
    Info << "[ThermalPhaseChange] Entering correctInterfaceThermo() " << endl;

    forAllConstIter
    (
        saturationModelTable,
        saturationModels_,
        saturationModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[saturationModelIter.key()];

        // Compute Nucleation Mass Transfer
        const phaseModel& phase = pair.phase1().name() == condensePhase_ ? pair.phase1() : pair.phase2();
        const phaseModel& otherPhase = pair.phase1().name() != condensePhase_ ? pair.phase1() : pair.phase2();
        const rhoReactionThermo& thermo
        (
          phase.mesh().lookupObjectRef<rhoReactionThermo>
          (
              IOobject::groupName(basicThermo::dictName, condensePhase_)
          )
        );
        const volScalarField& T(thermo.T());
        volScalarField qv(volatile_, phase.Y(volatile_));
        qv += VSMALL;
        const volScalarField es(saturationModelIter()->pSat(T));
        const volScalarField qs
        (
            thermo.composition().Wi(thermo.composition().index(qv)) 
            / (thermo.W() * thermo.p())  * es 
            / dimensionedScalar(dimensionSet(-1, 0, 0, 0, 1), 1.0)
        );
        Info << "phase name : " << phase.name() << endl;
        Info << "max(qs) = " << gMax(qs) << endl;
        Info << "Unit of qs: " << qs.dimensions() << endl;

        volScalarField& ndmdtf(*this->nDmdtfs_[pair]);
        volScalarField& dmdtf(*this->dmdtfs_[pair]);

        dimensionedScalar dt = phase.mesh().time().deltaT();
        volScalarField dq((qv - qs));
        volScalarField rhodt(thermo.rho() / dt);
        Info << "[atmThermalPhaseChangePhaseSystem] avg(qv) = " << average(qv) << endl;
        Info << "[atmThermalPhaseChangePhaseSystem] avg(qs) = " << average(qs) << endl;
        Info << "[atmThermalPhaseChangePhaseSystem] gMax(dq) = " << gMax(dq) << endl;
        Info << "[atmThermalPhaseChangePhaseSystem] gMin(dq) = " << gMin(dq) << endl;

        // Since ndmdtf is added to phase1 and subtracted from phase2, 
        // for nucleation, and condensePhase_, "air", should be subtracted, 
        // so if phase1 is condensePhase_, pos(qv-qs) * -1 = ndmdtfNew, else, 
        // pos(qv-qs) = ndmdtfNew
        // When qv - qs < 0, air is undersaturated, therefore we consider
        // only neg(qv - qs), we here assume that evaporation is only 5% 
        // of the deficiency in saturation
        
        /* old version 
        volScalarField ndmdtfNew((qv - qs) * thermo.rho() / dt);
        ndmdtfNew = dm * (pair.phase1().name() == condensePhase_ ? -1 : 1)
        ndmdtfNew *= (pair.phase1().name() == condensePhase_ ? -1 : 1);
        ndmdtfNew *= pos(ndmdtfNew);
        */

        const scalar sign = pair.phase1().name() == condensePhase_ ? -1 : 1;

        volScalarField ndmdtfNew = pos(dq) * dq * sign * rhodt; // nucleation of droplets
        // Evaporation is a function of mass diffusivity, surface area
        // (diameter), also, the mass transfer rate should also be limited 
        // by the total water liquid content. A dynamic way to limit it is
        // through total surface area, which can be done through the
        // interfaceComposition model
        
        // evaporation of liquid phase
        // 0.01 is hard coded diffusivity and the rate is limited by phase
        // fraction
        // Changed to 0.5 to test
        volScalarField dmdtfNew
        (
            neg(dq) 
          * dq * otherPhase 
          * 0.5 * rhodt * sign
        );

        // mass transfer update
        {

            // [Note] This is an ad hoc implementation. Should take R directly
            // from thermo info

            const scalar dmdtfRelax = 1.0;
                // this->mesh().fieldRelaxationFactor(dmdtf.member());

            ndmdtf = (1 - dmdtfRelax)*ndmdtf + dmdtfRelax*ndmdtfNew;
            // ndmdtf *= 0;
            dmdtf =  (1 - dmdtfRelax) * dmdtf + dmdtfRelax * dmdtfNew;

            Info<< ndmdtf.name()
                << ": min = " << gMin(ndmdtf.primitiveField())
                << ", mean = " << gAverage(ndmdtf.primitiveField())
                << ", max = " << gMax(ndmdtf.primitiveField())
                << ", integral = " << fvc::domainIntegrate(ndmdtf).value()
                << endl;
            Info<< dmdtf.name()
                << ": min = " << gMin(dmdtf.primitiveField())
                << ", mean = " << gAverage(dmdtf.primitiveField())
                << ", max = " << gMax(dmdtf.primitiveField())
                << ", integral = " << fvc::domainIntegrate(dmdtf).value()
                << endl;
        }

    }
}


template<class BasePhaseSystem>
bool Foam::atmThermalPhaseChangePhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
