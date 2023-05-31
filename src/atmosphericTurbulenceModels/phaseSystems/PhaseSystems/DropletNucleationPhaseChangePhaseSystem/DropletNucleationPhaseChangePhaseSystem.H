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

Class
    Foam::DropletNucleationPhaseChangePhaseSystem

Description
    Adapted from thermalPhaseChange. Focusing on nucleation only.


SourceFiles
    DropletNucleationPhaseChangePhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef DropletNucleationPhaseChangePhaseSystem_H
#define DropletNucleationPhaseChangePhaseSystem_H

#include "phaseSystem.H"
#include "saturationModel.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class DropletNucleationPhaseChangePhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class DropletNucleationPhaseChangePhaseSystem
:
    public BasePhaseSystem
{
    // Private typedefs

        typedef HashTable
        <
            autoPtr<saturationModel>,
            phasePairKey,
            phasePairKey::hash
        > saturationModelTable;

        using latentHeatScheme = typename BasePhaseSystem::latentHeatScheme;

        using latentHeatTransfer = typename BasePhaseSystem::latentHeatTransfer;


    // Private data

        //- Saturation models used to evaluate Tsat = Tf
        saturationModelTable saturationModels_;

        //- Mass transfer rates
        phaseSystem::dmdtfTable dmdtfs_;

        //- Mass transfer linearisation coeffs
        phaseSystem::dmdtfTable d2mdtdpfs_;

        //- Interface temperatures
        phaseSystem::dmdtfTable Tfs_;

        //- Nucleate Mass transfer rates
        phaseSystem::dmdtfTable nDmdtfs_;

        //- Nucleate thermal energy transfer rates
        phaseSystem::dmdtfTable nDmdtLfs_;

        //- Previous continuity error update phase dmdts for the heat transfer
        //  function
        PtrList<volScalarField> dmdt0s_;

    // Private Member Functions

        //- Sum the mass transfer rates for each phase into the given list
        void addDmdts(PtrList<volScalarField>&) const;


public:

    // Constructors

        //- Construct from fvMesh
        DropletNucleationPhaseChangePhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~DropletNucleationPhaseChangePhaseSystem();


    // Member Functions

        //- Return the saturationModel
        const saturationModel& saturation(const phasePairKey& key) const;

        //- Return the mass transfer rate for an interface
        virtual tmp<volScalarField> dmdtf(const phasePairKey& key) const;

        //- Return the mass transfer rates for each phase
        virtual PtrList<volScalarField> dmdts() const;

        //- Return the mass transfer linearisation coeffs for each phase
        virtual PtrList<volScalarField> d2mdtdps() const;

        //- Return the momentum transfer matrices for the cell-based algorithm
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransfer();

        //- Return the momentum transfer matrices for the face-based algorithm
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransferf();

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;

        //- Return the specie transfer matrices
        virtual autoPtr<phaseSystem::specieTransferTable>
            specieTransfer() const;

        //- Store phase dmdts at the during the continuity error update
        virtual void correctContinuityError();

        //- Correct the interface thermodynamics
        virtual void correctInterfaceThermo();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DropletNucleationPhaseChangePhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
