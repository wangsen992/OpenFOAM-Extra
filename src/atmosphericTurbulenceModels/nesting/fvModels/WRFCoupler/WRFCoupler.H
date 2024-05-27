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
    Foam::fv::WRFCoupler

Description
    WRF coupler to be used for the atm*Foam solvers loading input 

Usage
    Example usage:
    \verbatim
    buoyancyForce1
    {
        type        WRFCoupler;

        U           U;            // Name of the velocity field
        Ug          (10 0 0);     // Value of geostrophic velocity
        f           1e-4;         // Value of Coriolis Parameter
    }
    \endverbatim

SourceFiles
    WRFCoupler.C

\*---------------------------------------------------------------------------*/

#ifndef WRFCoupler_H
#define WRFCoupler_H

#include "HashTable.H"
#include "fvModel.H"
#include "uniformDimensionedFields.H"
#include "dimensionedVector.H"
#include "fluidAtmThermo.H"
#include "WRF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class buoyancyForce Declaration
\*---------------------------------------------------------------------------*/
template<typename Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> smooth(GeometricField<Type, fvPatchField, volMesh>& psi, label smooth_n=5);
 
class WRFCoupler
:
    public fvModel
{
    // Private Data

        //- Host case fvMesh
        const fvMesh& mesh_;
        const fluidAtmThermo& thermo_;

        //- WRF driver
        WRF& wrf_;

        //- CellList for applying nesting
        DynamicList<label, 20> nestingCells_;
        HashTable<scalar, label> nestingCellTbl_;
        pointField nestingCellCentres_;
        volScalarField cellWeights_;
        scalar nestingDist_;
        scalar relaxationFactor_;

        //- Flags to control update
        label currTimeInd_;

        //- Optional phase name
        word phaseName_;

    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- static function for smoothing

public:

    //- Runtime type information
    TypeName("WRFCoupler");

    // Constructors

        //- Construct from explicit source name and mesh
        WRFCoupler
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        WRFCoupler(const WRFCoupler&) = delete;

    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluate

            //- Add explicit contribution to incompressible momentum equation
            template<typename T>
            void addSup
            (
                fvMatrix<T>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to compressible momentum equation
            template<typename T>
            void addSup
            (
                const volScalarField& rho,
                fvMatrix<T>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to phase momentum equation
            template<typename T>
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<T>& eqn,
                const word& fieldName
            ) const;
        
        // Correct fvModel
        void correct();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const WRFCoupler&) = delete;
        
        void addSup
        (
            const volScalarField& alpha,
            fvMatrix<Foam::scalar>& eqn,
            const word& fieldName
        ) const;

        void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<Foam::vector>& eqn,
            const word& fieldName
        ) const;
        void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<Foam::scalar>& eqn,
            const word& fieldName
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<typename T>
void Foam::fv::WRFCoupler::addSup
(
    fvMatrix<T>& eqn,
    const word& fieldName
) const
{
}


template<typename T>
void Foam::fv::WRFCoupler::addSup
(
    const volScalarField& rho,
    fvMatrix<T>& eqn,
    const word& fieldName
) const
{
  typedef GeometricField<T, fvPatchField, volMesh> psiType;
  auto psi = wrf_.time().lookupObjectRef<psiType>(fieldName);
  auto psi_foam = mesh().time().lookupObjectRef<psiType>(fieldName);

  const scalarField& V = mesh().V();
  auto source = eqn.source();
  for(auto i : nestingCells_)
  {
    source[i] -= V[i] * rho[i] * relaxationFactor_*(psi[i] - psi_foam[i]);
  }
}

template<typename Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> smooth(GeometricField<Type, fvPatchField, volMesh>& psi, label smooth_n)
{
    typedef GeometricField<Type, fvPatchField, volMesh> psiType;
    psiType f0(psi);
    psiType f(psi);
    label smoothingPasses = smooth_n;
    for (label pass=1; pass<=smoothingPasses; pass++)
    {
        f = fvc::average(fvc::interpolate(f0));
        f.correctBoundaryConditions();
        f0 = f;
    }
    return tmp<psiType> ( new psiType(f0));
}


// template<typename T>
// void Foam::fv::WRFCoupler::addSup
// (
//     const volScalarField& alpha,
//     const volScalarField& rho,
//     fvMatrix<T>& eqn,
//     const word& fieldName
// ) const
// {
//   Info << "[fvModel] adding wrf field " << fieldName << endl;
//   typedef GeometricField<T, fvPatchField, volMesh> psiType;
//   auto psi = wrf_.time().lookupObjectRef<psiType>(fieldName);
//   auto psi_foam = mesh().time().lookupObjectRef<psiType>(fieldName);
// 
//   auto source = eqn.source();
//   const scalarField& V = mesh().V();
//   for(auto i : nestingCells_)
//   {
//     source[i] += V[i] * alpha[i] * rho[i] * relaxationFactor_*(psi[i] - psi_foam[i]);
//   }
// }

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
