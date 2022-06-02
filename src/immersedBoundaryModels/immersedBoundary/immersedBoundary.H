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
    Foam::fv::immersedBoundary

Description
    Base class for immersedBoundary. Offers interface for source term 
    calculations and surface analysis functions

SourceFiles
    immersedBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundary_H
#define immersedBoundary_H

#include "cellSetHashTables.H"
#include "volFieldsFwd.H"
#include "cellSet.H"
#include "triSurfaceMesh.H"
#include "triSurfaceSearch.H"
#include "polyMesh.H"
#include "meshSearch.H"

namespace Foam
{

class immersedBoundary
{
    // Private Member Data
    
    // References to mesh and surface
    const polyMesh&         mesh_;
    const triSurfaceMesh&   surf_;
    
    // geometries relevant to the IBM. 
    scalar              dpp_;
    cellSet             forcingCells_;
    vectorCellSet       surfaceFaceCenters_;
    vectorCellSet       surfaceFaceNormals_;
    vectorCellSet       forcingPoints_;
    vectorCellSet       forcingPointAreas_;
    vectorCellSet       forcingPointNormals_;
    vectorCellSet       probePoints_;

    // parameters for computing source terms (elasticity formulation)
    dimensionedScalar   Kp_;
    dimensionedScalar   KI_;
    dimensionedScalar   ustar_;

    // Storage of source terms
    dimensionedVectorCellSet     fU_;

protected:
    
    // Protected Member Functions

    // Compute the forcing point (center of gravity of surface within the 
    // forcing cell) and forcing normal.
    static void calcForcingPointAndNormal
    (
      const polyMesh& mesh, 
      const triSurface& surf, 
      const cellSet& forcingCells,
      vectorCellSet& forcingPoints,
      vectorCellSet& forcingPointNormals
    );

    // compute all attributes
    void init();

    // Read coefficients
    void readCoeffs();

public:
    
    // Constructors
      
      //- Construct from surface and mesh
      immersedBoundary
      (
          const polyMesh& mesh, 
          const triSurfaceMesh& surf
      );
    // Destructor
    virtual ~immersedBoundary(){};

    // Member Functions
    
      //- Access

        const polyMesh& mesh(){return mesh_;};

        const triSurfaceMesh& surfaceMesh(){return surf_;};

        const cellSet& forcingCells(){return forcingCells_;};

        const vectorCellSet& surfaceFaceCenters(){return surfaceFaceCenters_;};

        const vectorCellSet& surfaceFaceNormals(){return surfaceFaceNormals_;};

        const vectorCellSet& forcingPoints(){return forcingPoints_;};

        const vectorCellSet& forcingPointNormals(){return forcingPointNormals_;};

        const vectorCellSet& probePoints(){return probePoints_;};

    //- Computing source terms
    
        virtual tmp<dimensionedVectorCellSet> fFp(volVectorField& U);
        virtual tmp<dimensionedVectorCellSet> fnFp(volVectorField& U);
        dimensionedVectorCellSet& fU();
        virtual scalarField kSource();
        virtual scalarField epsilonSource();
        virtual tensorField RSource();

    //- interface to call correcting everything
    void correct();

};
}

#endif
