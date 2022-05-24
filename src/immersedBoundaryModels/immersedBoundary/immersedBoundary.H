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
    calculations.

SourceFiles
    immersedBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundary_H
#define immersedBoundary_H

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
    
    const polyMesh&         mesh_;
    const triSurfaceMesh&   surf_;
    PtrList<
    
    cellSet           forcingCells_;
    pointField        surfaceFaceCenters_;
    vectorField       surfaceFaceNormals_;
    pointField        forcingPoints_;
    vectorField       forcingPointNormals_;
    pointField        probePoints_;

protected:
    
    // Protected Member Functions
    static labelHashSet calcForcingCells
    (
      const polyMesh& mesh,
      const triSurfaceMesh& surf
     );

    void calcForcingPointAndNormal
    (
      const polyMesh& mesh, 
      const triSurface& surf, 
      const labelHashSet& forcingCells,
      pointField& forcingPoints,
      vectorField& forcingPointNormals
    );

    // compute all attributes
    void calc();

public:
    
    // Constructors
      
      //- Construct from surface and mesh
      immersedBoundary
      (
          const polyMesh& mesh, 
          const triSurfaceMesh& surf
      );

    // Member Functions
    
      //- Access

        const polyMesh& mesh(){return mesh_;};

        const triSurfaceMesh& surfaceMesh(){return surf_;};

        const cellSet& forcingCells(){return forcingCells_;};

        const pointField& surfaceFaceCenters(){return surfaceFaceCenters_;};

        const vectorField& surfaceFaceNormals(){return surfaceFaceNormals_;};

        const pointField& forcingPoints(){return forcingPoints_;};

        const vectorField& forcingPointNormals(){return forcingPointNormals_;};

        const pointField probePoints(){return probePoints_;};

    //- Computing source terms
    
        tmp<volVectorField> fFp(volVectorField& U);
        tmp<volVectorField> fnFp(volVectorField& U);
        tmp<volVectorField> fpFp();
        tmp<volScalarField> kSource();
        tmp<volScalarField> epsilonSource();
        tmp<volTensorField> RSource();

};
}

#endif
