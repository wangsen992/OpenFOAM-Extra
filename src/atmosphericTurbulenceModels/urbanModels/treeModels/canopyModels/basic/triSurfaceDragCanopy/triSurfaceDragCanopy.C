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

Class
    Foam::triSurfaceDragCanopy

Description
    A test implementation that does not include a drag model yet

SourceFiles
    triSurfaceDragCanopy.C

\*---------------------------------------------------------------------------*/

#include "triSurfaceDragCanopy.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicDragCanopy>
void Foam::triSurfaceDragCanopy<BasicDragCanopy>::calculate()
{
    if (debug)
    {
        Info << "Running calculate..." << endl;
    }

    const fvMesh& mesh_ (this->mesh());

    // Get velocity field
    const volVectorField& UCells 
      = mesh_.lookupObjectRef<volVectorField>("U");
    
    if (debug)
    {
        Info << "Size of U: " << UCells.size() << endl;
    }

    // Prepare for assigning values
    labelList cells = this->canopyCells().sortedToc();
    scalarCellSet& CdCells = this->Cd_;
    dimensionedScalarCellSet& ladCells = this->lad_;

    // Assign fU values based on simple drag model
    forAll(cells, i)
    {
        this->fU_.insert
        (
          cells[i],
          CdCells[cells[i]] * mag(UCells[cells[i]]) * ladCells[cells[i]] * UCells[cells[i]]
        );
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicDragCanopy>
Foam::triSurfaceDragCanopy<BasicDragCanopy>::triSurfaceDragCanopy
(
    const fvMesh& mesh
)
:
    triSurfaceCanopy<BasicDragCanopy>(mesh)
{
    this->fU_.resize(this->canopyCells().toc().size());
    calculate();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicDragCanopy>
void Foam::triSurfaceDragCanopy<BasicDragCanopy>::correct()
{
    
    if (debug)
    {
        Info << "Running correct..." << endl;
    }

    calculate();

}
