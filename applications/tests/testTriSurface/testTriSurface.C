/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

Application
    testTriSurface

Description
    This simple test shows that triSurface loads obj file as triangle faces

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "cellShape.H"
#include "triSurfaceMesh.H"
#include "fvCFD.H"

int main(int argc, char *argv[])
{
    
    triSurface surf("cube.obj");
    Info << "Surface No. of faces: " << surf.faces().size() << endl;
    const faceList& surfaceFaces(surf.faces());
    forAll(surfaceFaces, i)
    {
        Info << "Face " << i << ": " << surfaceFaces[i] << endl;
    }

    OFstream outFile("write_cube.obj");
    forAll(surf.points(), i)
    {
        outFile << "v " 
              << surf.points()[i].x() << " "
              << surf.points()[i].y() << " "
              << surf.points()[i].z() << endl;
    }
    forAll(surf.faces(), i)
    {
        face facei = surf.faces()[i];
        outFile << "f ";
        forAll(facei, plabel)
        {
          outFile << facei[plabel]+1 << " ";
        }
        outFile << endl;
    }
              


    // Create a simple blockMesh without boundaries
    return 0;
}
