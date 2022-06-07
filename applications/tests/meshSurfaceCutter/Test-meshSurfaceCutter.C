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
#include "polyMesh.H"
#include "meshTools.H"
#include "argList.H"
#include "DynamicList.H"
#include "meshSearch.H"

#include "indexedOctree.H"
#include "treeDataFace.H"
#include "treeBoundBox.H"

#include "surfaceMeshTools.H"
#include "meshSurfaceCutter.H"


int main(int argc, char *argv[])
{

    argList args(argc, argv);
    Time runTime(Time::controlDictName, args);
    
    polyMesh mesh
    (
      IOobject
      (
        polyMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
      )
    );
    triSurface surf(runTime.path()/"tree.obj");

    // Visualize the obj and mesh as seen in OpenFOAM
    
    // surfaceMeshTools::writeObjFile
    // (
    //   OFstream(runTime.path()/"triSphere.obj"), 
    //   surf
    // );
    // surfaceMeshTools::writeObjFile
    // (
    //     OFstream(runTime.path()/"mesh.obj"),
    //     mesh
    // );

    // Slicing up the surface faces
    // Step 1: Prepare the variables
    const faceList& surfFaces(surf.faces());

    // Step 2: Checking each face of surfFaces
    // For each surfFace, locate all the intersection points with the mesh face
    // Test: Find all intersection faces and write to obj for viewing

    DynamicList<face> cutFaces(0);
    boolList cutFaceMask = surfaceMeshTools::markFaces(mesh, surf);
    forAll(cutFaceMask, i)
    {
        if(cutFaceMask[i])
        {
            cutFaces.append(mesh.faces()[i]);
        }
    }
    OFstream cutFaceObj(runTime.path()/"cutFaces.obj");
    surfaceMeshTools::writeObjFile(cutFaceObj, cutFaces, mesh.points());

    // Test: Find all intersection points and write to obj for viewing
    // Approach: using indexedOctree to find intersections
    // labelList meshFacesList(mesh.nInternalFaces());
    // DynamicList<point> intersectionPoints(0);
    // forAll(meshFacesList, i)
    // {
    //     meshFacesList[i] = i;
    // }
    // indexedOctree<treeDataFace> faceTree
    // (
    //   treeDataFace
    //   (false, mesh, meshFacesList),
    //   treeBoundBox(mesh.points()),
    //   8, 10, 3.00
    // );
    // forAll(surfFaces, i)
    // {
    //     face facei = surfFaces[i];
    //     const pointField& facePoints = facei.points(surf.points());
    //     forAll(facePoints, j)
    //     {
    //         point p1 = facePoints[j];
    //         point p2 = facePoints[facePoints.fcIndex(j)];
    //         pointIndexHit hit = faceTree.findLine(p1, p2);
    //         if (hit.hit())
    //         {
    //             intersectionPoints.append(hit.hitPoint());
    //         }
    //     }
    // }

    surfaceMeshTools::meshSurfaceCutter cutter(mesh, surf);
    pointField& intersectionPoints(cutter.intersectionPoints());
            
    OFstream intersectObj(runTime.path()/"intersectionPts.obj");
    forAll(intersectionPoints, i)
    {
        point p = intersectionPoints[i];
        intersectObj << "v " 
                  << p.x() << " " 
                  << p.y() << " "
                  << p.z() << endl;
    }

    // Outcome: a HashTable containing a key of mesh cells and value of list of
    // surface faces


    return 0;
}
