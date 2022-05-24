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
    Immersed Boundary Method implementation using fvModel containing the analytical engine for surface analysis. Multiple methods can be included, and should support multiple variable updates. Support for turbulence model should include both k-eps and reynold-stress models

Usage
    Example usage:
    \verbatim
    immersedBoundary
    {
        type        immersedBoundary;

        surface           U;            // Name of the velocity field
        Ug          (10 0 0);     // Value of geostrophic velocity
        f           1e-4;         // Value of Coriolis Parameter
    }
    \endverbatim

SourceFiles
    geostrophicForce.C

\*---------------------------------------------------------------------------*/
