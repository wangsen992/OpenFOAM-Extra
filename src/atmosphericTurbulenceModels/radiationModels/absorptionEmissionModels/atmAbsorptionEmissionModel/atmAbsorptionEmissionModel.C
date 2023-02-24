/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "atmAbsorptionEmissionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace atmRadiationModels
    {
        defineTypeNameAndDebug(atmAbsorptionEmissionModel, 0);
        defineRunTimeSelectionTable(atmAbsorptionEmissionModel, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmRadiationModels::atmAbsorptionEmissionModel::atmAbsorptionEmissionModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::atmRadiationModels::atmAbsorptionEmissionModel::~atmAbsorptionEmissionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::atmRadiationModels::atmAbsorptionEmissionModel::a
(
  const label bandI,
  const vector dir
) const
{
    return aDisp(bandI, dir) + aCont(bandI, dir);
}


Foam::tmp<Foam::volScalarField>
Foam::atmRadiationModels::atmAbsorptionEmissionModel::aCont
(
  const label bandI,
  const vector dir
) const
{
    return volScalarField::New
    (
        "aCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::atmRadiationModels::atmAbsorptionEmissionModel::aDisp
(
  const label bandI,
  const vector dir
) const
{
    return volScalarField::New
    (
        "aDisp",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::atmRadiationModels::atmAbsorptionEmissionModel::e
(
  const label bandI,
  const vector dir
) const
{
    return eDisp(bandI,dir) + eCont(bandI,dir);
}


Foam::tmp<Foam::volScalarField>
Foam::atmRadiationModels::atmAbsorptionEmissionModel::eCont
(
  const label bandI,
  const vector dir
) const
{
    return volScalarField::New
    (
        "eCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::atmRadiationModels::atmAbsorptionEmissionModel::eDisp
(
  const label bandI,
  const vector dir
) const
{
    return volScalarField::New
    (
        "eDisp",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::atmRadiationModels::atmAbsorptionEmissionModel::E
(
  const label bandI,
  const vector dir
) const
{
    return EDisp(bandI,dir) + ECont(bandI,dir);
}


Foam::tmp<Foam::volScalarField>
Foam::atmRadiationModels::atmAbsorptionEmissionModel::ECont
(
  const label bandI,
  const vector dir
) const
{
    return volScalarField::New
    (
        "ECont",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::atmRadiationModels::atmAbsorptionEmissionModel::EDisp
(
  const label bandI,
  const vector dir
) const
{
    return volScalarField::New
    (
        "EDisp",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    );
}


Foam::label Foam::atmRadiationModels::atmAbsorptionEmissionModel::nBands() const
{
    return pTraits<label>::one;
}


const Foam::Vector2D<Foam::scalar>&
Foam::atmRadiationModels::atmAbsorptionEmissionModel::bands(const label n) const
{
    return Vector2D<scalar>::one;
}


bool Foam::atmRadiationModels::atmAbsorptionEmissionModel::isGrey() const
{
    return false;
}


void Foam::atmRadiationModels::atmAbsorptionEmissionModel::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aj
) const
{
    a = this->a();
    aj[0] =  a;
}


// ************************************************************************* //
