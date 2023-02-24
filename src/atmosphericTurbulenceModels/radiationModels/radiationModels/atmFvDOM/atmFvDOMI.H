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

inline const Foam::radiationModels::atmRadiativeIntensityRay&
Foam::radiationModels::atmFvDOM::IRay(const label rayI) const
{
    return  IRay_[rayI];
}


inline const Foam::volScalarField&
Foam::radiationModels::atmFvDOM::IRayLambda
(
    const label rayI,
    const label lambdaI
) const
{
    return IRay_[rayI].ILambda(lambdaI);
}


inline Foam::label Foam::radiationModels::atmFvDOM::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::radiationModels::atmFvDOM::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::radiationModels::atmFvDOM::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::radiationModels::atmFvDOM::nLambda() const
{
    return nLambda_;
}


inline const Foam::volScalarField& Foam::radiationModels::atmFvDOM::a() const
{
    return a_;
}


inline const Foam::volScalarField& Foam::radiationModels::atmFvDOM::aLambda
(
    const label lambdaI
) const
{
    return aLambda_[lambdaI];
}


inline const Foam::volScalarField& Foam::radiationModels::atmFvDOM::G() const
{
    return G_;
}


inline const Foam::volScalarField& Foam::radiationModels::atmFvDOM::qr() const
{
    return qr_;
}

inline const Foam::volScalarField& Foam::radiationModels::atmFvDOM::qin() const
{
    return qin_;
}


inline const Foam::volScalarField& Foam::radiationModels::atmFvDOM::qem() const
{
    return qem_;
}


inline const Foam::radiationModels::blackBodyEmission&
Foam::radiationModels::atmFvDOM::blackBody() const
{
    return blackBody_;
}


inline Foam::scalar Foam::radiationModels::atmFvDOM::omegaMax() const
{
    return omegaMax_;
}


// ************************************************************************* //
