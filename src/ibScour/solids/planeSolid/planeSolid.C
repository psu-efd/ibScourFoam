/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 ibScourFoam developers (Penn State)
-------------------------------------------------------------------------------
License
    This file is part of ibScourFoam, built on OpenFOAM.

    ibScourFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ibScourFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ibScourFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "planeSolid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(planeSolid, 0);
    addToRunTimeSelectionTable(ibSolid, planeSolid, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::planeSolid::planeSolid
(
    const word& name,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    ibSolid(name, dict),
    point_(dict.get<point>("point")),
    normal_(normalised(dict.get<vector>("normal")))
{
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::planeSolid::bodyDistance
(
    const pointField& bodyPoints,
    scalarField& distance,
    vectorField& normal
) const
{
    distance.resize(bodyPoints.size());
    normal.resize(bodyPoints.size());

    forAll(bodyPoints, i)
    {
        distance[i] = ((bodyPoints[i] - point_) & normal_);
        normal[i] = normal_;
    }
}


// ************************************************************************* //
