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

#include "sphereSolid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sphereSolid, 0);
    addToRunTimeSelectionTable(ibSolid, sphereSolid, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sphereSolid::sphereSolid
(
    const word& name,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    ibSolid(name, dict),
    centre_(dict.get<point>("centre")),
    radius_(dict.get<scalar>("radius"))
{
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::sphereSolid::bodyDistance
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
        const vector r(bodyPoints[i] - centre_);
        const scalar magR = mag(r);

        distance[i] = magR - radius_;
        normal[i] = (magR > VSMALL ? r/magR : vector(0, 0, 1));
    }
}


// ************************************************************************* //
