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

#include "triSurfaceSolid.H"
#include "addToRunTimeSelectionTable.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(triSurfaceSolid, 0);
    addToRunTimeSelectionTable(ibSolid, triSurfaceSolid, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceSolid::triSurfaceSolid
(
    const word& name,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    ibSolid(name, dict),
    surface_
    (
        IOobject
        (
            dict.get<fileName>("file"),
            obr.time().constant(),
            "triSurface",
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    ),
    closed_(dict.getOrDefault<bool>("closed", surface_.hasVolumeType())),
    normalSign_(dict.getOrDefault<bool>("flipNormals", false) ? -1 : 1)
{
    Info<< "        " << surface_.size() << " triangles, "
        << (closed_ ? "closed" : "open") << " surface" << endl;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::triSurfaceSolid::bodyDistance
(
    const pointField& bodyPoints,
    scalarField& distance,
    vectorField& normal
) const
{
    distance.resize(bodyPoints.size());
    normal.resize(bodyPoints.size());

    // Unlimited search radius: every point gets its nearest triangle
    const scalarField nearestDistSqr(bodyPoints.size(), GREAT);

    List<pointIndexHit> hits;
    surface_.findNearest(bodyPoints, nearestDistSqr, hits);

    vectorField faceNormals;
    surface_.getNormal(hits, faceNormals);

    List<volumeType> volType;
    if (closed_)
    {
        surface_.getVolumeType(bodyPoints, volType);
    }

    forAll(bodyPoints, i)
    {
        const vector toPoint(bodyPoints[i] - hits[i].point());
        const scalar d = mag(toPoint);

        // Face normal oriented into the fluid
        vector n(normalSign_*faceNormals[i]);

        scalar sign = 1;
        if (closed_)
        {
            sign = (volType[i] == volumeType::INSIDE ? -1 : 1);

            // Closed surfaces: make n point away from the interior
            if ((toPoint & n)*sign < 0)
            {
                n = -n;
            }
        }
        else
        {
            sign = ((toPoint & n) >= 0 ? 1 : -1);
        }

        distance[i] = sign*d;
        normal[i] = n;
    }
}


// ************************************************************************* //
