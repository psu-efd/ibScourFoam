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

#include "cylinderSolid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cylinderSolid, 0);
    addToRunTimeSelectionTable(ibSolid, cylinderSolid, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylinderSolid::cylinderSolid
(
    const word& name,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    ibSolid(name, dict),
    point1_(dict.get<point>("point1")),
    axis_(dict.get<point>("point2") - point1_),
    length_(mag(axis_)),
    radius_(dict.get<scalar>("radius")),
    finite_(dict.getOrDefault<bool>("finite", false))
{
    if (length_ < VSMALL || radius_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Cylinder " << name << " needs distinct axis points and a"
            << " positive radius" << exit(FatalIOError);
    }

    axis_ /= length_;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::cylinderSolid::bodyDistance
(
    const pointField& bodyPoints,
    scalarField& distance,
    vectorField& normal
) const
{
    distance.resize(bodyPoints.size());
    normal.resize(bodyPoints.size());

    // Any unit vector normal to the axis, for points on the axis itself
    vector perp(axis_ ^ vector(1, 0, 0));
    if (mag(perp) < 0.1)
    {
        perp = axis_ ^ vector(0, 1, 0);
    }
    perp.normalise();

    forAll(bodyPoints, i)
    {
        const vector v(bodyPoints[i] - point1_);
        const scalar t = (v & axis_);
        const vector rVec(v - t*axis_);
        const scalar r = mag(rVec);
        const vector rHat(r > VSMALL ? rVec/r : perp);

        // Radial signed distance
        const scalar dr = r - radius_;

        if (!finite_)
        {
            distance[i] = dr;
            normal[i] = rHat;
            continue;
        }

        // Axial signed distance to the caps and outward axial direction
        const scalar da = max(-t, t - length_);
        const vector aHat(t < 0.5*length_ ? -axis_ : axis_);

        if (dr <= 0 && da <= 0)
        {
            // Inside: nearest of the mantle and the caps
            if (dr > da)
            {
                distance[i] = dr;
                normal[i] = rHat;
            }
            else
            {
                distance[i] = da;
                normal[i] = aHat;
            }
        }
        else
        {
            // Outside: distance to the mantle, cap or rim
            const scalar pr = max(dr, scalar(0));
            const scalar pa = max(da, scalar(0));
            distance[i] = Foam::sqrt(sqr(pr) + sqr(pa));
            normal[i] = normalised(pr*rHat + pa*aHat);
        }
    }
}


// ************************************************************************* //
