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

#include "ibSolid.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ibSolid, 0);
    defineRunTimeSelectionTable(ibSolid, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ibSolid::ibSolid(const word& name, const dictionary& dict)
:
    name_(name),
    origin_(dict.getOrDefault<point>("origin", Zero)),
    R_(I),
    velocity_(dict.getOrDefault<vector>("velocity", Zero)),
    omega_(dict.getOrDefault<vector>("omega", Zero)),
    origin0_(origin_),
    displacement_(nullptr),
    dict_(dict)
{
    if (dict.found("displacement"))
    {
        displacement_ = Function1<vector>::New("displacement", dict);
        origin_ = origin0_ + displacement_->value(0);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ibSolid::applyPrescribedMotion(const scalar t, const scalar dt)
{
    if (!displacement_)
    {
        return false;
    }

    const vector d1(displacement_->value(t));
    const vector d0(displacement_->value(max(t - dt, scalar(0))));
    const point o(origin0_ + d1);
    const bool moved = (magSqr(o - origin_) > VSMALL);

    origin_ = o;
    velocity_ = (dt > VSMALL ? (d1 - d0)/dt : vector(Zero));
    return moved;
}


bool Foam::ibSolid::moving() const
{
    return (magSqr(velocity_) > VSMALL || magSqr(omega_) > VSMALL);
}


void Foam::ibSolid::distance
(
    const pointField& points,
    scalarField& distance,
    vectorField& normal
) const
{
    // Global -> body frame: x_b = R^T (x - origin)
    const pointField bodyPoints(R_.T() & (points - origin_));

    bodyDistance(bodyPoints, distance, normal);

    // Body -> global frame for the normals
    normal = (R_ & normal);
}


void Foam::ibSolid::distance
(
    const pointField& points,
    scalarField& distance,
    vectorField& normal,
    const point& origin,
    const tensor& R
) const
{
    const pointField bodyPoints(R.T() & (points - origin));
    bodyDistance(bodyPoints, distance, normal);
    normal = (R & normal);
}


Foam::tmp<Foam::vectorField>
Foam::ibSolid::wallVelocity(const pointField& points) const
{
    return velocity_ + (omega_ ^ (points - origin_));
}


void Foam::ibSolid::setState
(
    const point& origin,
    const tensor& R,
    const vector& velocity,
    const vector& omega
)
{
    origin_ = origin;
    R_ = R;
    velocity_ = velocity;
    omega_ = omega;
}


// ************************************************************************* //
