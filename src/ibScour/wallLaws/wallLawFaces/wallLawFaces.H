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

Namespace
    Foam::wallLawFaces

Description
    Wall law (ibWallLaw) evaluated at the faces of a body-fitted wall patch
    from the velocity of the first cell, for the body-fitted wall
    conditions (nutWallLaw, wallLawTurbulence):

        - y: normal distance of the cell centre from the face, at least
          0.1 of the centre-to-face distance (distorted cells);
        - U_t: tangential velocity of the cell relative to the wall;
        - u_tau from |U_t| at y (ibWallLaw::uTau).

SourceFiles
    wallLawFaces.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_wallLawFaces_H
#define Foam_wallLawFaces_H

#include "ibWallLaw.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLawFaces
{

//- Friction velocity uTau [m/s], wall distance y [m] and tangential
//- relative velocity magnitude magUt [m/s] of the faces of the patch of Up,
//- with kinematic viscosity nuw [m2/s] per face. Optionally counts the
//- faces whose distance was bounded (n.d < minRatio |d|, distorted cells)
//- and the smallest ratio n.d/|d| (local, not reduced).
void evaluate
(
    const ibWallLaw& law,
    const fvPatchVectorField& Up,
    const scalarField& nuw,
    scalarField& uTau,
    scalarField& y,
    scalarField& magUt,
    label* nClamped = nullptr,
    scalar* minRatio = nullptr
);

//- Lower bound of the wall distance as a fraction of the centre-to-face
//- distance |d|
constexpr scalar yMinFraction = 0.1;

} // End namespace wallLawFaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
