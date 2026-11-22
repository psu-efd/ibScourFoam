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
    Foam::cutGeometry

Description
    Fractions of faces and cells on the positive side of a level set phi
    given at the mesh points (cut-face immersed boundary, FAVOR-type
    apertures).

    Faces are split into triangles (face centre = average of the face
    points, phi there = average of the point values) and cells into
    tetrahedra (these triangles and the cell centre, phi there = average of
    the cell points). phi is linear on every triangle and tetrahedron, and
    the positive part is computed exactly on each:

        - triangle: area fraction from the edge crossings;
        - tetrahedron: one positive vertex -> corner tetrahedron, three ->
          its complement, two -> the clipped prism (three tetrahedra).

    A point with phi = 0 counts as negative, so a plane through a layer of
    mesh points gives open fractions 0 or 1 exactly.

    The centroid of the positive part of a cell (cellFluidCentre) follows
    from the first moments of the same positive tetrahedra.

SourceFiles
    cutGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cutGeometry_H
#define Foam_cutGeometry_H

#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cutGeometry
{

//- Area fraction of a triangle on the positive side (phi linear)
scalar triangleFraction(const scalar pa, const scalar pb, const scalar pc);

//- Volume of a tetrahedron on the positive side (phi linear) [m3]
scalar tetPositiveVolume
(
    const point& a,
    const point& b,
    const point& c,
    const point& d,
    const scalar pa,
    const scalar pb,
    const scalar pc,
    const scalar pd
);

//- Volume [m3] and first moment [m4] of the positive part of a tetrahedron
//- (phi linear), added to V and M
void tetPositiveMoment
(
    const point& a,
    const point& b,
    const point& c,
    const point& d,
    const scalar pa,
    const scalar pb,
    const scalar pc,
    const scalar pd,
    scalar& V,
    vector& M
);

//- Open (positive) area fraction of face facei
scalar faceFraction
(
    const primitiveMesh& mesh,
    const label facei,
    const scalarField& pointPhi
);

//- Positive volume fraction of cell celli
scalar cellFraction
(
    const primitiveMesh& mesh,
    const label celli,
    const scalarField& pointPhi
);

//- Centroid of the positive part of cell celli (the cell centre for a
//- cell that is not cut)
point cellFluidCentre
(
    const primitiveMesh& mesh,
    const label celli,
    const scalarField& pointPhi
);

//- Centroid of the positive (open) part of face facei (the face centre
//- for a face that is not cut), from the positive parts of the same
//- triangle fan as faceFraction
point faceOpenCentre
(
    const primitiveMesh& mesh,
    const label facei,
    const scalarField& pointPhi
);

//- Wall element of a cut cell: area vector Sw of the zero level set inside
//- the cell, pointing from the positive (fluid) side into the negative
//- (solid) side, and its area-weighted centre xw. The facets of the
//- positive tetrahedra are summed, so for a planar level set Sw equals the
//- closure of the open faces, -sum_f alpha_f S_f. Returns false (Sw = 0,
//- xw = cell centre) for a cell that is not cut.
bool cellWallElement
(
    const primitiveMesh& mesh,
    const label celli,
    const scalarField& pointPhi,
    vector& Sw,
    point& xw
);

} // End namespace cutGeometry
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
