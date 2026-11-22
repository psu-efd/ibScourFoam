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

#include "cutGeometry.H"

// * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

using namespace Foam;

//- Fraction t along a -> b of the zero crossing (pa > 0 >= pb)
inline scalar crossing(const scalar pa, const scalar pb)
{
    return pa/(pa - pb);
}

//- Signed volume of a tetrahedron [m3]
inline scalar tetVolume
(
    const point& a,
    const point& b,
    const point& c,
    const point& d
)
{
    return ((b - a) & ((c - a) ^ (d - a)))/6.0;
}

//- Add the volume and first moment of a tetrahedron to V and M
inline void addTetMoment
(
    const point& a,
    const point& b,
    const point& c,
    const point& d,
    scalar& V,
    vector& M
)
{
    const scalar v = mag(tetVolume(a, b, c, d));
    V += v;
    M += 0.25*v*(a + b + c + d);
}


//- Average of the point values
inline scalar pointAverage
(
    const labelUList& labels,
    const scalarField& values
)
{
    scalar sum = 0;
    for (const label i : labels)
    {
        sum += values[i];
    }
    return sum/labels.size();
}

//- Add the facet of the zero level set inside a tetrahedron (phi linear)
//- to the area-vector sum S (pointing from positive to negative side), the
//- area-weighted centre sum M and the area sum A
inline void addTetCutFacet
(
    const point& a,
    const point& b,
    const point& c,
    const point& d,
    const scalar pa,
    const scalar pb,
    const scalar pc,
    const scalar pd,
    vector& S,
    vector& M,
    scalar& A
)
{
    const point x[4] = {a, b, c, d};
    const scalar p[4] = {pa, pb, pc, pd};

    label pos[4];
    label neg[4];
    label nPos = 0;
    label nNeg = 0;

    for (label i = 0; i < 4; ++i)
    {
        if (p[i] > 0)
        {
            pos[nPos++] = i;
        }
        else
        {
            neg[nNeg++] = i;
        }
    }

    if (nPos == 0 || nPos == 4)
    {
        return;
    }

    // Polygon of the facet (triangle, or the quadrilateral of the prism
    // case) and a vector from the positive to the negative side
    point y[4];
    label nY = 0;
    vector toNeg(Zero);

    if (nPos == 1 || nPos == 3)
    {
        // Crossings on the three edges from the odd vertex
        const label i = (nPos == 1 ? pos[0] : neg[0]);
        const label* others = (nPos == 1 ? neg : pos);

        if (nPos == 3 && p[i] >= 0)
        {
            return;
        }

        for (label k = 0; k < 3; ++k)
        {
            const scalar f =
            (
                nPos == 1
              ? crossing(p[i], p[others[k]])
              : crossing(-p[i], -p[others[k]])
            );
            y[nY++] = x[i] + f*(x[others[k]] - x[i]);
        }
        toNeg = (nPos == 1 ? x[neg[0]] - x[i] : x[i] - x[pos[0]]);
    }
    else
    {
        const label i = pos[0];
        const label j = pos[1];
        const label k = neg[0];
        const label l = neg[1];

        auto cross = [&](const label s, const label t) -> point
        {
            return x[s] + crossing(p[s], p[t])*(x[t] - x[s]);
        };

        // Cycle i-k, i-l, j-l, j-k: consecutive crossings share a vertex
        y[nY++] = cross(i, k);
        y[nY++] = cross(i, l);
        y[nY++] = cross(j, l);
        y[nY++] = cross(j, k);
        toNeg = x[k] - x[i];
    }

    // Area vector and area-weighted centre of the polygon (fan from y[0])
    vector Sp(Zero);
    vector Mp(Zero);
    scalar Ap = 0;

    for (label t = 1; t + 1 < nY; ++t)
    {
        const vector St(0.5*((y[t] - y[0]) ^ (y[t + 1] - y[0])));
        const scalar At = mag(St);
        Sp += St;
        Mp += At*(y[0] + y[t] + y[t + 1])/3.0;
        Ap += At;
    }

    if ((Sp & toNeg) < 0)
    {
        Sp = -Sp;
    }

    S += Sp;
    M += Mp;
    A += Ap;
}

//- Add the area and the area-weighted centroid of the positive part of a
//- triangle (phi linear) to A and M
inline void addTrianglePositiveMoment
(
    const point& a,
    const point& b,
    const point& c,
    const scalar pa,
    const scalar pb,
    const scalar pc,
    scalar& A,
    vector& M
)
{
    const point x[3] = {a, b, c};
    const scalar p[3] = {pa, pb, pc};

    // Clipped polygon: positive vertices and the edge crossings
    point y[5];
    label nY = 0;

    for (label i = 0; i < 3; ++i)
    {
        const label j = (i + 1) % 3;

        if (p[i] > 0)
        {
            y[nY++] = x[i];
        }
        if ((p[i] > 0) != (p[j] > 0))
        {
            const scalar f = p[i]/(p[i] - p[j]);
            y[nY++] = x[i] + f*(x[j] - x[i]);
        }
    }

    for (label t = 1; t + 1 < nY; ++t)
    {
        const scalar At = 0.5*mag((y[t] - y[0]) ^ (y[t + 1] - y[0]));
        A += At;
        M += At*(y[0] + y[t] + y[t + 1])/3.0;
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::cutGeometry::triangleFraction
(
    const scalar pa,
    const scalar pb,
    const scalar pc
)
{
    const scalar p[3] = {pa, pb, pc};
    label nPos = 0;
    label iPos = -1;
    label iNeg = -1;

    for (label i = 0; i < 3; ++i)
    {
        if (p[i] > 0)
        {
            ++nPos;
            iPos = i;
        }
        else
        {
            iNeg = i;
        }
    }

    if (nPos == 0)
    {
        return 0;
    }
    if (nPos == 3)
    {
        return 1;
    }

    if (nPos == 1)
    {
        // Corner triangle at the positive vertex
        const scalar q = p[iPos];
        const scalar r1 = p[(iPos + 1) % 3];
        const scalar r2 = p[(iPos + 2) % 3];
        return crossing(q, r1)*crossing(q, r2);
    }

    // Two positive vertices: complement of the corner at the negative one
    const scalar q = -p[iNeg];
    const scalar r1 = -p[(iNeg + 1) % 3];
    const scalar r2 = -p[(iNeg + 2) % 3];

    // q >= 0 > r1, r2 (the negative vertex may be exactly 0)
    if (q <= 0)
    {
        return 1;
    }
    return 1 - crossing(q, r1)*crossing(q, r2);
}


Foam::scalar Foam::cutGeometry::tetPositiveVolume
(
    const point& a,
    const point& b,
    const point& c,
    const point& d,
    const scalar pa,
    const scalar pb,
    const scalar pc,
    const scalar pd
)
{
    const point x[4] = {a, b, c, d};
    const scalar p[4] = {pa, pb, pc, pd};
    const scalar V = mag(tetVolume(a, b, c, d));

    label pos[4];
    label neg[4];
    label nPos = 0;
    label nNeg = 0;

    for (label i = 0; i < 4; ++i)
    {
        if (p[i] > 0)
        {
            pos[nPos++] = i;
        }
        else
        {
            neg[nNeg++] = i;
        }
    }

    if (nPos == 0)
    {
        return 0;
    }
    if (nPos == 4)
    {
        return V;
    }

    if (nPos == 1)
    {
        const label i = pos[0];
        scalar f = 1;
        for (label k = 0; k < 3; ++k)
        {
            f *= crossing(p[i], p[neg[k]]);
        }
        return f*V;
    }

    if (nPos == 3)
    {
        const label i = neg[0];
        if (p[i] >= 0)
        {
            return V;
        }
        scalar f = 1;
        for (label k = 0; k < 3; ++k)
        {
            f *= crossing(-p[i], -p[pos[k]]);
        }
        return (1 - f)*V;
    }

    // Two positive vertices (i, j), two negative (k, l): the positive part
    // is the prism with the triangles (i, x_ik, x_il) and (j, x_jk, x_jl)
    const label i = pos[0];
    const label j = pos[1];
    const label k = neg[0];
    const label l = neg[1];

    const point xik = x[i] + crossing(p[i], p[k])*(x[k] - x[i]);
    const point xil = x[i] + crossing(p[i], p[l])*(x[l] - x[i]);
    const point xjk = x[j] + crossing(p[j], p[k])*(x[k] - x[j]);
    const point xjl = x[j] + crossing(p[j], p[l])*(x[l] - x[j]);

    // Prism (0, 1, 2 | 3, 4, 5) = (i, xik, xil | j, xjk, xjl) as three
    // tetrahedra (0,1,2,3), (1,2,3,4), (2,3,4,5)
    return
        mag(tetVolume(x[i], xik, xil, x[j]))
      + mag(tetVolume(xik, xil, x[j], xjk))
      + mag(tetVolume(xil, x[j], xjk, xjl));
}


void Foam::cutGeometry::tetPositiveMoment
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
)
{
    const point x[4] = {a, b, c, d};
    const scalar p[4] = {pa, pb, pc, pd};

    label pos[4];
    label neg[4];
    label nPos = 0;
    label nNeg = 0;

    for (label i = 0; i < 4; ++i)
    {
        if (p[i] > 0)
        {
            pos[nPos++] = i;
        }
        else
        {
            neg[nNeg++] = i;
        }
    }

    if (nPos == 0)
    {
        return;
    }
    if (nPos == 4)
    {
        addTetMoment(a, b, c, d, V, M);
        return;
    }

    if (nPos == 1)
    {
        // Corner tetrahedron at the positive vertex
        const label i = pos[0];
        point y[3];
        for (label k = 0; k < 3; ++k)
        {
            y[k] = x[i] + crossing(p[i], p[neg[k]])*(x[neg[k]] - x[i]);
        }
        addTetMoment(x[i], y[0], y[1], y[2], V, M);
        return;
    }

    if (nPos == 3)
    {
        // Whole tetrahedron minus the corner at the negative vertex
        const label i = neg[0];
        addTetMoment(a, b, c, d, V, M);
        if (p[i] >= 0)
        {
            return;
        }
        point y[3];
        for (label k = 0; k < 3; ++k)
        {
            y[k] = x[i] + crossing(-p[i], -p[pos[k]])*(x[pos[k]] - x[i]);
        }
        scalar Vc = 0;
        vector Mc(Zero);
        addTetMoment(x[i], y[0], y[1], y[2], Vc, Mc);
        V -= Vc;
        M -= Mc;
        return;
    }

    // Two positive vertices: the prism of tetPositiveVolume
    const label i = pos[0];
    const label j = pos[1];
    const label k = neg[0];
    const label l = neg[1];

    const point xik = x[i] + crossing(p[i], p[k])*(x[k] - x[i]);
    const point xil = x[i] + crossing(p[i], p[l])*(x[l] - x[i]);
    const point xjk = x[j] + crossing(p[j], p[k])*(x[k] - x[j]);
    const point xjl = x[j] + crossing(p[j], p[l])*(x[l] - x[j]);

    addTetMoment(x[i], xik, xil, x[j], V, M);
    addTetMoment(xik, xil, x[j], xjk, V, M);
    addTetMoment(xil, x[j], xjk, xjl, V, M);
}


Foam::scalar Foam::cutGeometry::faceFraction
(
    const primitiveMesh& mesh,
    const label facei,
    const scalarField& pointPhi
)
{
    const face& f = mesh.faces()[facei];
    const pointField& pts = mesh.points();

    bool anyPos = false;
    bool anyNeg = false;
    for (const label pointi : f)
    {
        (pointPhi[pointi] > 0 ? anyPos : anyNeg) = true;
    }

    if (!anyNeg)
    {
        return 1;
    }
    if (!anyPos)
    {
        return 0;
    }

    const point fc(f.average(pts));
    const scalar pc = pointAverage(f, pointPhi);

    scalar sumA = 0;
    scalar sumPosA = 0;

    forAll(f, i)
    {
        const label a = f[i];
        const label b = f.nextLabel(i);
        const scalar A = 0.5*mag((pts[a] - fc) ^ (pts[b] - fc));

        sumA += A;
        sumPosA += A*triangleFraction(pc, pointPhi[a], pointPhi[b]);
    }

    return (sumA > VSMALL ? sumPosA/sumA : 0);
}


Foam::scalar Foam::cutGeometry::cellFraction
(
    const primitiveMesh& mesh,
    const label celli,
    const scalarField& pointPhi
)
{
    const labelList& cPoints = mesh.cellPoints()[celli];

    bool anyPos = false;
    bool anyNeg = false;
    for (const label pointi : cPoints)
    {
        (pointPhi[pointi] > 0 ? anyPos : anyNeg) = true;
    }

    if (!anyNeg)
    {
        return 1;
    }
    if (!anyPos)
    {
        return 0;
    }

    const pointField& pts = mesh.points();
    const faceList& faces = mesh.faces();
    const point cc(mesh.cellCentres()[celli]);
    const scalar pcc = pointAverage(cPoints, pointPhi);

    scalar sumV = 0;
    scalar sumPosV = 0;

    for (const label facei : mesh.cells()[celli])
    {
        const face& f = faces[facei];
        const point fc(f.average(pts));
        const scalar pfc = pointAverage(f, pointPhi);

        forAll(f, i)
        {
            const label a = f[i];
            const label b = f.nextLabel(i);

            sumV += mag(tetVolume(cc, fc, pts[a], pts[b]));
            sumPosV +=
                tetPositiveVolume
                (
                    cc, fc, pts[a], pts[b],
                    pcc, pfc, pointPhi[a], pointPhi[b]
                );
        }
    }

    return (sumV > VSMALL ? min(sumPosV/sumV, scalar(1)) : 0);
}


Foam::point Foam::cutGeometry::cellFluidCentre
(
    const primitiveMesh& mesh,
    const label celli,
    const scalarField& pointPhi
)
{
    const labelList& cPoints = mesh.cellPoints()[celli];
    const point cc(mesh.cellCentres()[celli]);

    bool anyPos = false;
    bool anyNeg = false;
    for (const label pointi : cPoints)
    {
        (pointPhi[pointi] > 0 ? anyPos : anyNeg) = true;
    }

    if (!anyNeg || !anyPos)
    {
        return cc;
    }

    const pointField& pts = mesh.points();
    const faceList& faces = mesh.faces();
    const scalar pcc = pointAverage(cPoints, pointPhi);

    scalar V = 0;
    vector M(Zero);

    for (const label facei : mesh.cells()[celli])
    {
        const face& f = faces[facei];
        const point fc(f.average(pts));
        const scalar pfc = pointAverage(f, pointPhi);

        forAll(f, i)
        {
            const label a = f[i];
            const label b = f.nextLabel(i);

            tetPositiveMoment
            (
                cc, fc, pts[a], pts[b],
                pcc, pfc, pointPhi[a], pointPhi[b],
                V, M
            );
        }
    }

    return (V > VSMALL ? M/V : cc);
}


bool Foam::cutGeometry::cellWallElement
(
    const primitiveMesh& mesh,
    const label celli,
    const scalarField& pointPhi,
    vector& Sw,
    point& xw
)
{
    const labelList& cPoints = mesh.cellPoints()[celli];
    const point cc(mesh.cellCentres()[celli]);

    Sw = Zero;
    xw = cc;

    bool anyPos = false;
    bool anyNeg = false;
    for (const label pointi : cPoints)
    {
        (pointPhi[pointi] > 0 ? anyPos : anyNeg) = true;
    }

    if (!anyNeg || !anyPos)
    {
        return false;
    }

    const pointField& pts = mesh.points();
    const faceList& faces = mesh.faces();
    const scalar pcc = pointAverage(cPoints, pointPhi);

    vector M(Zero);
    scalar A = 0;

    for (const label facei : mesh.cells()[celli])
    {
        const face& f = faces[facei];
        const point fc(f.average(pts));
        const scalar pfc = pointAverage(f, pointPhi);

        forAll(f, i)
        {
            const label a = f[i];
            const label b = f.nextLabel(i);

            addTetCutFacet
            (
                cc, fc, pts[a], pts[b],
                pcc, pfc, pointPhi[a], pointPhi[b],
                Sw, M, A
            );
        }
    }

    if (A > VSMALL)
    {
        xw = M/A;
        return true;
    }

    Sw = Zero;
    return false;
}


Foam::point Foam::cutGeometry::faceOpenCentre
(
    const primitiveMesh& mesh,
    const label facei,
    const scalarField& pointPhi
)
{
    const face& f = mesh.faces()[facei];
    const pointField& pts = mesh.points();
    const point fc(f.average(pts));

    bool anyPos = false;
    bool anyNeg = false;
    for (const label pointi : f)
    {
        (pointPhi[pointi] > 0 ? anyPos : anyNeg) = true;
    }

    if (!anyNeg || !anyPos)
    {
        return fc;
    }

    const scalar pc = pointAverage(f, pointPhi);

    scalar A = 0;
    vector M(Zero);

    forAll(f, i)
    {
        const label a = f[i];
        const label b = f.nextLabel(i);

        addTrianglePositiveMoment
        (
            fc, pts[a], pts[b], pc, pointPhi[a], pointPhi[b], A, M
        );
    }

    return (A > VSMALL ? M/A : fc);
}


// ************************************************************************* //
