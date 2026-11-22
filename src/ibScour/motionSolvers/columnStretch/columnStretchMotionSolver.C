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

#include "columnStretchMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "valuePointPatchFields.H"
#include "SortList.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(columnStretchMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        columnStretchMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        columnStretchMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::columnStretchMotionSolver::bedTable
(
    scalarField& x,
    scalarField& zb0,
    scalarField& dz
) const
{
    const polyPatch& pp = mesh().boundaryMesh()[bedPatchID_];
    const labelList& mp = pp.meshPoints();

    const pointPatchVectorField& bcBase =
        pointDisplacement_.boundaryField()[bedPatchID_];

    if (!isA<valuePointPatchVectorField>(bcBase))
    {
        FatalErrorInFunction
            << "pointDisplacement of patch " << bedPatchName_
            << " must hold values (fixedValue), found type "
            << bcBase.type() << exit(FatalError);
    }

    const vectorField& disp = refCast<const valuePointPatchVectorField>(bcBase);

    scalarField xl(mp.size());
    scalarField zl(mp.size());
    scalarField dl(mp.size());
    forAll(mp, i)
    {
        xl[i] = points0()[mp[i]].x();
        zl[i] = points0()[mp[i]].z();
        dl[i] = disp[i].z();
    }

    // All processors need the whole bed
    List<scalarField> xAll(Pstream::nProcs());
    List<scalarField> zAll(Pstream::nProcs());
    List<scalarField> dAll(Pstream::nProcs());
    xAll[Pstream::myProcNo()] = xl;
    zAll[Pstream::myProcNo()] = zl;
    dAll[Pstream::myProcNo()] = dl;
    Pstream::allGatherList(xAll);
    Pstream::allGatherList(zAll);
    Pstream::allGatherList(dAll);

    label n = 0;
    forAll(xAll, proci)
    {
        n += xAll[proci].size();
    }

    scalarField xg(n);
    scalarField zg(n);
    scalarField dg(n);
    n = 0;
    forAll(xAll, proci)
    {
        forAll(xAll[proci], i)
        {
            xg[n] = xAll[proci][i];
            zg[n] = zAll[proci][i];
            dg[n] = dAll[proci][i];
            ++n;
        }
    }

    if (n == 0)
    {
        FatalErrorInFunction
            << "Bed patch " << bedPatchName_ << " has no points"
            << exit(FatalError);
    }

    const SortList<scalar> order(xg);
    const labelList& idx = order.indices();

    // Merge points of the same column (front and back of the 2-D mesh)
    DynamicList<scalar> xs(n);
    DynamicList<scalar> zs(n);
    DynamicList<scalar> ds(n);
    forAll(idx, k)
    {
        const label i = idx[k];
        if (xs.size() && mag(xg[i] - xs.last()) < xTol_)
        {
            continue;
        }
        xs.append(xg[i]);
        zs.append(zg[i]);
        ds.append(dg[i]);
    }

    x.transfer(xs);
    zb0.transfer(zs);
    dz.transfer(ds);
}


Foam::scalar Foam::columnStretchMotionSolver::interpolate
(
    const scalar x,
    const scalarField& xs,
    const scalarField& ys
)
{
    if (x <= xs.first())
    {
        return ys.first();
    }
    if (x >= xs.last())
    {
        return ys.last();
    }

    // Binary search for the interval
    label lo = 0;
    label hi = xs.size() - 1;
    while (hi - lo > 1)
    {
        const label mid = (lo + hi)/2;
        if (xs[mid] <= x)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }
    const scalar w = (x - xs[lo])/max(xs[hi] - xs[lo], VSMALL);
    return (1 - w)*ys[lo] + w*ys[hi];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::columnStretchMotionSolver::columnStretchMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    bedPatchName_(coeffDict().get<word>("bedPatch")),
    bedPatchID_(mesh.boundaryMesh().findPatchID(bedPatchName_)),
    xFace_(coeffDict().get<scalar>("xFace")),
    apronLevel_(coeffDict().get<scalar>("apronLevel")),
    blendLength_(coeffDict().get<scalar>("blendLength")),
    minGapFraction_(coeffDict().getOrDefault<scalar>("minGapFraction", 0.25)),
    zTop_(gMax(points0().component(vector::Z))),
    xTol_(1e-6*mesh.bounds().mag())
{
    if (bedPatchID_ < 0)
    {
        FatalIOErrorInFunction(coeffDict())
            << "Bed patch " << bedPatchName_ << " not found"
            << exit(FatalIOError);
    }
}


Foam::columnStretchMotionSolver::columnStretchMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    bedPatchName_(coeffDict().get<word>("bedPatch")),
    bedPatchID_(mesh.boundaryMesh().findPatchID(bedPatchName_)),
    xFace_(coeffDict().get<scalar>("xFace")),
    apronLevel_(coeffDict().get<scalar>("apronLevel")),
    blendLength_(coeffDict().get<scalar>("blendLength")),
    minGapFraction_(coeffDict().getOrDefault<scalar>("minGapFraction", 0.25)),
    zTop_(gMax(this->points0().component(vector::Z))),
    xTol_(1e-6*mesh.bounds().mag())
{
    if (bedPatchID_ < 0)
    {
        FatalIOErrorInFunction(coeffDict())
            << "Bed patch " << bedPatchName_ << " not found"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::columnStretchMotionSolver::curPoints() const
{
    auto tcurPoints = tmp<pointField>::New
    (
        points0() + pointDisplacement_.primitiveField()
    );
    twoDCorrectPoints(tcurPoints.ref());
    return tcurPoints;
}


void Foam::columnStretchMotionSolver::solve()
{
    scalarField xb, zb0, dzb;
    bedTable(xb, zb0, dzb);

    vectorField& disp = pointDisplacement_.primitiveFieldRef();
    const pointField& p0 = points0();

    forAll(disp, pointi)
    {
        const scalar x = p0[pointi].x();
        const scalar z0 = p0[pointi].z();

        disp[pointi] = Zero;

        if (x < xFace_ - xTol_)
        {
            continue;
        }

        const scalar db = interpolate(x, xb, dzb);
        const scalar zbed0 = interpolate(x, xb, zb0);
        const scalar gap = apronLevel_ - zbed0;

        scalar w = 1;
        if (blendLength_ > 0)
        {
            w = min(max((x - xFace_)/blendLength_, scalar(0)), scalar(1));
        }
        else if (x < xFace_ + xTol_)
        {
            w = 0;
        }

        auto full = [&](const scalar z) -> scalar
        {
            return db*(zTop_ - z)/max(zTop_ - zbed0, VSMALL);
        };

        if (z0 >= apronLevel_ - xTol_)
        {
            disp[pointi].z() = w*full(z0);
        }
        else
        {
            // Gap block: between the moved bed and the moved apron level
            // of the column, never thinner than minGapFraction of the gap
            scalar dTop = w*full(apronLevel_);
            if (gap > 0)
            {
                dTop = max(dTop, db - (1 - minGapFraction_)*gap);
                const scalar s = (z0 - zbed0)/gap;
                disp[pointi].z() = db + s*(dTop - db);
            }
            else
            {
                disp[pointi].z() = db;
            }
        }
    }

    pointDisplacement_.correctBoundaryConditions();
}


// ************************************************************************* //
