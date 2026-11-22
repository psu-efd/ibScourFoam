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

#include "profiling.H"
#include "immersedBoundary.H"
#include "quaternion.H"
#include "ibWallModel.H"
#include "fvMatrices.H"
#include "surfaceInterpolate.H"
#include "syncTools.H"
#include "wallDist.H"
#include "regionSplit.H"
#include "bitSet.H"
#include "fvcGrad.H"
#include "pimpleControl.H"
#include "fvm.H"
#include "fvc.H"
#include "cutGeometry.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundary, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::immersedBoundary::readSolids()
{
    const dictionary& solidsDict = dict_.subOrEmptyDict("solids");

    solids_.resize(solidsDict.size());

    label i = 0;
    for (const entry& e : solidsDict)
    {
        if (!e.isDict())
        {
            continue;
        }

        solids_.set(i++, ibSolid::New(e.keyword(), e.dict(), mesh()));
    }

    solids_.resize(i);

    if (dict_.found("bed"))
    {
        bed_ = bedModel::New(mesh(), dict_.subDict("bed"));
    }

    if (solids_.empty() && !bed_)
    {
        FatalIOErrorInFunction(dict_)
            << "No solids defined in " << dict_.objectPath()
            << exit(FatalIOError);
    }
}


void Foam::immersedBoundary::classify()
{
    const fvMesh& m = mesh();
    const pointField& C = m.C().primitiveField();

    scalarField& d = distance_.primitiveFieldRef();
    vectorField& n = normal_.primitiveFieldRef();

    d = GREAT;
    n = vector(0, 0, 1);
    nearestSolid_ = -1;

    // Nearest solid surface of every cell centre. A point inside any solid
    // gets the most negative distance (deepest inside).
    scalarField ds;
    vectorField ns;

    forAll(solids_, solidi)
    {
        solids_[solidi].distance(C, ds, ns);

        forAll(d, celli)
        {
            if (ds[celli] < d[celli])
            {
                d[celli] = ds[celli];
                n[celli] = ns[celli];
                nearestSolid_[celli] = solidi;
            }
        }
    }

    // Distance to the solids alone (cut faces: the bed is classified by
    // its volume fraction instead)
    const scalarField dSolids(d);

    // The bed: distance from the bed model (nearest-solid index -1, i.e.
    // a bed at rest)
    if (bed_)
    {
        bed_->cellDistance(ds, ns);

        forAll(d, celli)
        {
            if (ds[celli] < d[celli])
            {
                d[celli] = ds[celli];
                n[celli] = ns[celli];
                nearestSolid_[celli] = -1;
            }
        }
    }

    distance_.correctBoundaryConditions();
    normal_.correctBoundaryConditions();

    const labelList oldStatus(status_);

    // Solid cells
    status_.resize(m.nCells());

    // A relative tolerance of 1e-6 h makes centres that lie on a surface
    // (to round-off) consistently solid, e.g. a plane bed through a layer
    // of cell centres
    const scalar tol = max(solidThreshold_, 1e-6);

    const bool cut = (cutFaces_ && bed_);

    if (cut)
    {
        computeCutGeometry();
    }

    const scalarField& theta = theta_.primitiveField();

    forAll(status_, celli)
    {
        const scalar threshold = tol*cellExtent(celli, n[celli]);

        if (cut)
        {
            status_[celli] =
            (
                theta[celli] < cutThetaMin_ || dSolids[celli] <= threshold
              ? SOLID
              : LIVE
            );
        }
        else
        {
            status_[celli] = (d[celli] <= threshold ? SOLID : LIVE);
        }
    }

    removeFluidPockets();

    // IB cells: fluid cells with a solid face neighbour, including
    // neighbours across processor boundaries
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();

    labelList nbrStatus;
    syncTools::swapBoundaryCellList(m, status_, nbrStatus);

    labelList newStatus(status_);

    forAll(nei, facei)
    {
        const label o = own[facei];
        const label nb = nei[facei];

        if (status_[o] == SOLID && status_[nb] != SOLID)
        {
            newStatus[nb] = IB;
        }
        else if (status_[nb] == SOLID && status_[o] != SOLID)
        {
            newStatus[o] = IB;
        }
    }

    for (const polyPatch& pp : m.boundaryMesh())
    {
        if (!pp.coupled())
        {
            continue;
        }

        const labelUList& fc = pp.faceCells();

        forAll(fc, i)
        {
            const label bFacei = pp.start() + i - m.nInternalFaces();

            if (nbrStatus[bFacei] == SOLID && status_[fc[i]] != SOLID)
            {
                newStatus[fc[i]] = IB;
            }
        }
    }

    // Cut faces: every fluid cell cut by the bed carries a wall element
    if (cut)
    {
        forAll(newStatus, celli)
        {
            if (newStatus[celli] == LIVE && theta[celli] < 1 - 1e-6)
            {
                newStatus[celli] = IB;
            }
        }
    }

    status_.transfer(newStatus);

    // Cell lists and the post-processing indicator
    DynamicList<label> ib;
    DynamicList<label> solid;

    scalarField& g = gamma_.primitiveFieldRef();

    forAll(status_, celli)
    {
        switch (status_[celli])
        {
            case SOLID:
            {
                solid.push_back(celli);
                g[celli] = -1;
                break;
            }
            case IB:
            {
                ib.push_back(celli);
                g[celli] = 0.5;
                break;
            }
            default:
            {
                g[celli] = 1;
            }
        }
    }

    gamma_.correctBoundaryConditions();

    ibCells_.transfer(ib);
    solidCells_.transfer(solid);

    // Cells that were solid and are fluid now (bed eroded, body moved)
    {
        DynamicList<label> fresh;
        DynamicList<label> buried;
        if (oldStatus.size() == status_.size())
        {
            forAll(status_, celli)
            {
                if (oldStatus[celli] == SOLID && status_[celli] != SOLID)
                {
                    fresh.push_back(celli);
                }
                else if (oldStatus[celli] != SOLID && status_[celli] == SOLID)
                {
                    buried.push_back(celli);
                }
            }
        }
        freshCells_.transfer(fresh);
        buriedCells_.transfer(buried);
    }

    makeFaceMask();

    if (cut && (centroidDiffusion_ || centroidPressure_))
    {
        makeCutDeltaRatio();
    }

    if (cut && (skewCorrection_ || faceInterpolation_))
    {
        makeCutFaceData();
    }

    // Wall area of the IB cells: faces to solid cells projected on the
    // wall normal (for a stair-stepped plane the sum is the true area)
    {
        const vectorField& Sf = m.faceAreas();

        labelList ibIndex(m.nCells(), -1);
        forAll(ibCells_, i)
        {
            ibIndex[ibCells_[i]] = i;
        }

        ibWallArea_.resize(ibCells_.size());
        ibWallArea_ = 0;

        forAll(nei, facei)
        {
            const label o = own[facei];
            const label nb = nei[facei];

            if (ibIndex[o] >= 0 && status_[nb] == SOLID)
            {
                ibWallArea_[ibIndex[o]] += mag(Sf[facei] & n[o]);
            }
            else if (ibIndex[nb] >= 0 && status_[o] == SOLID)
            {
                ibWallArea_[ibIndex[nb]] += mag(Sf[facei] & n[nb]);
            }
        }

        labelList nbrStatus2;
        syncTools::swapBoundaryCellList(m, status_, nbrStatus2);

        for (const polyPatch& pp : m.boundaryMesh())
        {
            if (!pp.coupled())
            {
                continue;
            }

            const labelUList& fc = pp.faceCells();

            forAll(fc, i)
            {
                const label facei = pp.start() + i;
                const label bFacei = facei - m.nInternalFaces();

                if (ibIndex[fc[i]] >= 0 && nbrStatus2[bFacei] == SOLID)
                {
                    ibWallArea_[ibIndex[fc[i]]] +=
                        mag(Sf[facei] & n[fc[i]]);
                }
            }
        }

        // Cut faces: the wall element closes the open faces of the cell,
        // A_w n_w = sum_f (1 - alpha_f) S_f (outward). The wall distance is
        // that of the cell centre where the velocity lives (centre in the
        // fluid), or that of the fluid part, y = theta V/(2 A_w), for a cell
        // whose centre lies below the bed
        if (cut)
        {
            const scalarField& mi = faceMask_.primitiveField();
            vectorField wallVec(ibCells_.size(), Zero);

            // Geometric closure from the apertures alone (faces to solid
            // cells keep their open fraction): the reference for the
            // wall-element consistency check
            vectorField geomVec(ibCells_.size(), Zero);

            forAll(nei, facei)
            {
                const vector closed((1 - mi[facei])*Sf[facei]);
                const vector geom((1 - faceAlpha_[facei])*Sf[facei]);

                if (ibIndex[own[facei]] >= 0)
                {
                    wallVec[ibIndex[own[facei]]] += closed;
                    geomVec[ibIndex[own[facei]]] += geom;
                }
                if (ibIndex[nei[facei]] >= 0)
                {
                    wallVec[ibIndex[nei[facei]]] -= closed;
                    geomVec[ibIndex[nei[facei]]] -= geom;
                }
            }

            forAll(m.boundary(), patchi)
            {
                const fvPatch& fp = m.boundary()[patchi];
                const labelUList& fc = fp.faceCells();
                const scalarField& pm = faceMask_.boundaryField()[patchi];
                const vectorField& pSf = fp.Sf();

                forAll(fc, i)
                {
                    if (ibIndex[fc[i]] >= 0)
                    {
                        wallVec[ibIndex[fc[i]]] += (1 - pm[i])*pSf[i];
                        geomVec[ibIndex[fc[i]]] +=
                            (1 - faceAlpha_[fp.start() + i])*pSf[i];
                    }
                }
            }

            const scalarField& V = m.V();
            const pointField& C = m.cellCentres();
            const scalarField dColumn(d);
            const vectorField nColumn(n);
            nFlooredCutDistance_ = 0;
            maxWallClosureError_ = 0;

            forAll(ibCells_, i)
            {
                const label celli = ibCells_[i];
                const scalar Aw = mag(wallVec[i]);

                ibWallArea_[i] = Aw;

                const scalar magSw = mag(cutWallVec_[celli]);

                if
                (
                    centroidWallDistance_
                 && nearestSolid_[celli] == -1
                 && magSw > VSMALL
                )
                {
                    // Bed, cut cell: normal and distance from the wall
                    // element of the cut surface inside the cell. The
                    // normal points into the fluid; the distance is that
                    // of the fluid centroid from the wall-element plane,
                    // exact for a planar cut (the fluid part is convex, so
                    // the distance is positive). A bed kinked inside the
                    // cell can bring the centroid close to or below the
                    // mean plane: fall back to the column-plane distance
                    // of the centroid, bounded below, and count
                    const vector nw(-cutWallVec_[celli]/magSw);

                    // Consistency: for a planar cut the facet sum equals
                    // the closure of the open faces
                    maxWallClosureError_ = max
                    (
                        maxWallClosureError_,
                        mag(cutWallVec_[celli] - geomVec[i])/magSw
                    );
                    const scalar dc
                    (
                        (fluidCentreOffset_[celli] - cutWallOffset_[celli])
                      & nw
                    );
                    const scalar dMin = 0.05*cellExtent(celli, nw);

                    n[celli] = nw;
                    d[celli] = dc;

                    if (dc < dMin)
                    {
                        d[celli] = max
                        (
                            dColumn[celli]
                          + (fluidCentreOffset_[celli] & nColumn[celli]),
                            dMin
                        );
                        ++nFlooredCutDistance_;

                        if (debug)
                        {
                            Pout<< "immersedBoundary: floored cut cell "
                                << celli << " C " << C[celli]
                                << " theta " << theta[celli]
                                << " d_c " << dc << " dMin " << dMin
                                << " |Sw| " << magSw << " Aw " << Aw
                                << " nw " << nw
                                << " xc*-C " << fluidCentreOffset_[celli]
                                << " xw-C " << cutWallOffset_[celli]
                                << " columnDist " << dColumn[celli]
                                << " used " << d[celli] << endl;
                        }
                    }
                }
                else if (centroidWallDistance_ && nearestSolid_[celli] == -1)
                {
                    // Bed, IB cell not cut: distance of the cell centre to
                    // the column bed plane, positive by classification
                    d[celli] = max
                    (
                        d[celli] + (fluidCentreOffset_[celli] & n[celli]),
                        1e-6*cbrt(V[celli])
                    );
                }
                else if (d[celli] <= 0 && Aw > VSMALL)
                {
                    d[celli] = 0.5*theta[celli]*V[celli]/Aw;
                }
            }

            reduce(nFlooredCutDistance_, sumOp<label>());
            reduce(maxWallClosureError_, maxOp<scalar>());

            distance_.correctBoundaryConditions();
            normal_.correctBoundaryConditions();
        }
    }
}


void Foam::immersedBoundary::computeCutGeometry()
{
    const fvMesh& m = mesh();

    // Point levels at the round-off level of the bed (|phi| ~ 1e-18 m on
    // flat faces) are set to zero, i.e. solid, so that the face and cell
    // fractions do not depend on the sign of round-off
    scalarField phiP(bed_->pointLevel());

    scalar hMin = GREAT;
    forAll(m.V(), celli)
    {
        hMin = min(hMin, cbrt(m.V()[celli]));
    }
    reduce(hMin, minOp<scalar>());

    const scalar phiTol = 1e-9*hMin;

    forAll(phiP, pointi)
    {
        if (mag(phiP[pointi]) <= phiTol)
        {
            phiP[pointi] = 0;
        }
    }

    scalarField& theta = theta_.primitiveFieldRef();
    forAll(theta, celli)
    {
        theta[celli] = cutGeometry::cellFraction(m, celli, phiP);
    }
    theta_.correctBoundaryConditions();

    // Fluid centroids of the cut cells
    fluidCentreOffset_.resize(m.nCells());
    fluidCentreOffset_ = Zero;

    cutWallVec_.resize(m.nCells());
    cutWallVec_ = Zero;
    cutWallOffset_.resize(m.nCells());
    cutWallOffset_ = Zero;

    if
    (
        centroidDiffusion_ || centroidPressure_ || centroidWallDistance_
     || skewCorrection_ || faceInterpolation_
    )
    {
        const pointField& C = m.cellCentres();

        forAll(theta, celli)
        {
            if (theta[celli] > 0 && theta[celli] < 1)
            {
                fluidCentreOffset_[celli] =
                    cutGeometry::cellFluidCentre(m, celli, phiP) - C[celli];

                // Wall element of the same cut surface: the bed normal and
                // the wall distance of the cut cell come from it, so that
                // fraction, apertures, centroid, normal and distance all
                // describe one surface
                vector Sw;
                point xw;
                if (cutGeometry::cellWallElement(m, celli, phiP, Sw, xw))
                {
                    cutWallVec_[celli] = Sw;
                    cutWallOffset_[celli] = xw - C[celli];
                }
            }
        }
    }

    faceAlpha_.resize(m.nFaces());
    forAll(faceAlpha_, facei)
    {
        faceAlpha_[facei] = cutGeometry::faceFraction(m, facei, phiP);
    }

    // Centroids of the open parts of the cut faces (open-face interpolation)
    faceOpenOffset_.resize(m.nFaces());
    faceOpenOffset_ = Zero;

    if (skewCorrection_ || faceInterpolation_)
    {
        const vectorField& Cf = m.faceCentres();

        forAll(faceAlpha_, facei)
        {
            if (faceAlpha_[facei] > 0 && faceAlpha_[facei] < 1)
            {
                faceOpenOffset_[facei] =
                    cutGeometry::faceOpenCentre(m, facei, phiP) - Cf[facei];
            }
        }
    }
}


void Foam::immersedBoundary::removeFluidPockets()
{
    // Fluid cells enclosed by solid cells form a pressure block without a
    // reference (pure Neumann); with an incompatible source the pressure
    // solver stalls. Find the connected fluid regions and turn small
    // isolated ones into solid cells. With cut faces a face whose open
    // fraction is zero is closed as well: a cut cell whose fluid part
    // touches only solid cells (e.g. the bed scoured under a structure
    // member) is then a region of its own, although its cell neighbours
    // are fluid. Left alone, its pressure row is empty and the solver
    // divides by zero.
    const fvMesh& m = mesh();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();

    labelList nbrStatus;
    syncTools::swapBoundaryCellList(m, status_, nbrStatus);

    const bool cut = (cutFaces_ && faceAlpha_.size() == m.nFaces());
    auto closed = [&](const label facei) -> bool
    {
        return cut && faceAlpha_[facei] <= 0;
    };

    // Block the faces between fluid and solid cells and the closed faces
    bitSet blocked(m.nFaces());

    forAll(nei, facei)
    {
        if
        (
            (status_[own[facei]] == SOLID) != (status_[nei[facei]] == SOLID)
         || closed(facei)
        )
        {
            blocked.set(facei);
        }
    }

    for (const polyPatch& pp : m.boundaryMesh())
    {
        if (!pp.coupled())
        {
            continue;
        }

        const labelUList& fc = pp.faceCells();

        forAll(fc, i)
        {
            const label facei = pp.start() + i;
            const label bFacei = facei - m.nInternalFaces();

            if
            (
                (status_[fc[i]] == SOLID) != (nbrStatus[bFacei] == SOLID)
             || closed(facei)
            )
            {
                blocked.set(facei);
            }
        }
    }

    const regionSplit regions(m, blocked);

    // Fluid cells per region
    labelList nFluid(regions.nRegions(), Zero);

    forAll(status_, celli)
    {
        if (status_[celli] != SOLID)
        {
            ++nFluid[regions[celli]];
        }
    }

    Pstream::listCombineReduce(nFluid, plusEqOp<label>());

    const label nFluidTotal = sum(nFluid);
    const label nKeep = minPocketFraction_*nFluidTotal;

    label nPockets = 0;
    label nPocketCells = 0;

    forAll(nFluid, regioni)
    {
        if (nFluid[regioni] > 0 && nFluid[regioni] < nKeep)
        {
            ++nPockets;
            nPocketCells += nFluid[regioni];
        }
    }

    if (nPockets)
    {
        forAll(status_, celli)
        {
            const label n = nFluid[regions[celli]];

            if (status_[celli] != SOLID && n < nKeep)
            {
                status_[celli] = SOLID;
            }
        }

        Info<< "Immersed boundary: " << nPockets
            << " isolated fluid pocket(s) with " << nPocketCells
            << " cells set to solid" << endl;
    }
}


void Foam::immersedBoundary::makeFaceMask()
{
    const fvMesh& m = mesh();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();

    scalarField& mi = faceMask_.primitiveFieldRef();

    scalarField& li = laplacianMask_.primitiveFieldRef();

    // Open fraction of a face between two fluid cells: 1, or the
    // geometric aperture with cut faces
    const bool cut = (cutFaces_ && faceAlpha_.size() == m.nFaces());
    auto open = [&](const label facei) -> scalar
    {
        return (cut ? faceAlpha_[facei] : scalar(1));
    };

    forAll(nei, facei)
    {
        const bool solidO = (status_[own[facei]] == SOLID);
        const bool solidN = (status_[nei[facei]] == SOLID);

        mi[facei] = (!solidO && !solidN ? open(facei) : 0);
        li[facei] = (solidO == solidN ? (solidO ? 1 : open(facei)) : 0);
    }

    labelList nbrStatus;
    syncTools::swapBoundaryCellList(m, status_, nbrStatus);

    auto& mb = faceMask_.boundaryFieldRef();
    auto& lb = laplacianMask_.boundaryFieldRef();

    forAll(mb, patchi)
    {
        const fvPatch& p = m.boundary()[patchi];
        const labelUList& fc = p.faceCells();
        fvsPatchScalarField& pm = mb[patchi];
        fvsPatchScalarField& pl = lb[patchi];

        forAll(pm, i)
        {
            const bool solidO = (status_[fc[i]] == SOLID);

            const label facei = p.start() + i;

            if (p.coupled())
            {
                const label bFacei = facei - m.nInternalFaces();
                const bool solidN = (nbrStatus[bFacei] == SOLID);

                pm[i] = (!solidO && !solidN ? open(facei) : 0);
                pl[i] = (solidO == solidN ? (solidO ? 1 : open(facei)) : 0);
            }
            else
            {
                pm[i] = (solidO ? 0 : open(facei));
                pl[i] = pm[i];
            }
        }
    }
}


void Foam::immersedBoundary::makeCutDeltaRatio()
{
    const fvMesh& m = mesh();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const pointField& C = m.cellCentres();
    const vectorField& Sf = m.faceAreas();
    const scalarField& magSf = m.magSf();

    // Ratio d/d_c of the cell-centre and fluid-centroid normal distances
    // across a face, d_c >= minDistanceFraction*d
    auto ratio = [&]
    (
        const vector& delta,
        const vector& offsetDelta,
        const vector& nf
    ) -> scalar
    {
        const scalar d = (delta & nf);
        const scalar dc = d + (offsetDelta & nf);
        return d/max(dc, minDistanceFraction_*d);
    };

    scalarField& ri = cutDeltaRatio_.primitiveFieldRef();
    ri = 1;

    forAll(nei, facei)
    {
        const label o = own[facei];
        const label nb = nei[facei];

        if (status_[o] == SOLID || status_[nb] == SOLID)
        {
            continue;
        }

        const vector dOff(fluidCentreOffset_[nb] - fluidCentreOffset_[o]);

        if (magSqr(dOff) > 0)
        {
            ri[facei] =
                ratio(C[nb] - C[o], dOff, Sf[facei]/magSf[facei]);
        }
    }

    labelList nbrStatus;
    syncTools::swapBoundaryCellList(m, status_, nbrStatus);

    pointField nbrC;
    syncTools::swapBoundaryCellPositions(m, C, nbrC);

    vectorField nbrOffset;
    syncTools::swapBoundaryCellList(m, fluidCentreOffset_, nbrOffset);

    auto& rb = cutDeltaRatio_.boundaryFieldRef();

    forAll(rb, patchi)
    {
        const fvPatch& p = m.boundary()[patchi];
        fvsPatchScalarField& pr = rb[patchi];
        pr = 1;

        if (!p.coupled())
        {
            continue;
        }

        const labelUList& fc = p.faceCells();

        forAll(pr, i)
        {
            const label facei = p.start() + i;
            const label bFacei = facei - m.nInternalFaces();

            if (status_[fc[i]] == SOLID || nbrStatus[bFacei] == SOLID)
            {
                continue;
            }

            const vector dOff(nbrOffset[bFacei] - fluidCentreOffset_[fc[i]]);

            if (magSqr(dOff) > 0)
            {
                // magSf() has internal faces only: use the face area vector
                pr[i] =
                    ratio
                    (
                        nbrC[bFacei] - C[fc[i]],
                        dOff,
                        Sf[facei]/mag(Sf[facei])
                    );
            }
        }
    }

    Info<< "    cut-face distance ratio d/d_c: min "
        << gMin(ri) << ", max " << gMax(ri) << endl;
}


void Foam::immersedBoundary::report() const
{
    Info<< "Immersed boundary: "
        << returnReduce(ibCells_.size(), sumOp<label>()) << " IB cells, "
        << returnReduce(solidCells_.size(), sumOp<label>()) << " solid cells, "
        << returnReduce
           (
               status_.size() - ibCells_.size() - solidCells_.size(),
               sumOp<label>()
           )
        << " live cells" << endl;

    if (cutFaces_ && bed_)
    {
        label nCut = 0;
        forAll(status_, celli)
        {
            if (status_[celli] != SOLID && theta_[celli] < 1 - SMALL)
            {
                ++nCut;
            }
        }
        Info<< "    cut faces: " << returnReduce(nCut, sumOp<label>())
            << " fluid cells cut by the bed (theta_min " << cutThetaMin_
            << "), " << nFlooredCutDistance_
            << " with the wall-element distance floored, max closure error "
            << maxWallClosureError_ << ", "
            << returnReduce(nTruncatedGrad_, sumOp<label>())
            << " with an unresolved gradient direction" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundary::immersedBoundary(const fvMesh& mesh)
:
    MeshObject_type(mesh),
    dict_
    (
        IOobject
        (
            "immersedBoundaryDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    solids_(),
    bed_(nullptr),
    bedMoved_(false),
    freshCells_(),
    solidThreshold_(dict_.getOrDefault<scalar>("solidThreshold", 0)),
    minPocketFraction_
    (
        dict_.getOrDefault<scalar>("minFluidRegionFraction", 0.01)
    ),
    penaltyFactor_(dict_.getOrDefault<scalar>("solidPenaltyFactor", 1e3)),
    status_(),
    nearestSolid_(mesh.nCells(), -1),
    ibCells_(),
    solidCells_(),
    distance_
    (
        IOobject
        (
            "ibDistance",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            dict_.getOrDefault<bool>("writeFields", true)
          ? IOobject::AUTO_WRITE
          : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimLength, GREAT)
    ),
    meshWallDistance_(),
    normal_
    (
        IOobject
        (
            "ibNormal",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimless, vector(0, 0, 1))
    ),
    gamma_
    (
        IOobject
        (
            "Gamma",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            dict_.getOrDefault<bool>("writeFields", true)
          ? IOobject::AUTO_WRITE
          : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 1)
    ),
    faceMask_
    (
        IOobject
        (
            "ibFaceMask",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 1)
    ),
    laplacianMask_
    (
        IOobject
        (
            "ibLaplacianMask",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 1)
    ),
    wallModel_(nullptr),
    pressureSubset_
    (
        dict_.getOrDefault<word>("pressureSolve", "subset") == "subset"
    ),
    fluidSubset_(nullptr),
    subsetCell_(),
    cutFaces_(dict_.getOrDefault<bool>("cutFaces", false)),
    cutThetaMin_(dict_.getOrDefault<scalar>("cutThetaMin", 0.05)),
    theta_
    (
        IOobject
        (
            "ibTheta",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            dict_.getOrDefault<bool>("writeFields", true)
         && dict_.getOrDefault<bool>("cutFaces", false)
          ? IOobject::AUTO_WRITE
          : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 1)
    ),
    faceAlpha_(),
    centroidDiffusion_
    (
        dict_.subOrEmptyDict("cutCentroids")
       .getOrDefault<bool>("diffusion", false)
    ),
    centroidPressure_
    (
        dict_.subOrEmptyDict("cutCentroids")
       .getOrDefault<bool>("pressure", false)
    ),
    centroidWallDistance_
    (
        dict_.subOrEmptyDict("cutCentroids")
       .getOrDefault<bool>("wallDistance", false)
    ),
    minDistanceFraction_
    (
        dict_.subOrEmptyDict("cutCentroids")
       .getOrDefault<scalar>("minDistanceFraction", 0.5)
    ),
    fluidCentreOffset_(mesh.nCells(), Zero),
    nFlooredCutDistance_(0),
    nTruncatedGrad_(0),
    maxWallClosureError_(0),
    skewCorrection_
    (
        dict_.subOrEmptyDict("cutCentroids")
       .getOrDefault<bool>("skewCorrection", true)
    ),
    faceInterpolation_
    (
        dict_.subOrEmptyDict("cutCentroids")
       .getOrDefault<bool>("faceInterpolation", true)
    ),
    faceOpenOffset_(),
    cutFaceFlag_
    (
        IOobject
        (
            "ibCutFaceFlag",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    ),
    cutWeights_
    (
        IOobject
        (
            "ibCutWeights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh.weights()
    ),
    cutInterpOffset_
    (
        IOobject
        (
            "ibCutInterpOffset",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimLength, Zero)
    ),
    cutSkewVec_
    (
        IOobject
        (
            "ibCutSkewVec",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimless, Zero)
    ),
    cutDeltaRatio_
    (
        IOobject
        (
            "ibCutDeltaRatio",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 1)
    ),
    gradMode_(0),
    tangentialSink_(false),
    harmonicRAU_(false),
    maskNonOrthCorr_(dict_.getOrDefault<bool>("maskNonOrthCorrection", true))
{
    Info<< "Constructing immersed boundary" << endl;

    {
        const word mode(dict_.getOrDefault<word>("pressureGradient", "gauss"));

        if (mode == "faceIB")
        {
            gradMode_ = 1;
        }
        else if (mode == "face")
        {
            gradMode_ = 2;
        }
        else if (mode != "gauss")
        {
            FatalIOErrorInFunction(dict_)
                << "Unknown pressureGradient " << mode
                << ", valid: gauss faceIB face" << exit(FatalIOError);
        }

        Info<< "    pressure gradient: " << mode << endl;
    }

    {
        const word sink(dict_.getOrDefault<word>("wallSink", "isotropic"));

        if (sink != "isotropic" && sink != "tangential")
        {
            FatalIOErrorInFunction(dict_)
                << "Unknown wallSink " << sink
                << ", valid: isotropic tangential" << exit(FatalIOError);
        }

        tangentialSink_ = (sink == "tangential");

        Info<< "    wall sink: " << sink << endl;
    }

    Info<< "    non-orthogonal correction on the faces of IB and solid cells: "
        << (maskNonOrthCorr_ ? "removed" : "kept") << endl;

    {
        const word interp
        (
            dict_.getOrDefault<word>("rAUInterpolation", "linear")
        );

        if (interp != "linear" && interp != "harmonic")
        {
            FatalIOErrorInFunction(dict_)
                << "Unknown rAUInterpolation " << interp
                << ", valid: linear harmonic" << exit(FatalIOError);
        }

        harmonicRAU_ = (interp == "harmonic");

        Info<< "    rAU interpolation: " << interp << endl;
    }

    if (centroidDiffusion_ || centroidPressure_ || centroidWallDistance_)
    {
        Info<< "    cut cells at the fluid centroids:"
            << (centroidDiffusion_ ? " diffusion" : "")
            << (centroidPressure_ ? " pressure" : "")
            << (centroidWallDistance_ ? " wallDistance" : "")
            << ", minDistanceFraction " << minDistanceFraction_ << endl;
    }

    if (skewCorrection_ || faceInterpolation_)
    {
        Info<< "    cut faces:"
            << (skewCorrection_ ? " skew-corrected diffusion flux" : "")
            << (faceInterpolation_ ? " interpolation to the open-face centroid" : "")
            << endl;
    }

    readSolids();
    classify();
    report();

    wallModel_.reset(new ibWallModel(*this, dict_.subDict("wallModel")));
    wallModel_->updateGeometry();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundary::~immersedBoundary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::ibWallModel& Foam::immersedBoundary::wallModel() const
{
    return *wallModel_;
}


Foam::ibWallModel& Foam::immersedBoundary::wallModel()
{
    return *wallModel_;
}


bool Foam::immersedBoundary::moving() const
{
    if (bedMoved_)
    {
        return true;
    }

    for (const ibSolid& s : solids_)
    {
        if (s.moving())
        {
            return true;
        }
    }

    return false;
}


Foam::scalar Foam::immersedBoundary::cellExtent
(
    const label celli,
    const vector& n
) const
{
    const fvMesh& m = mesh();
    const vectorField& Sf = m.faceAreas();

    scalar sumProj = 0;

    for (const label facei : m.cells()[celli])
    {
        sumProj += mag(n & Sf[facei]);
    }

    return 2*m.V()[celli]/max(sumProj, VSMALL);
}


Foam::tmp<Foam::vectorField> Foam::immersedBoundary::solidVelocity
(
    const labelUList& cells,
    const pointField& points
) const
{
    auto tU = tmp<vectorField>::New(cells.size(), Zero);
    vectorField& U = tU.ref();

    forAll(cells, i)
    {
        const label solidi = nearestSolid_[cells[i]];

        if (solidi >= 0 && solids_[solidi].moving())
        {
            U[i] = solids_[solidi].wallVelocity(pointField(1, points[i]))()[0];
        }
    }

    return tU;
}


Foam::vector Foam::immersedBoundary::solidForce
(
    const label solidi,
    const scalarField& pIB,
    const scalarField& rhoIB
) const
{
    vector F(Zero);

    const vectorField& n = normal().primitiveField();
    const vectorField& tw = wallModel_->wallShearStress().primitiveField();
    const scalarField& Aw = ibWallArea();

    forAll(ibCells_, i)
    {
        const label celli = ibCells_[i];
        if (nearestSolid_[celli] == solidi)
        {
            F += (-pIB[i]*n[celli] + rhoIB[i]*tw[celli])*Aw[i];
        }
    }

    reduce(F, sumOp<vector>());

    return F;
}


Foam::vector Foam::immersedBoundary::solidForce
(
    const label solidi,
    const scalar rho
) const
{
    if (!mesh().foundObject<volScalarField>("p"))
    {
        return vector(Zero);
    }

    const volScalarField& p = mesh().lookupObject<volScalarField>("p");
    scalarField pIB(ibCells_.size());
    forAll(ibCells_, i)
    {
        pIB[i] = rho*p[ibCells_[i]];
    }

    return solidForce(solidi, pIB, scalarField(ibCells_.size(), rho));
}


Foam::vector Foam::immersedBoundary::solidMoment
(
    const label solidi,
    const scalarField& pIB,
    const scalarField& rhoIB,
    const point& xRef
) const
{
    vector M(Zero);

    const vectorField& n = normal().primitiveField();
    const scalarField& d = distance().primitiveField();
    const vectorField& tw = wallModel_->wallShearStress().primitiveField();
    const scalarField& Aw = ibWallArea();
    const pointField& C = mesh().cellCentres();

    forAll(ibCells_, i)
    {
        const label celli = ibCells_[i];
        if (nearestSolid_[celli] == solidi)
        {
            const point xw(C[celli] - d[celli]*n[celli]);
            M += (xw - xRef) ^ ((-pIB[i]*n[celli] + rhoIB[i]*tw[celli])*Aw[i]);
        }
    }

    reduce(M, sumOp<vector>());

    return M;
}


Foam::vector Foam::immersedBoundary::solidMoment
(
    const label solidi,
    const scalar rho,
    const point& xRef
) const
{
    if (!mesh().foundObject<volScalarField>("p"))
    {
        return vector(Zero);
    }

    const volScalarField& p = mesh().lookupObject<volScalarField>("p");
    scalarField pIB(ibCells_.size());
    forAll(ibCells_, i)
    {
        pIB[i] = rho*p[ibCells_[i]];
    }

    return solidMoment(solidi, pIB, scalarField(ibCells_.size(), rho), xRef);
}


void Foam::immersedBoundary::update()
{
    addProfiling(updateProf, "immersedBoundary::update");

    // The fresh and buried cells belong to the step that re-classified;
    // a step without a re-classification has none (stale lists made
    // initialiseFreshCells overwrite the same fresh cells and move the
    // ghost water of the same buried cells up at every step: water
    // created and the flow reset at the bed, dam break 2026-11-02)
    freshCells_.clear();
    buriedCells_.clear();

    // Solids with a prescribed displacement (e.g. a lifted gate) move to
    // their position at the new time; a changed position calls for a new
    // classification even if the solid has come to rest
    bool prescribedMoved = false;
    {
        const scalar t = mesh_.time().value();
        const scalar dt = mesh_.time().deltaTValue();
        for (ibSolid& s : solids_)
        {
            if (s.prescribed() && s.applyPrescribedMotion(t, dt))
            {
                prescribedMoved = true;
            }
        }
    }

    if (!moving() && !prescribedMoved)
    {
        return;
    }

    classify();
    report();
    wallModel_->updateGeometry();
    correctWallDistance();

    // The fluid-only subset follows the classification
    fluidSubset_.reset(nullptr);

    if (bed_)
    {
        bed_->classified();
        if (prescribedMoved)
        {
            bed_->solidsMoved();
        }
    }
    bedMoved_ = false;

    Info<< "Immersed boundary: "
        << returnReduce(freshCells_.size(), sumOp<label>())
        << " fresh cells" << endl;
}


void Foam::immersedBoundary::evolveBed(const scalar dt)
{
    addProfiling(evolveBedProf, "immersedBoundary::evolveBed");

    if (bed_)
    {
        bedMoved_ = bed_->evolve(*this, dt) || bedMoved_;

        // Solids the bed moves (e.g. liquefied pipes): advance their
        // reference points with the prescribed velocity; any displacement
        // calls for a new classification at the next update (also when
        // the solid comes to rest, with zero velocity)
        HashTable<vector> velocities;
        HashTable<vector> omegas;
        bed_->solidVelocities(velocities);
        bed_->solidAngularVelocities(omegas);
        for (ibSolid& s : solids_)
        {
            if (velocities.found(s.name()) || omegas.found(s.name()))
            {
                const vector v
                (
                    velocities.found(s.name())
                  ? velocities[s.name()]
                  : vector(Zero)
                );
                const vector w
                (
                    omegas.found(s.name()) ? omegas[s.name()] : vector(Zero)
                );
                if (magSqr(v) > 0 || magSqr(w) > 0)
                {
                    bedMoved_ = true;
                }
                tensor R(s.R());
                const scalar angle = dt*mag(w);
                if (angle > VSMALL)
                {
                    R = quaternion(w/mag(w), angle).R() & R;
                }
                s.setState(s.origin() + dt*v, R, v, w);
            }
        }
    }
}


void Foam::immersedBoundary::initialiseFreshCells(volVectorField& U) const
{
    if (returnReduce(freshCells_.empty(), andOp<bool>()))
    {
        return;
    }

    // Fresh IB cells start at the wall-law velocity, the others at the
    // wall velocity; their faces start with zero flux and the next
    // pressure solve makes the flux divergence-free
    wallModel_->initialiseFreshCells(freshCells_, U);
}


void Foam::immersedBoundary::initialiseFreshCells
(
    volVectorField& U,
    volScalarField& alpha,
    volScalarField& p_rgh
) const
{
    initialiseFreshCells(U);

    const bool anyFresh = !returnReduce(freshCells_.empty(), andOp<bool>());
    const bool anyBuried = !returnReduce(buriedCells_.empty(), andOp<bool>());
    if (!anyFresh && !anyBuried)
    {
        return;
    }

    const fvMesh& m = mesh();
    const cellList& cells = m.cells();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const vectorField& Cc = m.cellCentres();
    const scalarField& V = m.V();
    const label nInternal = m.nInternalFaces();
    scalarField& a = alpha.primitiveFieldRef();
    scalarField& pr = p_rgh.primitiveFieldRef();

    // The fluid face neighbour with the highest centre (the cell above),
    // or -1 (processor neighbours are not considered: a bed moves by a
    // fraction of a cell per step and the cell above is almost always
    // local)
    auto fluidAbove = [&](const label celli) -> label
    {
        label best = -1;
        scalar zBest = -GREAT;
        for (const label facei : cells[celli])
        {
            if (facei >= nInternal)
            {
                continue;
            }
            const label other = (own[facei] == celli ? nei[facei] : own[facei]);
            if (status_[other] != SOLID && Cc[other].z() > zBest)
            {
                zBest = Cc[other].z();
                best = other;
            }
        }
        return best;
    };

    // Fresh cells: volume fraction and pressure of the cell above (or of
    // any fluid neighbour). The water filling a fresh cell is drawn from
    // the top of the water column above it (an eroding bed lowers the free
    // surface); without the draw every cell the bed uncovers created water,
    // and a cell flipping between buried and fresh created it each cycle
    label nNoNeighbour = 0;
    scalar drawn = 0;
    scalar shortfall = 0;
    for (const label celli : freshCells_)
    {
        const label above = fluidAbove(celli);
        if (above < 0)
        {
            ++nNoNeighbour;
            continue;
        }
        a[celli] = a[above];
        pr[celli] = pr[above];

        scalar need = a[celli]*V[celli];
        if (need <= 0)
        {
            continue;
        }

        // The column above, bottom-up
        DynamicList<label> column;
        label cur = celli;
        for (label climb = 0; climb < 1000; ++climb)
        {
            const label up = fluidAbove(cur);
            if (up < 0 || Cc[up].z() <= Cc[cur].z())
            {
                break;
            }
            column.push_back(up);
            cur = up;
        }

        // Take the water from the top of the column down
        for (label k = column.size() - 1; k >= 0 && need > 0; --k)
        {
            const label c = column[k];
            const scalar take = min(max(a[c], scalar(0))*V[c], need);
            a[c] -= take/V[c];
            need -= take;
            drawn += take;
        }
        shortfall += need;
    }

    // Buried cells: their water goes up the column to the first cells with
    // room (a bed rising under water displaces it to the free surface)
    scalar moved = 0;
    scalar lost = 0;
    for (const label celli : buriedCells_)
    {
        scalar water = a[celli]*V[celli];
        if (water <= 0)
        {
            continue;
        }
        label cur = celli;
        for (label climb = 0; climb < 1000 && water > 0; ++climb)
        {
            const label above = fluidAbove(cur);
            if (above < 0 || Cc[above].z() <= Cc[cur].z())
            {
                break;
            }
            const scalar room = max(scalar(1) - a[above], scalar(0))*V[above];
            const scalar take = min(water, room);
            a[above] += take/V[above];
            moved += take;
            water -= take;
            cur = above;
        }
        lost += water;
    }

    reduce(moved, sumOp<scalar>());
    reduce(lost, sumOp<scalar>());
    reduce(drawn, sumOp<scalar>());
    reduce(shortfall, sumOp<scalar>());
    reduce(nNoNeighbour, sumOp<label>());
    Info<< "Immersed boundary: " << returnReduce(freshCells_.size(), sumOp<label>())
        << " fresh cells took the volume fraction of the cell above ("
        << nNoNeighbour << " without a fluid neighbour; " << drawn
        << " m3 drawn from the column above, " << shortfall
        << " m3 not found), "
        << returnReduce(buriedCells_.size(), sumOp<label>())
        << " buried cells moved " << moved << " m3 of water up, lost "
        << lost << " m3" << endl;

    // The boundary of alpha is evaluated here; that of p_rgh is left to
    // the pressure solve (a fixedFluxPressure patch cannot be evaluated
    // before the gradient of the step is set, pier run 2026-10-29)
    alpha.correctBoundaryConditions();
}


Foam::scalar Foam::immersedBoundary::fluidVolume
(
    const volScalarField& alpha
) const
{
    const scalarField& V = mesh().V();
    scalar vol = 0;
    forAll(status_, celli)
    {
        if (status_[celli] != SOLID)
        {
            vol += alpha[celli]*V[celli];
        }
    }
    return returnReduce(vol, sumOp<scalar>());
}


void Foam::immersedBoundary::fluidOnly
(
    const volScalarField& alpha,
    volScalarField& result
) const
{
    result.primitiveFieldRef() = alpha.primitiveField();
    scalarField& r = result.primitiveFieldRef();
    for (const label celli : solidCells_)
    {
        r[celli] = 0;
    }
    result.correctBoundaryConditions();
}


void Foam::immersedBoundary::correctWallDistance() const
{
    // The wall distance is a const MeshObject member of the turbulence
    // models; it is computed from wall patches only and does not know the
    // immersed solids. Limit it by the distance to the solids.
    if (!mesh().foundObject<wallDist>(wallDist::typeName))
    {
        return;
    }

    volScalarField& y =
        const_cast<volScalarField&>(wallDist::New(mesh()).y());

    scalarField& yi = y.primitiveFieldRef();
    const scalarField& d = distance_.primitiveField();

    // The mesh-wall distance is taken once, before any limiting, so that
    // the limit can be rebuilt (and increase) when the bed erodes
    if (meshWallDistance_.size() != yi.size())
    {
        meshWallDistance_ = yi;
    }

    // meshWave leaves the cells no wall patch reaches at -GREAT (a mesh
    // without wall patches has them everywhere): no mesh wall there, so
    // the immersed distance is the wall distance. Taking min() with the
    // invalid value gave y = SMALL in the whole domain, i.e. the SST
    // blending F1 = F2 = 1 (Bradshaw limiter) in every free shear layer.
    label nNoWall = 0;

    forAll(yi, celli)
    {
        const scalar yMesh = meshWallDistance_[celli];

        if (yMesh > 0)
        {
            yi[celli] = min(yMesh, d[celli]);
        }
        else
        {
            yi[celli] = d[celli];
            ++nNoWall;
        }

        yi[celli] = max(yi[celli], SMALL);
    }

    y.correctBoundaryConditions();

    Info<< "Immersed boundary: wall distance limited by the solids, "
        << returnReduce(nNoWall, sumOp<label>())
        << " cells without a mesh wall; y min "
        << gMin(yi) << " max " << gMax(yi) << endl;
}


void Foam::immersedBoundary::maskFlux(surfaceScalarField& phi) const
{
    phi *= faceMask_;
}


Foam::tmp<Foam::surfaceScalarField> Foam::immersedBoundary::ibTouchMask() const
{
    const fvMesh& m = mesh();
    auto tmask = surfaceScalarField::New
    (
        "ibTouchMask",
        IOobject::NO_REGISTER,
        m,
        dimensionedScalar(dimless, Zero)
    );
    surfaceScalarField& mask = tmask.ref();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    scalarField& mi = mask.primitiveFieldRef();
    forAll(nei, facei)
    {
        if (status_[own[facei]] != LIVE || status_[nei[facei]] != LIVE)
        {
            mi[facei] = 1;
        }
    }
    labelList nbrStatus;
    syncTools::swapBoundaryCellList(m, status_, nbrStatus);
    auto& mb = mask.boundaryFieldRef();
    forAll(mb, patchi)
    {
        const fvPatch& p = m.boundary()[patchi];
        if (!p.coupled())
        {
            continue;
        }
        const labelUList& fc = p.faceCells();
        forAll(fc, i)
        {
            const label bFacei = p.start() + i - m.nInternalFaces();
            if (status_[fc[i]] != LIVE || nbrStatus[bFacei] != LIVE)
            {
                mb[patchi][i] = 1;
            }
        }
    }
    return tmask;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::immersedBoundary::maskedLaplacianCoeff(const volScalarField& rAU) const
{
    tmp<surfaceScalarField> trAUf
    (
        faceInterpolation_ ? cutInterpolate(rAU) : fvc::interpolate(rAU)
    );

    if (harmonicRAU_)
    {
        // Faces of IB cells: harmonic mean with the linear weights,
        // 1/rAU_f = w/rAU_P + (1 - w)/rAU_N
        const fvMesh& m = mesh();
        const labelUList& own = m.owner();
        const labelUList& nei = m.neighbour();
        const scalarField& w = m.weights().primitiveField();
        const scalarField& r = rAU.primitiveField();
        scalarField& rf = trAUf.ref().primitiveFieldRef();

        forAll(nei, facei)
        {
            const label o = own[facei];
            const label nb = nei[facei];

            if (status_[o] == IB || status_[nb] == IB)
            {
                rf[facei] =
                    1/(w[facei]/max(r[o], VSMALL)
                  + (1 - w[facei])/max(r[nb], VSMALL));
            }
        }

        labelList nbrStatus;
        syncTools::swapBoundaryCellList(m, status_, nbrStatus);

        auto& rb = trAUf.ref().boundaryFieldRef();

        forAll(rb, patchi)
        {
            const fvPatch& p = m.boundary()[patchi];

            if (!p.coupled())
            {
                continue;
            }

            const scalarField rN(rAU.boundaryField()[patchi].patchNeighbourField());
            const scalarField& pw = p.weights();
            const labelUList& fc = p.faceCells();

            forAll(fc, i)
            {
                const label bFacei = p.start() + i - m.nInternalFaces();

                if (status_[fc[i]] == IB || nbrStatus[bFacei] == IB)
                {
                    rb[patchi][i] =
                        1/(pw[i]/max(r[fc[i]], VSMALL)
                      + (1 - pw[i])/max(rN[i], VSMALL));
                }
            }
        }
    }

    if (centroidPressure_ && cutFaces_)
    {
        return laplacianMask_*cutDeltaRatio_*trAUf;
    }

    return laplacianMask_*trAUf;
}


void Foam::immersedBoundary::reportEmptyRow
(
    const label celli,
    const surfaceScalarField& rAUf
) const
{
    const fvMesh& m = mesh();
    const scalarField& li = laplacianMask_.primitiveField();
    const scalarField& ri = rAUf.primitiveField();

    Pout<< "Empty pressure row: cell " << celli
        << " centre " << m.cellCentres()[celli]
        << " status " << status_[celli]
        << " theta " << theta_[celli]
        << " nearestSolid " << nearestSolid_[celli] << nl;

    for (const label facei : m.cells()[celli])
    {
        const scalar alpha = (faceAlpha_.size() ? faceAlpha_[facei] : -1);

        if (m.isInternalFace(facei))
        {
            const label nb =
            (
                m.owner()[facei] == celli
              ? m.neighbour()[facei]
              : m.owner()[facei]
            );
            Pout<< "    face " << facei << " internal, nbr " << nb
                << " status " << status_[nb] << " theta " << theta_[nb]
                << " alpha " << alpha << " mask " << li[facei]
                << " rAUf " << ri[facei] << " Sf " << m.faceAreas()[facei]
                << nl;
        }
        else
        {
            const label patchi = m.boundaryMesh().whichPatch(facei);
            const label pf = facei - m.boundaryMesh()[patchi].start();
            Pout<< "    face " << facei << " patch "
                << m.boundaryMesh()[patchi].name()
                << " alpha " << alpha
                << " mask " << laplacianMask_.boundaryField()[patchi][pf]
                << " rAUf " << rAUf.boundaryField()[patchi][pf]
                << " Sf " << m.faceAreas()[facei] << nl;
        }
    }
}


void Foam::immersedBoundary::solvePressure
(
    volScalarField& p,
    const surfaceScalarField& rAUf,
    const surfaceScalarField& phiHbyA,
    surfaceScalarField& phi,
    pimpleControl& pimple,
    const label pRefCell,
    const scalar pRefValue
) const
{
    addProfiling(solvePressureProf, "immersedBoundary::solvePressure");

    const fvMesh& m = mesh();

    // Fluid-only subset; the exposed faces (to solid cells) form an empty
    // patch, i.e. zero flux and zero normal gradient of p
    if (!fluidSubset_)
    {
        bitSet fluidCells(m.nCells());
        forAll(status_, celli)
        {
            if (status_[celli] != SOLID)
            {
                fluidCells.set(celli);
            }
        }

        fluidSubset_.reset(new fvMeshSubset(m, fluidCells));

        subsetCell_.resize(m.nCells());
        subsetCell_ = -1;

        const labelList& cellMap = fluidSubset_->cellMap();
        forAll(cellMap, i)
        {
            subsetCell_[cellMap[i]] = i;
        }

        Info<< "Immersed boundary: pressure solved on a fluid subset of "
            << returnReduce(cellMap.size(), sumOp<label>()) << " cells"
            << endl;
    }

    const fvMeshSubset& sub = *fluidSubset_;
    const fvMesh& sm = sub.subMesh();
    const labelList& cellMap = sub.cellMap();
    const labelList& faceMap = sub.faceMap();

    tmp<volScalarField> tpSub(sub.interpolate(p));

    // Boundary conditions of the subset copy: derived types that look up
    // other fields (totalPressure) or need their gradient set each solve
    // (fixedFluxPressure) are replaced by plain fixedValue or
    // fixedGradient conditions with the values and gradients of the
    // full-mesh field mapped face by face (one corrector behind for the
    // dynamic part of a total pressure)
    {
        volScalarField& pSub = tpSub.ref();
        const labelList& patchMap = sub.patchMap();
        const labelList& faceMap = sub.faceMap();
        forAll(pSub.boundaryField(), patchi)
        {
            const label origi = patchMap[patchi];
            if (origi < 0)
            {
                continue;
            }
            const fvPatchScalarField& pOrig = p.boundaryField()[origi];
            const fvPatch& fp = pSub.mesh().boundary()[patchi];
            const label startOrig = m.boundaryMesh()[origi].start();
            const label startSub = fp.start();

            const bool gradType = isA<fixedGradientFvPatchScalarField>(pOrig);
            const bool valueType =
                pOrig.fixesValue() && !pOrig.coupled()
             && pOrig.type() != fixedValueFvPatchScalarField::typeName;
            if (!gradType && !valueType)
            {
                continue;
            }

            scalarField mapped(fp.size());
            if (gradType)
            {
                const scalarField& g =
                    refCast<const fixedGradientFvPatchScalarField>(pOrig)
                   .gradient();
                forAll(mapped, i)
                {
                    mapped[i] = g[faceMap[startSub + i] - startOrig];
                }
                tmp<fvPatchScalarField> tnew
                (
                    fvPatchScalarField::New
                    (
                        fixedGradientFvPatchScalarField::typeName,
                        fp,
                        pSub
                    )
                );
                refCast<fixedGradientFvPatchScalarField>(tnew.ref())
                   .gradient() = mapped;
                tnew.ref() == pSub.boundaryField()[patchi];
                pSub.boundaryFieldRef().set(patchi, tnew);
            }
            else
            {
                forAll(mapped, i)
                {
                    mapped[i] = pOrig[faceMap[startSub + i] - startOrig];
                }
                tmp<fvPatchScalarField> tnew
                (
                    fvPatchScalarField::New
                    (
                        fixedValueFvPatchScalarField::typeName,
                        fp,
                        pSub
                    )
                );
                tnew.ref() == mapped;
                pSub.boundaryFieldRef().set(patchi, tnew);
            }
        }
    }
    volScalarField& pSub = tpSub.ref();

    // The subset mesh has its own fluxRequired list
    const_cast<fvMesh&>(sm).setFluxRequired(pSub.name());

    const surfaceScalarField rAUfSub(sub.interpolate(rAUf));
    const surfaceScalarField phiHbyASub(sub.interpolate(phiHbyA));
    const surfaceScalarField touchSub(sub.interpolate(ibTouchMask()()));

    surfaceScalarField phiSub(phiHbyASub);

    // Reference cell in the subset (only needed without fixed-value p)
    const label refSub =
    (
        pRefCell >= 0 && pRefCell < m.nCells() ? subsetCell_[pRefCell] : -1
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUfSub, pSub) == fvc::div(phiHbyASub)
        );
        maskNonOrthCorrection(pEqn, touchSub);

        // A fluid cell without an open face has an empty row (the
        // classification removes such cells; this guards the solver
        // against a division by zero if one slips through): hold it at
        // the reference pressure. The coefficients of the coupled
        // (processor) faces are not in diag() but are added to it before
        // the solve, so they count as open faces here
        {
            scalarField diagAll(pEqn.diag());
            const FieldField<Field, scalar>& ic = pEqn.internalCoeffs();

            forAll(sm.boundary(), patchi)
            {
                if (sm.boundary()[patchi].coupled())
                {
                    const labelUList& fc = sm.boundary()[patchi].faceCells();

                    forAll(fc, i)
                    {
                        diagAll[fc[i]] += ic[patchi][i];
                    }
                }
            }

            scalarField& diag = pEqn.diag();
            scalarField& source = pEqn.source();
            label nEmpty = 0;

            forAll(diagAll, i)
            {
                if (mag(diagAll[i]) < VSMALL)
                {
                    diag[i] = -1;
                    source[i] = -pRefValue;
                    ++nEmpty;

                    // Report the first few: the cell, its status and
                    // fraction, and every face with its aperture, mask,
                    // coefficient and neighbour status
                    static label nReported = 0;
                    if (nReported < 3)
                    {
                        ++nReported;
                        reportEmptyRow(cellMap[i], rAUf);
                    }
                }
            }

            if (returnReduce(nEmpty, sumOp<label>()))
            {
                WarningInFunction
                    << returnReduce(nEmpty, sumOp<label>())
                    << " fluid cells without an open face in the pressure "
                    << "equation; held at the reference pressure" << endl;
            }
        }

        if (pSub.needReference())
        {
            pEqn.setReference(max(refSub, 0), pRefValue, refSub >= 0);
        }

        pEqn.solve(m.solverDict(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phiSub = phiHbyASub - pEqn.flux();
        }
    }

    // Pressure back to the base mesh; solid cells keep their values (they
    // do not enter the fluid: blocked faces, IB-aware gradient)
    scalarField& pi = p.primitiveFieldRef();
    forAll(cellMap, i)
    {
        pi[cellMap[i]] = pSub[i];
    }
    p.correctBoundaryConditions();

    // Flux: zero on blocked faces, subset values elsewhere
    phi = phiHbyA;
    scalarField& phii = phi.primitiveFieldRef();

    forAll(phiSub, facei)
    {
        phii[faceMap[facei]] = phiSub[facei];
    }

    auto& phib = phi.boundaryFieldRef();

    forAll(sm.boundary(), patchi)
    {
        const fvPatch& sp = sm.boundary()[patchi];
        const fvsPatchScalarField& psub = phiSub.boundaryField()[patchi];

        forAll(psub, i)
        {
            const label baseFacei = faceMap[sp.start() + i];

            // Exposed faces were internal faces of the base mesh: blocked
            if (baseFacei < m.nInternalFaces())
            {
                continue;
            }

            const label basePatchi =
                m.boundaryMesh().whichPatch(baseFacei);
            const label locali =
                baseFacei - m.boundaryMesh()[basePatchi].start();

            phib[basePatchi][locali] = psub[i];
        }
    }

    (void)sm;
}


void Foam::immersedBoundary::fixSolidCells(fvScalarMatrix& pEqn) const
{
    // Solid cells are cut off from the fluid but coupled to each other. A
    // diagonal penalty holds them at a fixed pressure: the solid block is
    // well conditioned and the fluid solution is not affected.
    //
    // The penalty of each solid cell is solidPenaltyFactor times its own
    // Laplacian diagonal, so the balance between penalty and coupling does
    // not depend on the cell size (a global scale fails on meshes with
    // strong refinement: coarse solid cells become isolated rows and the
    // GAMG hierarchy degrades).
    //
    // The fixed pressure is the mean fluid pressure. OpenFOAM normalises
    // residuals with sum |A(p - <p>)|: solid rows with a large penalty and
    // a pressure far from <p> would inflate the normalisation and make the
    // relative tolerance unreachable (round-off).
    scalarField& D = pEqn.diag();
    scalarField& b = pEqn.source();
    scalarField& psi =
        const_cast<volScalarField&>(pEqn.psi()).primitiveFieldRef();

    scalar sumD = 0;
    scalar sumP = 0;
    label nD = 0;

    forAll(D, celli)
    {
        if (status_[celli] != SOLID)
        {
            sumD += mag(D[celli]);
            sumP += psi[celli];
            ++nD;
        }
    }

    reduce(sumD, sumOp<scalar>());
    reduce(sumP, sumOp<scalar>());
    reduce(nD, sumOp<label>());

    // Fallback for solid cells without solid neighbours
    const scalar dFallback = sumD/max(nD, 1) + VSMALL;
    const scalar pMean = sumP/max(nD, 1);

    // The Laplacian diagonal is negative: add -penalty (p - pMean) = 0
    for (const label celli : solidCells_)
    {
        const scalar dOwn = mag(D[celli]);
        const scalar penalty =
            penaltyFactor_*(dOwn > SMALL*dFallback ? dOwn : dFallback);

        psi[celli] = pMean;
        D[celli] -= penalty;
        b[celli] -= penalty*pMean;
    }
}


void Foam::immersedBoundary::constrainU
(
    fvVectorMatrix& UEqn,
    const bool addWallSink
) const
{
    constrainUImpl(UEqn, addWallSink, nullptr);
}


void Foam::immersedBoundary::constrainU
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const bool addWallSink
) const
{
    constrainUImpl(UEqn, addWallSink, &rho.primitiveField());
}


void Foam::immersedBoundary::constrainUImpl
(
    fvVectorMatrix& UEqn,
    const bool addWallSink,
    const scalarField* rhoPtr
) const
{
    addProfiling(constrainUProf, "immersedBoundary::constrainU");

    const pointField& C = mesh().C().primitiveField();

    if (wallModel_->stressMode())
    {
        // IB cells solve their momentum; the wall acts through the wall
        // shear stress on the blocked faces (implicit linearisation
        // tau_w = c (U - U_w), c = u_tau^2/|u_t| times the wall area)
        decoupleSolid(UEqn);

        if (addWallSink)
        {
            scalarField& D = UEqn.diag();
            vectorField& b = UEqn.source();

            const scalarField& coeff = wallModel_->wallCoeff();
            const vectorField& Uw = wallModel_->wallVelocity();
            const vectorField& U = UEqn.psi().primitiveField();
            const vectorField& nw = normal_.primitiveField();

            forAll(ibCells_, i)
            {
                const label celli = ibCells_[i];
                scalar c = coeff[i]*ibWallArea_[i];
                if (rhoPtr)
                {
                    c *= (*rhoPtr)[celli];
                }

                D[celli] += c;
                b[celli] += c*Uw[i];

                if (tangentialSink_)
                {
                    // Return the normal part, c ((U - U_w).n) n, explicitly
                    const vector& n = nw[celli];
                    b[celli] += c*(((U[celli] - Uw[i]) & n)*n);
                }
            }
        }
    }
    else
    {
        UEqn.setValues(ibCells_, wallModel_->Uib());
    }

    UEqn.setValues
    (
        solidCells_,
        solidVelocity(solidCells_, pointField(C, solidCells_))
    );
}


Foam::tmp<Foam::volVectorField>
Foam::immersedBoundary::wallSinkTerm(const volVectorField& U) const
{
    auto tS = volVectorField::New
    (
        "wallSinkTerm",
        IOobject::NO_REGISTER,
        mesh(),
        dimensionedVector(U.dimensions()/dimTime, Zero)
    );

    if (!wallModel_->stressMode())
    {
        return tS;
    }

    vectorField& S = tS.ref().primitiveFieldRef();
    const scalarField& coeff = wallModel_->wallCoeff();
    const vectorField& Uw = wallModel_->wallVelocity();
    const vectorField& nw = normal_.primitiveField();
    const scalarField& V = mesh().V();

    forAll(ibCells_, i)
    {
        const label celli = ibCells_[i];
        const scalar c = coeff[i]*ibWallArea_[i];
        vector dU(U[celli] - Uw[i]);

        if (tangentialSink_)
        {
            dU -= (dU & nw[celli])*nw[celli];
        }

        S[celli] = c*dU/V[celli];
    }

    return tS;
}


void Foam::immersedBoundary::correctU(volVectorField& U) const
{
    addProfiling(correctUProf, "immersedBoundary::correctU");

    const pointField& C = mesh().C().primitiveField();
    vectorField& Ui = U.primitiveFieldRef();

    if (!wallModel_->stressMode())
    {
        const vectorField& Uib = wallModel_->Uib();

        forAll(ibCells_, i)
        {
            Ui[ibCells_[i]] = Uib[i];
        }
    }

    const vectorField Us
    (
        solidVelocity(solidCells_, pointField(C, solidCells_))
    );

    forAll(solidCells_, i)
    {
        Ui[solidCells_[i]] = Us[i];
    }

    U.correctBoundaryConditions();
}


Foam::tmp<Foam::volVectorField>
Foam::immersedBoundary::grad(const volScalarField& p) const
{
    return (gradMode_ > 0 ? faceGrad(p) : gaussGrad(p));
}


Foam::tmp<Foam::volVectorField>
Foam::immersedBoundary::faceGrad(const volScalarField& p) const
{
    addProfiling(faceGradProf, "immersedBoundary::faceGrad");

    const fvMesh& m = mesh();

    // Solid cells (and live cells for faceIB) keep the Gauss gradient
    tmp<volVectorField> tgradP(gaussGrad(p));
    volVectorField& gradP = tgradP.ref();

    // Open-face gradient fluxes, as in the flux correction of the pressure
    // equation: alpha_f r_f snGrad(p) |S_f|
    surfaceScalarField F(faceMask_*fvc::snGrad(p)*m.magSf());
    if (centroidPressure_ && cutFaces_)
    {
        F *= cutDeltaRatio_;
    }

    openFaceReconstruct(F, gradMode_, gradP.primitiveFieldRef());

    gradP.correctBoundaryConditions();

    return tgradP;
}


void Foam::immersedBoundary::openFaceReconstruct
(
    const surfaceScalarField& F,
    const label mode,
    vectorField& g
) const
{
    const fvMesh& m = mesh();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const vectorField& Sf = m.faceAreas();
    const scalarField& magSf = m.magSf();

    // Least-squares face reconstruction weighted by the open fraction:
    // T = sum alpha_f |S_f| n_f n_f, b = sum n_f F_f
    symmTensorField T(m.nCells(), Zero);
    vectorField b(m.nCells(), Zero);

    const scalarField& mi = faceMask_.primitiveField();
    const scalarField& Fi = F.primitiveField();

    // Velocity reconstruction (mode 2): the blocked faces enter the face
    // tensor as zero-gradient data (a wall has no normal pressure
    // gradient), so that the least-squares fit of a cell with an
    // asymmetric set of open faces is not biased. Fitted over the open
    // faces alone, the pressure flux that cancels the cell's own HbyA flux
    // reconstructed into an acceleration along the cell velocity, 2 % per
    // step in a pit cell of the wave basin (momentum budget, 2026-11-01)
    const bool allFaces = (mode == 2);
    forAll(nei, facei)
    {
        const vector nf(Sf[facei]/magSf[facei]);
        const symmTensor t((allFaces ? scalar(1) : mi[facei])*magSf[facei]*sqr(nf));

        T[own[facei]] += t;
        T[nei[facei]] += t;
        b[own[facei]] += nf*Fi[facei];
        b[nei[facei]] += nf*Fi[facei];
    }

    forAll(m.boundaryMesh(), patchi)
    {
        const polyPatch& pp = m.boundaryMesh()[patchi];
        const labelUList& fc = pp.faceCells();
        const fvsPatchScalarField& pm = faceMask_.boundaryField()[patchi];
        const fvsPatchScalarField& pF = F.boundaryField()[patchi];

        // Patches without face values (empty) close the cell with zero flux
        const bool valued = (pF.size() == pp.size());

        forAll(fc, i)
        {
            const label facei = pp.start() + i;
            const scalar a = mag(Sf[facei]);
            const vector nf(Sf[facei]/a);

            T[fc[i]] += (valued && !allFaces ? pm[i] : scalar(1))*a*sqr(nf);
            if (valued)
            {
                b[fc[i]] += nf*pF[i];
            }
        }
    }


    label nSingular = 0;

    forAll(status_, celli)
    {
        if
        (
            status_[celli] == SOLID
         || (mode == 1 && status_[celli] != IB)
        )
        {
            continue;
        }

        const symmTensor& t = T[celli];
        const scalar s = tr(t)/3;

        if (det(t) > 1e-9*pow3(s))
        {
            g[celli] = inv(t) & b[celli];
        }
        else
        {
            ++nSingular;
        }
    }

    if (returnReduce(nSingular, sumOp<label>()))
    {
        WarningInFunction
            << returnReduce(nSingular, sumOp<label>())
            << " cells with a singular face tensor keep the Gauss gradient"
            << endl;
    }

}


Foam::tmp<Foam::volVectorField>
Foam::immersedBoundary::reconstruct(const surfaceScalarField& ssf) const
{
    addProfiling(reconstructProf, "immersedBoundary::reconstruct");

    auto tv = volVectorField::New
    (
        "reconstruct(" + ssf.name() + ')',
        IOobject::NO_REGISTER,
        mesh(),
        dimensionedVector(ssf.dimensions()/dimArea, Zero),
        fvPatchFieldBase::extrapolatedCalculatedType()
    );

    const surfaceScalarField F(faceMask_*ssf);
    openFaceReconstruct(F, 2, tv.ref().primitiveFieldRef());
    tv.ref().correctBoundaryConditions();

    return tv;
}


void Foam::immersedBoundary::correctSolidAlpha(volScalarField& alpha) const
{
    if (returnReduce(solidCells_.empty(), andOp<bool>()))
    {
        return;
    }

    const fvMesh& m = mesh();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    scalarField& a = alpha.primitiveFieldRef();

    // Cells with a value: the fluid cells, then the solid layers in turn
    volScalarField assigned
    (
        IOobject
        (
            "ibAlphaAssigned",
            m.time().timeName(),
            m,
            IOobject::NO_REGISTER
        ),
        m,
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    );
    scalarField& as = assigned.primitiveFieldRef();
    forAll(status_, celli)
    {
        if (status_[celli] != SOLID)
        {
            as[celli] = 1;
        }
    }

    const label nSweeps = 3;
    for (label sweep = 0; sweep < nSweeps; ++sweep)
    {
        assigned.correctBoundaryConditions();
        alpha.correctBoundaryConditions();

        scalarField sum(m.nCells(), Zero);
        scalarField cnt(m.nCells(), Zero);
        forAll(nei, facei)
        {
            const label o = own[facei];
            const label n = nei[facei];
            if (status_[o] == SOLID && as[n] > 0.5)
            {
                sum[o] += a[n];
                cnt[o] += 1;
            }
            if (status_[n] == SOLID && as[o] > 0.5)
            {
                sum[n] += a[o];
                cnt[n] += 1;
            }
        }
        forAll(alpha.boundaryField(), patchi)
        {
            const fvPatchScalarField& pa = alpha.boundaryField()[patchi];
            if (!pa.coupled())
            {
                continue;
            }
            const scalarField an(pa.patchNeighbourField());
            const scalarField asn
            (
                assigned.boundaryField()[patchi].patchNeighbourField()
            );
            const labelUList& fc = pa.patch().faceCells();
            forAll(fc, i)
            {
                if (status_[fc[i]] == SOLID && asn[i] > 0.5)
                {
                    sum[fc[i]] += an[i];
                    cnt[fc[i]] += 1;
                }
            }
        }

        label nNew = 0;
        forAll(solidCells_, i)
        {
            const label celli = solidCells_[i];
            if (as[celli] < 0.5 && cnt[celli] > 0)
            {
                a[celli] = sum[celli]/cnt[celli];
                as[celli] = 1;
                ++nNew;
            }
        }
        if (!returnReduce(nNew, sumOp<label>()))
        {
            break;
        }
    }

    alpha.correctBoundaryConditions();
}


Foam::tmp<Foam::volVectorField>
Foam::immersedBoundary::gaussGrad(const volScalarField& p) const
{
    addProfiling(gradProf, "immersedBoundary::grad");

    const fvMesh& m = mesh();

    tmp<volVectorField> tgradP(fvc::grad(p));
    volVectorField& gradP = tgradP.ref();

    if (ibCells_.empty() && !Pstream::parRun())
    {
        return tgradP;
    }

    // Gauss gradient of the IB cells with the blocked-face values replaced
    // by the cell value (zero normal gradient towards the solid)
    const surfaceScalarField pf
    (
        faceInterpolation_ ? cutInterpolate(p) : fvc::interpolate(p)
    );
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const vectorField& Sf = m.faceAreas();
    const scalarField& V = m.V();

    labelList ibIndex(m.nCells(), -1);
    forAll(ibCells_, i)
    {
        ibIndex[ibCells_[i]] = i;
    }

    vectorField sumPSf(ibCells_.size(), Zero);

    forAll(nei, facei)
    {
        const label o = own[facei];
        const label nb = nei[facei];
        const bool open = (faceMask_[facei] > 0.5);

        if (cutFaces_)
        {
            // Open part at the face value, closed part (wall) at the cell
            // value (zero normal gradient)
            const scalar a = faceMask_[facei];

            if (ibIndex[o] >= 0)
            {
                sumPSf[ibIndex[o]] += (a*pf[facei] + (1 - a)*p[o])*Sf[facei];
            }
            if (ibIndex[nb] >= 0)
            {
                sumPSf[ibIndex[nb]] -=
                    (a*pf[facei] + (1 - a)*p[nb])*Sf[facei];
            }
        }
        else
        {
            if (ibIndex[o] >= 0)
            {
                sumPSf[ibIndex[o]] += (open ? pf[facei] : p[o])*Sf[facei];
            }
            if (ibIndex[nb] >= 0)
            {
                sumPSf[ibIndex[nb]] -= (open ? pf[facei] : p[nb])*Sf[facei];
            }
        }
    }

    forAll(m.boundary(), patchi)
    {
        const fvPatch& fp = m.boundary()[patchi];
        const labelUList& fc = fp.faceCells();
        const fvsPatchScalarField& ppf = pf.boundaryField()[patchi];
        const fvsPatchScalarField& pm = faceMask_.boundaryField()[patchi];
        const vectorField& pSf = fp.Sf();

        forAll(fc, i)
        {
            const label ibi = ibIndex[fc[i]];

            if (ibi >= 0)
            {
                sumPSf[ibi] +=
                (
                    cutFaces_
                  ? (pm[i]*ppf[i] + (1 - pm[i])*p[fc[i]])
                  : (pm[i] > 0.5 ? ppf[i] : p[fc[i]])
                )*pSf[i];
            }
        }
    }

    forAll(ibCells_, i)
    {
        gradP[ibCells_[i]] = sumPSf[i]/V[ibCells_[i]];
    }

    gradP.correctBoundaryConditions();

    return tgradP;
}


void Foam::immersedBoundary::checkFlux
(
    const surfaceScalarField& phi,
    const volVectorField& U
) const
{
    addProfiling(checkFluxProf, "immersedBoundary::checkFlux");

    // Flux on blocked internal faces (zero by construction)
    scalar blocked = 0;

    forAll(faceMask_, facei)
    {
        if (faceMask_[facei] < 0.5)
        {
            blocked += mag(phi[facei]);
        }
    }

    // Boundary velocity flux into solid cells (e.g. an inlet profile that
    // reaches below the bed; it is masked and therefore lost) and the total
    // boundary inflow of the fluid
    scalar solidBoundary = 0;
    scalar inflow = 0;

    forAll(phi.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pphi = phi.boundaryField()[patchi];
        const fvsPatchScalarField& pm = faceMask_.boundaryField()[patchi];
        const fvPatchVectorField& pU = U.boundaryField()[patchi];
        const vectorField& Sf = pphi.patch().Sf();

        if (pphi.coupled())
        {
            continue;
        }

        forAll(pphi, i)
        {
            if (pm[i] < 0.5)
            {
                solidBoundary += mag(pU[i] & Sf[i]);
            }
            else if (pphi[i] < 0)
            {
                inflow -= pphi[i];
            }
        }
    }

    reduce(blocked, sumOp<scalar>());
    reduce(solidBoundary, sumOp<scalar>());
    reduce(inflow, sumOp<scalar>());

    Info<< "IB flux check: blocked faces sum|phi| = " << blocked
        << ", boundary flux into solid cells sum|U.Sf| = " << solidBoundary
        << " (fluid inflow " << inflow << ")" << endl;
}


void Foam::immersedBoundary::makeCutFaceData()
{
    const fvMesh& m = mesh();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const pointField& C = m.cellCentres();
    const vectorField& Cf = m.faceCentres();
    const vectorField& Sf = m.faceAreas();
    const scalarField& theta = theta_.primitiveField();
    const surfaceScalarField& mw = m.weights();

    scalarField& flag = cutFaceFlag_.primitiveFieldRef();
    scalarField& w = cutWeights_.primitiveFieldRef();
    vectorField& io = cutInterpOffset_.primitiveFieldRef();
    vectorField& kv = cutSkewVec_.primitiveFieldRef();

    flag = 0;
    w = mw.primitiveField();
    io = Zero;
    kv = Zero;

    label nDegenerate = 0;

    // Weight, offset and skew vector of one face from the owner (P) and
    // neighbour (N) centres and fluid-centroid offsets; false for a
    // degenerate face
    auto faceData = [&]
    (
        const label facei,
        const vector& CP,
        const vector& offP,
        const vector& CN,
        const vector& offN,
        scalar& wf,
        vector& iof,
        vector& kf
    ) -> bool
    {
        // Sf covers all faces (boundary faces included), magSf() would not
        const vector nf(Sf[facei]/max(mag(Sf[facei]), VSMALL));
        const vector xP(CP + offP);
        const vector xN(CN + offN);
        const vector Delta(xN - xP);
        const scalar d = ((CN - CP) & nf);
        const scalar dn = (Delta & nf);
        const scalar dc = max(dn, minDistanceFraction_*d);
        const vector xfo(Cf[facei] + faceOpenOffset_[facei]);

        if (d <= VSMALL || dc <= VSMALL)
        {
            // Degenerate geometry (coincident centres or centroids): leave
            // the face uncorrected
            ++nDegenerate;
            iof = Zero;
            kf = Zero;
            return false;
        }

        if (dn > VSMALL)
        {
            wf = min(max(((xN - xfo) & nf)/dn, scalar(0)), scalar(1));
        }
        iof = xfo - (xP + (1 - wf)*Delta);
        kf = nf - Delta/dc;
        return true;
    };

    forAll(nei, facei)
    {
        const label o = own[facei];
        const label nb = nei[facei];

        if
        (
            status_[o] == SOLID || status_[nb] == SOLID
         || (theta[o] >= 1 && theta[nb] >= 1)
        )
        {
            continue;
        }

        flag[facei] =
        (
            faceData
            (
                facei,
                C[o], fluidCentreOffset_[o], C[nb], fluidCentreOffset_[nb],
                w[facei], io[facei], kv[facei]
            ) ? 1 : 0
        );
    }

    labelList nbrStatus;
    syncTools::swapBoundaryCellList(m, status_, nbrStatus);

    pointField nbrC;
    syncTools::swapBoundaryCellPositions(m, C, nbrC);

    vectorField nbrOffset;
    syncTools::swapBoundaryCellList(m, fluidCentreOffset_, nbrOffset);

    scalarField nbrTheta;
    syncTools::swapBoundaryCellList(m, theta, nbrTheta);

    forAll(m.boundary(), patchi)
    {
        const fvPatch& p = m.boundary()[patchi];
        fvsPatchScalarField& pflag = cutFaceFlag_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pw = cutWeights_.boundaryFieldRef()[patchi];
        fvsPatchVectorField& pio = cutInterpOffset_.boundaryFieldRef()[patchi];
        fvsPatchVectorField& pkv = cutSkewVec_.boundaryFieldRef()[patchi];

        pflag = 0;
        pw = mw.boundaryField()[patchi];
        pio = Zero;
        pkv = Zero;

        if (!p.coupled())
        {
            continue;
        }

        const labelUList& fc = p.faceCells();

        forAll(fc, i)
        {
            const label facei = p.start() + i;
            const label bFacei = facei - m.nInternalFaces();
            const label o = fc[i];

            if
            (
                status_[o] == SOLID || nbrStatus[bFacei] == SOLID
             || (theta[o] >= 1 && nbrTheta[bFacei] >= 1)
            )
            {
                continue;
            }

            pflag[i] =
            (
                faceData
                (
                    facei,
                    C[o], fluidCentreOffset_[o], nbrC[bFacei], nbrOffset[bFacei],
                    pw[i], pio[i], pkv[i]
                ) ? 1 : 0
            );
        }
    }

    reduce(nDegenerate, sumOp<label>());

    if (nDegenerate)
    {
        Info<< "    cut faces: " << nDegenerate
            << " degenerate faces left uncorrected" << endl;
    }
}


Foam::tmp<Foam::volVectorField> Foam::immersedBoundary::cutSkewCorrection
(
    const volVectorField& U,
    const volScalarField& gamma
) const
{
    const fvMesh& m = mesh();

    auto tS = volVectorField::New
    (
        "cutSkewCorrection",
        IOobject::NO_REGISTER,
        m,
        dimensionedVector(U.dimensions()/dimTime, Zero)
    );

    if (!(cutFaces_ && skewCorrection_))
    {
        return tS;
    }

    addProfiling(cutSkewProf, "immersedBoundary::cutSkewCorrection");

    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const scalarField& magSf = m.magSf();
    const scalarField& V = m.V();

    const tmp<volTensorField> tgradU(cutLeastSquaresGrad(U));
    const tensorField& gU = tgradU().primitiveField();
    const surfaceScalarField gammaf(fvc::interpolate(gamma));
    const scalarField& flag = cutFaceFlag_.primitiveField();
    const vectorField& kv = cutSkewVec_.primitiveField();
    const scalarField& mask = faceMask_.primitiveField();

    vectorField& S = tS.ref().primitiveFieldRef();

    forAll(nei, facei)
    {
        if (flag[facei] < 0.5)
        {
            continue;
        }

        const label o = own[facei];
        const label nb = nei[facei];

        // gamma_f alpha_f |S_f| k_f . grad(U)_f, into the owner, out of the
        // neighbour
        const vector c
        (
            gammaf[facei]*mask[facei]*magSf[facei]
           *(kv[facei] & (0.5*(gU[o] + gU[nb])))
        );

        S[o] += c;
        S[nb] -= c;
    }

    tensorField nbrG;
    syncTools::swapBoundaryCellList(m, gU, nbrG);

    forAll(m.boundary(), patchi)
    {
        const fvPatch& p = m.boundary()[patchi];

        if (!p.coupled())
        {
            continue;
        }

        const labelUList& fc = p.faceCells();
        const fvsPatchScalarField& pflag = cutFaceFlag_.boundaryField()[patchi];
        const fvsPatchVectorField& pkv = cutSkewVec_.boundaryField()[patchi];
        const scalarField& pg = gammaf.boundaryField()[patchi];
        const fvsPatchScalarField& pm = faceMask_.boundaryField()[patchi];
        const scalarField& pMagSf = p.magSf();

        forAll(fc, i)
        {
            if (pflag[i] < 0.5)
            {
                continue;
            }

            const label bFacei = p.start() + i - m.nInternalFaces();

            S[fc[i]] +=
                pg[i]*pm[i]*pMagSf[i]
               *(pkv[i] & (0.5*(gU[fc[i]] + nbrG[bFacei])));
        }
    }

    S /= V;

    return tS;
}


Foam::tmp<Foam::surfaceScalarField> Foam::immersedBoundary::cutInterpolate
(
    const volScalarField& vf
) const
{
    tmp<surfaceScalarField> tf(fvc::interpolate(vf));

    if (!(cutFaces_ && faceInterpolation_))
    {
        return tf;
    }

    const fvMesh& m = mesh();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const scalarField& flag = cutFaceFlag_.primitiveField();
    const scalarField& w = cutWeights_.primitiveField();
    const scalarField& v = vf.primitiveField();

    scalarField& f = tf.ref().primitiveFieldRef();

    forAll(nei, facei)
    {
        if (flag[facei] > 0.5)
        {
            f[facei] = w[facei]*v[own[facei]] + (1 - w[facei])*v[nei[facei]];
        }
    }

    scalarField nbrV;
    syncTools::swapBoundaryCellList(m, v, nbrV);

    forAll(m.boundary(), patchi)
    {
        const fvPatch& p = m.boundary()[patchi];

        if (!p.coupled())
        {
            continue;
        }

        const labelUList& fc = p.faceCells();
        const fvsPatchScalarField& pflag = cutFaceFlag_.boundaryField()[patchi];
        const fvsPatchScalarField& pw = cutWeights_.boundaryField()[patchi];
        fvsPatchScalarField& pf = tf.ref().boundaryFieldRef()[patchi];

        forAll(fc, i)
        {
            if (pflag[i] > 0.5)
            {
                const label bFacei = p.start() + i - m.nInternalFaces();
                pf[i] = pw[i]*v[fc[i]] + (1 - pw[i])*nbrV[bFacei];
            }
        }
    }

    return tf;
}


void Foam::immersedBoundary::correctCutFlux
(
    const volVectorField& vf,
    surfaceScalarField& phi
) const
{
    if (!(cutFaces_ && faceInterpolation_))
    {
        return;
    }

    addProfiling(cutFluxProf, "immersedBoundary::correctCutFlux");

    const fvMesh& m = mesh();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const vectorField& Sf = m.faceAreas();
    const scalarField& flag = cutFaceFlag_.primitiveField();
    const scalarField& w = cutWeights_.primitiveField();
    const vectorField& io = cutInterpOffset_.primitiveField();
    const vectorField& v = vf.primitiveField();

    const tmp<volTensorField> tg(cutLeastSquaresGrad(vf));
    const tensorField& g = tg().primitiveField();

    scalarField& phii = phi.primitiveFieldRef();

    forAll(nei, facei)
    {
        if (flag[facei] < 0.5)
        {
            continue;
        }

        const label o = own[facei];
        const label nb = nei[facei];

        // value at the open-face centroid: along the centroid line, plus
        // the tangential offset times the face gradient
        const vector vff
        (
            w[facei]*v[o] + (1 - w[facei])*v[nb]
          + (io[facei] & (0.5*(g[o] + g[nb])))
        );

        phii[facei] = (vff & Sf[facei]);
    }

    vectorField nbrV;
    syncTools::swapBoundaryCellList(m, v, nbrV);

    tensorField nbrG;
    syncTools::swapBoundaryCellList(m, g, nbrG);

    forAll(m.boundary(), patchi)
    {
        const fvPatch& p = m.boundary()[patchi];

        if (!p.coupled())
        {
            continue;
        }

        const labelUList& fc = p.faceCells();
        const fvsPatchScalarField& pflag = cutFaceFlag_.boundaryField()[patchi];
        const fvsPatchScalarField& pw = cutWeights_.boundaryField()[patchi];
        const fvsPatchVectorField& pio = cutInterpOffset_.boundaryField()[patchi];
        const vectorField& pSf = p.Sf();
        fvsPatchScalarField& pphi = phi.boundaryFieldRef()[patchi];

        forAll(fc, i)
        {
            if (pflag[i] < 0.5)
            {
                continue;
            }

            const label bFacei = p.start() + i - m.nInternalFaces();

            const vector vff
            (
                pw[i]*v[fc[i]] + (1 - pw[i])*nbrV[bFacei]
              + (pio[i] & (0.5*(g[fc[i]] + nbrG[bFacei])))
            );

            pphi[i] = (vff & pSf[i]);
        }
    }
}


Foam::tmp<Foam::volVectorField> Foam::immersedBoundary::cutConvectionCorrection
(
    const surfaceScalarField& phi,
    const volVectorField& U
) const
{
    const fvMesh& m = mesh();

    auto tS = volVectorField::New
    (
        "cutConvectionCorrection",
        IOobject::NO_REGISTER,
        m,
        dimensionedVector(U.dimensions()/dimTime, Zero)
    );

    if (!(cutFaces_ && faceInterpolation_))
    {
        return tS;
    }

    addProfiling(cutConvProf, "immersedBoundary::cutConvectionCorrection");

    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const scalarField& mw = m.weights().primitiveField();
    const scalarField& V = m.V();
    const scalarField& flag = cutFaceFlag_.primitiveField();
    const scalarField& w = cutWeights_.primitiveField();
    const vectorField& io = cutInterpOffset_.primitiveField();
    const vectorField& v = U.primitiveField();
    const scalarField& phii = phi.primitiveField();

    const tmp<volTensorField> tg(cutLeastSquaresGrad(U));
    const tensorField& g = tg().primitiveField();

    vectorField& S = tS.ref().primitiveFieldRef();

    forAll(nei, facei)
    {
        if (flag[facei] < 0.5)
        {
            continue;
        }

        const label o = own[facei];
        const label nb = nei[facei];

        const vector dU
        (
            (w[facei] - mw[facei])*(v[o] - v[nb])
          + (io[facei] & (0.5*(g[o] + g[nb])))
        );

        // div(phi U): + phi_f U_f for the owner, - for the neighbour
        const vector c(phii[facei]*dU);
        S[o] -= c;
        S[nb] += c;
    }

    vectorField nbrV;
    syncTools::swapBoundaryCellList(m, v, nbrV);

    tensorField nbrG;
    syncTools::swapBoundaryCellList(m, g, nbrG);

    forAll(m.boundary(), patchi)
    {
        const fvPatch& p = m.boundary()[patchi];

        if (!p.coupled())
        {
            continue;
        }

        const labelUList& fc = p.faceCells();
        const fvsPatchScalarField& pflag = cutFaceFlag_.boundaryField()[patchi];
        const fvsPatchScalarField& pw = cutWeights_.boundaryField()[patchi];
        const scalarField& pmw = m.weights().boundaryField()[patchi];
        const fvsPatchVectorField& pio = cutInterpOffset_.boundaryField()[patchi];
        const fvsPatchScalarField& pphi = phi.boundaryField()[patchi];

        forAll(fc, i)
        {
            if (pflag[i] < 0.5)
            {
                continue;
            }

            const label bFacei = p.start() + i - m.nInternalFaces();

            const vector dU
            (
                (pw[i] - pmw[i])*(v[fc[i]] - nbrV[bFacei])
              + (pio[i] & (0.5*(g[fc[i]] + nbrG[bFacei])))
            );

            S[fc[i]] -= pphi[i]*dU;
        }
    }

    S /= V;

    return tS;
}


// ************************************************************************* //
