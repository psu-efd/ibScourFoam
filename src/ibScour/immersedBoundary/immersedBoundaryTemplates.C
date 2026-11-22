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

#include "fvMatrix.H"
#include "fvcDiv.H"
#include "syncTools.H"
#include "surfaceInterpolate.H"
#include "bitSet.H"
#include "symmTensorField.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Eigen-decomposition of a symmetric tensor by cyclic Jacobi rotations:
// eigenvalues ascending, eigenvectors as the rows of evec.
// Foam::eigenValues solves the characteristic cubic analytically and reports
// non-real roots for a degenerate spectrum with round-off off-diagonals,
// e.g. the stencil matrix diag(2, 2, 1) of a cut cell with one missing
// neighbour; the truncated inverse then drops a resolved direction.
void symmEigen
(
    const Foam::symmTensor& S,
    Foam::vector& ev,
    Foam::tensor& evec
)
{
    Foam::scalar a[3][3] =
    {
        {S.xx(), S.xy(), S.xz()},
        {S.xy(), S.yy(), S.yz()},
        {S.xz(), S.yz(), S.zz()}
    };
    Foam::scalar v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    for (int sweep = 0; sweep < 50; ++sweep)
    {
        const Foam::scalar off =
            Foam::sqr(a[0][1]) + Foam::sqr(a[0][2]) + Foam::sqr(a[1][2]);
        const Foam::scalar dia =
            Foam::sqr(a[0][0]) + Foam::sqr(a[1][1]) + Foam::sqr(a[2][2]);

        if (off <= 1e-30*dia)
        {
            break;
        }

        for (int p = 0; p < 2; ++p)
        {
            for (int q = p + 1; q < 3; ++q)
            {
                if (Foam::mag(a[p][q]) < Foam::VSMALL)
                {
                    continue;
                }

                const Foam::scalar th = (a[q][q] - a[p][p])/(2*a[p][q]);
                const Foam::scalar t =
                    Foam::sign(th)
                   /(Foam::mag(th) + Foam::sqrt(Foam::sqr(th) + 1));
                const Foam::scalar c = 1/Foam::sqrt(Foam::sqr(t) + 1);
                const Foam::scalar sn = t*c;

                for (int k = 0; k < 3; ++k)
                {
                    const Foam::scalar akp = a[k][p];
                    const Foam::scalar akq = a[k][q];
                    a[k][p] = c*akp - sn*akq;
                    a[k][q] = sn*akp + c*akq;
                }
                for (int k = 0; k < 3; ++k)
                {
                    const Foam::scalar apk = a[p][k];
                    const Foam::scalar aqk = a[q][k];
                    a[p][k] = c*apk - sn*aqk;
                    a[q][k] = sn*apk + c*aqk;
                }
                for (int k = 0; k < 3; ++k)
                {
                    const Foam::scalar vkp = v[k][p];
                    const Foam::scalar vkq = v[k][q];
                    v[k][p] = c*vkp - sn*vkq;
                    v[k][q] = sn*vkp + c*vkq;
                }
            }
        }
    }

    // Ascending order; the eigenvectors are the columns of v
    Foam::label order[3] = {0, 1, 2};
    for (int i = 0; i < 2; ++i)
    {
        for (int j = i + 1; j < 3; ++j)
        {
            if (a[order[j]][order[j]] < a[order[i]][order[i]])
            {
                std::swap(order[i], order[j]);
            }
        }
    }

    for (int i = 0; i < 3; ++i)
    {
        const Foam::label c = order[i];
        ev[i] = a[c][c];
        evec(i, 0) = v[0][c];
        evec(i, 1) = v[1][c];
        evec(i, 2) = v[2][c];
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::immersedBoundary::maskNonOrthCorrection
(
    fvMatrix<Type>& eqn,
    const surfaceScalarField& mask
) const
{
    if (!maskNonOrthCorr_ || !eqn.faceFluxCorrectionPtr())
    {
        return;
    }
    auto& corr = *eqn.faceFluxCorrectionPtr();
    const GeometricField<Type, fvsPatchField, surfaceMesh> masked(mask*corr);
    // The Laplacian added -V div(corr) to the source
    eqn.source() += eqn.psi().mesh().V()*fvc::div(masked)().primitiveField();
    corr -= masked;
}


template<class Type>
void Foam::immersedBoundary::decoupleSolid(fvMatrix<Type>& eqn) const
{
    const fvMesh& m = mesh();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();

    scalarField& D = eqn.diag();
    scalarField& upper = eqn.upper();

    // Lower coefficients exist only for asymmetric matrices
    const bool asymmetric = eqn.hasLower();
    scalarField& lower = (asymmetric ? eqn.lower() : upper);

    forAll(nei, facei)
    {
        const bool solidO = (status_[own[facei]] == SOLID);
        const bool solidN = (status_[nei[facei]] == SOLID);

        if (solidO == solidN)
        {
            continue;
        }

        // A face contributes -offDiag to the diagonal of each row (no flux
        // through a blocked face, so only the diffusive part): remove it
        // from the fluid row; the solid row is fixed separately
        if (!solidO)
        {
            D[own[facei]] += upper[facei];
        }
        else
        {
            D[nei[facei]] += lower[facei];
        }

        upper[facei] = 0;
        if (asymmetric)
        {
            lower[facei] = 0;
        }
    }

    // Processor faces between a fluid and a solid cell
    labelList nbrStatus;
    syncTools::swapBoundaryCellList(m, status_, nbrStatus);

    forAll(m.boundary(), patchi)
    {
        const fvPatch& p = m.boundary()[patchi];

        if (!p.coupled())
        {
            continue;
        }

        Field<Type>& intCoeffs = eqn.internalCoeffs()[patchi];
        Field<Type>& bouCoeffs = eqn.boundaryCoeffs()[patchi];
        const labelUList& fc = p.faceCells();

        forAll(fc, i)
        {
            const label bFacei = p.start() + i - m.nInternalFaces();

            if (status_[fc[i]] != SOLID && nbrStatus[bFacei] == SOLID)
            {
                intCoeffs[i] = Zero;
                bouCoeffs[i] = Zero;
            }
        }
    }
}


template<class Type>
void Foam::immersedBoundary::scaleCutDiffusion
(
    fvMatrix<Type>& eqn,
    const volScalarField& gamma
) const
{
    if (!cutFaces_)
    {
        return;
    }

    const fvMesh& m = mesh();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();

    const surfaceScalarField gammaf(fvc::interpolate(gamma));
    const surfaceScalarField& dc = m.nonOrthDeltaCoeffs();
    const scalarField& magSf = m.magSf();

    scalarField& D = eqn.diag();
    scalarField& upper = eqn.upper();

    const bool asymmetric = eqn.hasLower();
    scalarField& lower = (asymmetric ? eqn.lower() : upper);

    // Distance ratio d/d_c of the fluid centroids (1 without cutCentroids)
    const bool centroids = centroidDiffusion_;
    const scalarField& r = cutDeltaRatio_.primitiveField();

    // -laplacian(gamma, psi) contributes -c to the off-diagonals and +c to
    // both diagonals of a face, c = gamma_f |S_f| delta_f; keep the part
    // alpha_f r_f c (open part, centroid distance), i.e. remove
    // (1 - alpha_f r_f) c. Faces to solid cells are handled by
    // decoupleSolid.
    forAll(nei, facei)
    {
        const scalar ar = faceMask_[facei]*(centroids ? r[facei] : 1);

        if
        (
            ar == 1
         || (!centroids && ar > 1)
         || status_[own[facei]] == SOLID
         || status_[nei[facei]] == SOLID
        )
        {
            continue;
        }

        const scalar c = (1 - ar)*gammaf[facei]*magSf[facei]*dc[facei];

        upper[facei] += c;
        if (asymmetric)
        {
            lower[facei] += c;
        }
        D[own[facei]] -= c;
        D[nei[facei]] -= c;
    }

    // Coupled faces: -laplacian contributes +c (times one) to both the
    // internal and the boundary coefficients
    forAll(m.boundary(), patchi)
    {
        const fvPatch& p = m.boundary()[patchi];

        if (!p.coupled())
        {
            continue;
        }

        const fvsPatchScalarField& pa = faceMask_.boundaryField()[patchi];
        const fvsPatchScalarField& pr = cutDeltaRatio_.boundaryField()[patchi];
        const scalarField& pg = gammaf.boundaryField()[patchi];
        const scalarField& pdc = dc.boundaryField()[patchi];
        const scalarField& pMagSf = p.magSf();
        const labelUList& fc = p.faceCells();

        Field<Type>& intCoeffs = eqn.internalCoeffs()[patchi];
        Field<Type>& bouCoeffs = eqn.boundaryCoeffs()[patchi];

        forAll(fc, i)
        {
            // Faces to solid cells are zeroed afterwards by decoupleSolid
            const scalar ar = pa[i]*(centroids ? pr[i] : 1);

            if
            (
                ar == 1
             || (!centroids && ar > 1)
             || status_[fc[i]] == SOLID
            )
            {
                continue;
            }

            const scalar c = (1 - ar)*pg[i]*pMagSf[i]*pdc[i];

            intCoeffs[i] -= c*pTraits<Type>::one;
            bouCoeffs[i] -= c*pTraits<Type>::one;
        }
    }
}


namespace Foam
{
    //- Least-squares gradient: eigenvalue ratio below which a stencil
    //  direction is treated as unresolved (stencil spread below 10 % of
    //  the resolved directions)
    static constexpr scalar lsEigenRatio = 0.01;
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::immersedBoundary::cutLeastSquaresGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradField;

    const fvMesh& m = mesh();

    auto tg = GradField::New
    (
        "cutGrad(" + vf.name() + ")",
        IOobject::NO_REGISTER,
        m,
        dimensioned<GradType>(vf.dimensions()/dimLength, Zero)
    );
    Field<GradType>& g = tg.ref().primitiveFieldRef();

    if (!cutFaces_ || !(skewCorrection_ || faceInterpolation_))
    {
        return tg;
    }

    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const pointField& C = m.cellCentres();
    const Field<Type>& v = vf.primitiveField();
    const scalarField& flag = cutFaceFlag_.primitiveField();

    // Weighted least squares over the fluid face neighbours, distances
    // between the fluid centroids
    symmTensorField G(m.nCells(), Zero);
    Field<GradType> R(m.nCells(), Zero);
    bitSet active(m.nCells());

    auto accumulate = [&]
    (
        const label celli,
        const vector& d,
        const Type& dv
    )
    {
        const scalar wgt = 1/max(magSqr(d), VSMALL);
        G[celli] += wgt*symm(d*d);
        R[celli] += wgt*(d*dv);
    };

    forAll(nei, facei)
    {
        const label o = own[facei];
        const label nb = nei[facei];

        if (status_[o] == SOLID || status_[nb] == SOLID)
        {
            continue;
        }

        const vector d
        (
            (C[nb] + fluidCentreOffset_[nb]) - (C[o] + fluidCentreOffset_[o])
        );

        accumulate(o, d, v[nb] - v[o]);
        accumulate(nb, -d, v[o] - v[nb]);

        if (flag[facei] > 0.5)
        {
            active.set(o);
            active.set(nb);
        }
    }

    labelList nbrStatus;
    syncTools::swapBoundaryCellList(m, status_, nbrStatus);

    pointField nbrC;
    syncTools::swapBoundaryCellPositions(m, C, nbrC);

    vectorField nbrOffset;
    syncTools::swapBoundaryCellList(m, fluidCentreOffset_, nbrOffset);

    Field<Type> nbrV;
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

        forAll(fc, i)
        {
            const label bFacei = p.start() + i - m.nInternalFaces();
            const label o = fc[i];

            if (status_[o] == SOLID || nbrStatus[bFacei] == SOLID)
            {
                continue;
            }

            const vector d
            (
                (nbrC[bFacei] + nbrOffset[bFacei])
              - (C[o] + fluidCentreOffset_[o])
            );

            accumulate(o, d, nbrV[bFacei] - v[o]);

            if (pflag[i] > 0.5)
            {
                active.set(o);
            }
        }
    }

    nTruncatedGrad_ = 0;

    for (const label celli : active)
    {
        const symmTensor& Gc = G[celli];
        const scalar tr = Gc.xx() + Gc.yy() + Gc.zz();

        if (tr < VSMALL)
        {
            continue;
        }

        // Truncated inverse: a direction in which the stencil has little
        // spread carries no gradient. 2-D: the empty direction. A cut
        // cell between solid cells: the direction of the missing
        // neighbours, where only the centroid offsets of the remaining
        // neighbours contribute. Tikhonov regularisation is not enough
        // there: a spread of 1 % of the cell size turns a small velocity
        // difference into a gradient of order dU/(0.01 h), which the
        // explicit corrections turned into a runaway on the pier mesh
        vector ev;
        tensor evec;
        symmEigen(Gc, ev, evec);                   // ascending, rows
        const scalar evMin = lsEigenRatio*ev.z();

        symmTensor Ginv(Zero);
        bool truncated = false;
        for (direction i = 0; i < vector::nComponents; ++i)
        {
            if (ev[i] > evMin)
            {
                Ginv += sqr(evec.row(i))/ev[i];
            }
            else
            {
                truncated = true;
            }
        }

        if (truncated)
        {
            ++nTruncatedGrad_;
        }

        g[celli] = (Ginv & R[celli]);
    }

    return tg;
}


// ************************************************************************* //
