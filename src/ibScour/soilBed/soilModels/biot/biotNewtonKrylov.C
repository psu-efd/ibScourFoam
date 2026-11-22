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

#include "biot.H"
#include "fvm.H"
#include "fvc.H"
#include "directionMixedFvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Global vector of the Newton-Krylov solve: the D block (3 per cell),
// the p block and the top-face block b (the displacement of the top faces,
// or the defect of their traction condition), with block scales for the
// norms and dot products
struct nkVector
{
    Foam::vectorField D;
    Foam::scalarField p;
    Foam::vectorField b;

    nkVector() = default;

    nkVector(const Foam::label n, const Foam::label nb)
    :
        D(n, Foam::Zero), p(n, Foam::Zero), b(nb, Foam::Zero)
    {}
};


inline Foam::scalar dot
(
    const nkVector& a,
    const nkVector& b,
    const Foam::scalar sD,
    const Foam::scalar sp,
    const Foam::scalar sb
)
{
    return
        Foam::gSum(a.D & b.D)/Foam::sqr(sD)
      + Foam::gSum(a.p*b.p)/Foam::sqr(sp)
      + Foam::gSum(a.b & b.b)/Foam::sqr(sb);
}


inline Foam::scalar norm
(
    const nkVector& a,
    const Foam::scalar sD,
    const Foam::scalar sp,
    const Foam::scalar sb
)
{
    return Foam::sqrt(dot(a, a, sD, sp, sb));
}


inline void axpy(nkVector& y, const Foam::scalar a, const nkVector& x)
{
    y.D += a*x.D;
    y.p += a*x.p;
    y.b += a*x.b;
}

} // End namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::biot::coupledResidual
(
    const volScalarField& impK,
    volScalarField& invM,
    const dimensionedScalar& invMsrc,
    const volScalarField::Internal& fGen,
    const label nSteps,
    const vectorField* bFace,
    vectorField& rD,
    scalarField& rp,
    vectorField& rb
)
{
    // Structure cells keep their prescribed displacement
    if (maskedCells_.size())
    {
        vectorField& D = D_.primitiveFieldRef();
        forAll(maskedCells_, i)
        {
            D[maskedCells_[i]] = maskedDisplacement_[i];
        }
    }

    const fvPatch& topPatch = mesh_.boundary()[topPatchi_];
    auto& Dp =
        refCast<directionMixedFvPatchVectorField>
        (
            D_.boundaryFieldRef()[topPatchi_]
        );

    if (bFace)
    {
        // The face displacement of the traction directions, through the
        // gradient that reproduces it (the prescribed directions ignore
        // the gradient)
        Dp.refGrad() =
            (*bFace - Dp.patchInternalField())*topPatch.deltaCoeffs();
    }

    p_.correctBoundaryConditions();
    D_.correctBoundaryConditions();

    // Stress of the law for this iterate; without face unknowns the
    // traction condition is refreshed from this stress (the fixed-point
    // update of the corrector) and the boundary values follow
    updateStress();
    if (!bFace)
    {
        updateTraction();
        D_.correctBoundaryConditions();
        updateStress();
    }

    fvVectorMatrix DEqn
    (
        fvm::laplacian(impK, D_, "laplacian(DD,D)")
      + divSigmaExp_
     ==
        fvc::grad(p_)
    );

    if (maskedCells_.size())
    {
        DEqn.setValues(maskedCells_, maskedDisplacement_);
    }

    // Residual b - A x in matrix units through the matrix-vector product
    // (M & x = (A x - b)/V), which counts the processor-interface
    // contributions once; fvMatrix::residual() adds them twice (through
    // the interfaces of the ldu residual and again through the coupled
    // boundary source), which breaks the parallel residual
    rD = -(DEqn & D_.internalField())().primitiveField()*mesh_.V();

    // Storage equation with the current iterate in the explicit terms:
    // the stabilisation cancels identically
    invM =
        dimensionedScalar(storage_.dimensions()*dimPressure, 1)
       /law_->stabilisationModulusField();

    fvScalarMatrix pEqn
    (
        fvm::ddt(storage_ + invM, p_)
      - fvm::laplacian(nSteps*mobilityFace(mobility_), p_)
     ==
      - fvc::ddt(divD_)
      + invM*fvc::ddt(p_)
    );

    if (fGen.size())
    {
        pEqn -= nSteps*(storage_ + invMsrc)*fGen;
    }

    rp = -(pEqn & p_.internalField())().primitiveField()*mesh_.V();

    // Traction condition of the top faces, in force units as the D block
    rb = tractionDefect()*topPatch.magSf();
}


Foam::tmp<Foam::vectorField> Foam::biot::tractionDefect() const
{
    const fvPatch& patch = mesh_.boundary()[topPatchi_];
    const auto& Dp =
        refCast<const directionMixedFvPatchVectorField>
        (
            D_.boundaryField()[topPatchi_]
        );
    const vectorField n(patch.nf());
    const symmTensorField sigmaP(sigma_.boundaryField()[topPatchi_]);

    vectorField t(patch.size(), Zero);
    if (topTraction_.size() == patch.size())
    {
        t = topTraction_;
    }
    if (sealedTop_ && pTop_.size() == patch.size())
    {
        t += n*(p_.boundaryField()[topPatchi_].patchInternalField() - pTop_);
    }

    // Traction faces only (the free directions of the value fraction)
    return (symmTensor::I - Dp.valueFraction()) & (t - (n & sigmaP));
}


Foam::scalar Foam::biot::tractionMismatch(const scalar pScale) const
{
    const fvPatch& patch = mesh_.boundary()[topPatchi_];
    const vectorField mismatch(tractionDefect());

    return
        Foam::sqrt(gSum(magSqr(mismatch)*patch.magSf()))
       /max(pScale*Foam::sqrt(gSum(patch.magSf())), SMALL);
}


bool Foam::biot::solveNewtonKrylov
(
    const volScalarField& impK,
    volScalarField& invM,
    const dimensionedScalar& invMsrc,
    const volScalarField::Internal& fGen,
    const label nSteps,
    const scalar pScale,
    label& nNewton,
    scalar& changeP,
    scalar& changeD
)
{
    const label n = mesh_.nCells();
    const scalar dt = mesh_.time().deltaTValue();
    const scalarField& V = mesh_.V();

    // Top faces: the displacement of the traction directions is unknown
    const fvPatch& topPatch = mesh_.boundary()[topPatchi_];
    const label nb = topPatch.size();
    const labelUList& topCells = topPatch.faceCells();
    const scalarField magSfTop(topPatch.magSf());
    const scalarField dcTop(topPatch.deltaCoeffs());
    const scalarField impKTop
    (
        impKScale_*law_->impK().boundaryField()[topPatchi_]
    );
    const auto& DpTop =
        refCast<const directionMixedFvPatchVectorField>
        (
            D_.boundaryField()[topPatchi_]
        );

    // Increment fields with homogeneous versions of the D and p conditions
    volVectorField dD("dD", D_);
    dD.primitiveFieldRef() = Zero;
    {
        auto& top =
            refCast<directionMixedFvPatchVectorField>
            (
                dD.boundaryFieldRef()[topPatchi_]
            );
        top.refValue() = Zero;
        top.refGrad() = Zero;
        forAll(dD.boundaryField(), patchi)
        {
            if (isA<fixedValueFvPatchVectorField>(dD.boundaryField()[patchi]))
            {
                dD.boundaryFieldRef()[patchi] == Zero;
            }
        }
    }
    volScalarField dp("dp", p_);
    dp.primitiveFieldRef() = Zero;
    forAll(dp.boundaryField(), patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(dp.boundaryField()[patchi]))
        {
            dp.boundaryFieldRef()[patchi] == Zero;
        }
    }

    const surfaceScalarField mobFace(nSteps*mobilityFace(mobility_));

    // Inner solves of the preconditioner tightened: the case's relative
    // tolerances (1-5 %) leave a preconditioner that varies from one
    // Krylov vector to the next by more than the residual GMRES is asked
    // to reach (in parallel the DIC solves stop at that tolerance)
    dictionary solverDictD(mesh_.solverDict(D_.name()));
    dictionary solverDictP(mesh_.solverDict(p_.name()));
    for (dictionary* d : {&solverDictD, &solverDictP})
    {
        d->set("relTol", 0.0);
        d->set("tolerance", 1e-12);
        d->set("maxIter", 2000);
    }

    // Preconditioner: one linear segregated sweep on the increment,
    // A_lap dD = r_D with the top-face block eliminated (the traction
    // defect r_b as the gradient increment of the traction faces through
    // the implicit stiffness, then the face displacement from the face
    // cell and that gradient), then (storage/dt - lap) dp =
    // r_p - div(dD)/dt; r in matrix units (integrated over the cells)
    auto precondition = [&](const nkVector& r, nkVector& z)
    {
        dD.primitiveFieldRef() = Zero;
        const symmTensorField freeDir
        (
            symmTensor::I - DpTop.valueFraction()
        );
        {
            auto& top =
                refCast<directionMixedFvPatchVectorField>
                (
                    dD.boundaryFieldRef()[topPatchi_]
                );
            top.valueFraction() = DpTop.valueFraction();
            top.refGrad() = (freeDir & r.b)/(impKTop*magSfTop);
        }
        dD.correctBoundaryConditions();
        {
            // fvMatrix convention: A psi = source (volume-integrated)
            fvVectorMatrix M(fvm::laplacian(impK, dD, "laplacian(DD,D)"));
            M.source() += r.D;
            if (maskedCells_.size())
            {
                M.setValues
                (
                    maskedCells_,
                    vectorField(maskedCells_.size(), Zero)
                );
            }
            M.solve(solverDictD);
        }
        z.D = dD.primitiveField();
        z.b =
            freeDir
          & (
                vectorField(z.D, topCells)
              + r.b/(impKTop*dcTop*magSfTop)
            );

        dp.primitiveFieldRef() = Zero;
        dp.correctBoundaryConditions();
        {
            volScalarField::Internal divdD
            (
                IOobject("divdD", mesh_.time().timeName(), mesh_),
                fvc::div(dD, "div(D)")()/dt
            );
            const dimensionedScalar dtDim("dt", dimTime, dt);
            const volScalarField::Internal coef
            (
                "nkStorage",
                (storage_ + invM)()/dtDim
            );
            fvScalarMatrix M
            (
                fvm::Sp(coef, dp)
              - fvm::laplacian(mobFace, dp)
            );
            M.source() += r.p - divdD.field()*V;
            M.solve(solverDictP);
        }
        z.p = dp.primitiveField();
    };

    // The state as a vector, and residual evaluation at a state
    auto getState = [&](nkVector& x)
    {
        x.D = D_.primitiveField();
        x.p = p_.primitiveField();
        x.b = D_.boundaryField()[topPatchi_];
    };
    auto setState = [&](const nkVector& x)
    {
        D_.primitiveFieldRef() = x.D;
        p_.primitiveFieldRef() = x.p;
    };
    // initial: refresh the traction condition from the state (the fixed
    // point of the corrector) instead of taking the face displacement
    auto residual = [&](const nkVector& x, nkVector& r, const bool initial)
    {
        setState(x);
        law_->unfreeze();
        coupledResidual
        (
            impK, invM, invMsrc, fGen, nSteps,
            initial ? nullptr : &x.b,
            r.D, r.p, r.b
        );
    };

    // Starting point: the traction condition refreshed once from the
    // committed state gives the value fractions and a first face
    // displacement, then the residual of that state
    nkVector x(n, nb), r(n, nb);
    getState(x);
    residual(x, r, true);
    getState(x);
    residual(x, r, false);
    scalar mismatch = tractionMismatch(pScale);

    // Block scales (the blocks have different units): the residuals by
    // the magnitude of a typical term of each equation over the mesh (a
    // pressure gradient pScale over a cell for D, the storage term of
    // pScale for p), the state blocks by their magnitude
    const scalar one = 1;
    const scalar cellSize = Foam::cbrt(gAverage(V));
    const scalar sRD = max(pScale*gSum(pow(V, 2.0/3.0)), SMALL);
    const scalar sRp =
        max(gSum((storage_ + invM)().primitiveField()*V)*pScale/dt, SMALL);
    const scalar sXD =
        max
        (
            gMax(mag(D_.primitiveField())),
            pScale*cellSize/max(gAverage(impK.primitiveField()), SMALL)
        );
    const scalar sXp = max(pScale, SMALL);
    const scalar sRb = max(pScale*gSum(magSfTop), SMALL);

    const scalar r0 = norm(r, sRD, sRp, sRb);
    scalar rNorm = r0;

    if (debug > 2)
    {
        // Consistency of the residual: the same state twice, and again
        // after a perturbation (any difference is state carried between
        // evaluations, which the finite-difference Jacobian cannot see)
        nkVector r2(n, nb), r3(n, nb), xq(x);
        residual(x, r2, false);
        xq.D *= 1.001;
        xq.p *= 1.001;
        xq.b *= 1.001;
        residual(xq, r3, false);
        residual(x, r3, false);
        nkVector d2(r2), d3(r3);
        axpy(d2, -one, r);
        axpy(d3, -one, r);
        // Where the residual sits: the blocks, and the cells next to
        // coupled (processor) patches
        scalarField rD2(magSqr(r.D)), rp2(sqr(r.p));
        boolList coupledCell(n, false);
        forAll(mesh_.boundary(), patchi)
        {
            if (mesh_.boundary()[patchi].coupled())
            {
                for (const label celli : mesh_.boundary()[patchi].faceCells())
                {
                    coupledCell[celli] = true;
                }
            }
        }
        scalar rD2c = 0, rp2c = 0;
        label nc = 0;
        forAll(coupledCell, celli)
        {
            if (coupledCell[celli])
            {
                rD2c += rD2[celli];
                rp2c += rp2[celli];
                ++nc;
            }
        }
        Info<< "        residual blocks: D " << Foam::sqrt(gSum(rD2))/sRD
            << " p " << Foam::sqrt(gSum(rp2))/sRp
            << " b " << Foam::sqrt(gSum(magSqr(r.b)))/sRb
            << "; next to coupled patches ("
            << returnReduce(nc, sumOp<label>())
            << " cells): D "
            << Foam::sqrt(returnReduce(rD2c, sumOp<scalar>()))/sRD << " p "
            << Foam::sqrt(returnReduce(rp2c, sumOp<scalar>()))/sRp
            << "; max |rD| " << gMax(mag(r.D)) << " max |rp| " << gMax(mag(r.p))
            << endl;
        // The pressure residual term by term, matrix against explicit
        {
            fvScalarMatrix lapM(fvm::laplacian(mobFace, p_));
            const scalarField rLapM
            (
                -(lapM & p_.internalField())().primitiveField()*V
            );
            const scalarField rLapC
            (
                fvc::laplacian(mobFace, p_)().primitiveField()*V
            );
            fvScalarMatrix ddtM(fvm::ddt(storage_ + invM, p_));
            const scalarField rDdtM
            (
                ddtM.source() - ddtM.diag()*p_.primitiveField()
            );
            const scalarField rDdtC
            (
                ((storage_ + invM)*fvc::ddt(p_))().primitiveField()*V
            );
            const scalarField rDivD(fvc::ddt(divD_)().primitiveField()*V);
            Info<< "        p terms /sRp: laplacian matrix "
                << Foam::sqrt(gSum(sqr(rLapM)))/sRp << " explicit "
                << Foam::sqrt(gSum(sqr(rLapC)))/sRp << " diff "
                << Foam::sqrt(gSum(sqr(rLapM + rLapC)))/sRp
                << "; ddt matrix " << Foam::sqrt(gSum(sqr(rDdtM)))/sRp
                << " explicit " << Foam::sqrt(gSum(sqr(rDdtC)))/sRp
                << " diff " << Foam::sqrt(gSum(sqr(rDdtM + rDdtC)))/sRp
                << "; ddt(divD) " << Foam::sqrt(gSum(sqr(rDivD)))/sRp
                << endl;
        }
        Info<< "        residual consistency: |R| " << r0
            << " repeat " << norm(d2, sRD, sRp, sRb)/r0
            << " after a perturbation " << norm(d3, sRD, sRp, sRb)/r0
            << " (D " << Foam::sqrt(gSum(magSqr(d3.D)))/sRD/r0
            << " p " << Foam::sqrt(gSum(sqr(d3.p)))/sRp/r0
            << " b " << Foam::sqrt(gSum(magSqr(d3.b)))/sRb/r0 << ")"
            << endl;
    }
    label iNewton = 0;
    label nResidual = 2;
    bool lineSearchFailed = false;

    const vectorField D0(D_.primitiveField());
    const scalarField p0(p_.primitiveField());

    for
    (
        ;
        iNewton < maxNewton_ && rNorm > newtonTol_*r0;
        ++iNewton
    )
    {
        // Right-preconditioned restarted GMRES on J M^-1 y = -r
        nkVector dx(n, nb);
        scalar gmresNorm = rNorm;

        for (label restart = 0; restart < krylovRestarts_; ++restart)
        {
            // Residual of the linear system at the current dx
            nkVector b(r);
            if (restart > 0)
            {
                nkVector xp(x);
                axpy(xp, one, dx);
                residual(xp, b, false);
                ++nResidual;
            }
            // b = -(r + J dx) approximated by the full residual at x + dx
            b.D = -b.D;
            b.p = -b.p;
            b.b = -b.b;
            const scalar beta = norm(b, sRD, sRp, sRb);
            if (beta < newtonTol_*r0)
            {
                gmresNorm = beta;
                break;
            }

            List<nkVector> Vk;
            Vk.setSize(krylovDim_ + 1, nkVector(n, nb));
            List<nkVector> Zk;
            Zk.setSize(krylovDim_, nkVector(n, nb));
            scalarSquareMatrix H(krylovDim_ + 1, Zero);
            scalarField cs(krylovDim_, Zero), sn(krylovDim_, Zero);
            scalarField g(krylovDim_ + 1, Zero);
            g[0] = beta;
            Vk[0].D = b.D/beta;
            Vk[0].p = b.p/beta;
            Vk[0].b = b.b/beta;

            label k = 0;
            for (; k < krylovDim_; ++k)
            {
                // z = M^-1 v, w = J z by finite differences
                precondition(Vk[k], Zk[k]);
                const scalar zn = norm(Zk[k], sXD, sXp, sXD);
                if (zn < VSMALL)
                {
                    break;
                }
                const scalar eps =
                    fdEpsilon_*(1 + norm(x, sXD, sXp, sXD))/zn;
                nkVector xp(x);
                axpy(xp, one, dx);
                axpy(xp, eps, Zk[k]);
                nkVector rp(n, nb);
                residual(xp, rp, false);
                ++nResidual;
                // The linear residual at x + dx is -b;
                // J z = (R(x + dx + eps z) - R(x + dx))/eps
                nkVector w(n, nb);
                w.D = (rp.D - (-b.D))/eps;
                w.p = (rp.p - (-b.p))/eps;
                w.b = (rp.b - (-b.b))/eps;

                if (debug > 2 && k < 3)
                {
                    const scalar cosang =
                        dot(w, Vk[k], sRD, sRp, sRb)
                       /max
                        (
                            norm(w, sRD, sRp, sRb)
                           *norm(Vk[k], sRD, sRp, sRb),
                            VSMALL
                        );
                    Info<< "        GMRES k " << k << ": |v| "
                        << norm(Vk[k], sRD, sRp, sRb) << " (D "
                        << Foam::sqrt(gSum(magSqr(Vk[k].D)))/sRD << ", p "
                        << Foam::sqrt(gSum(sqr(Vk[k].p)))/sRp << ") |z| " << zn
                        << " (D " << gMax(mag(Zk[k].D)) << " m, p "
                        << gMax(mag(Zk[k].p)) << " Pa) eps " << eps
                        << " |Jz| " << norm(w, sRD, sRp, sRb) << " (D "
                        << Foam::sqrt(gSum(magSqr(w.D)))/sRD << ", p "
                        << Foam::sqrt(gSum(sqr(w.p)))/sRp << ") cos "
                        << cosang << endl;
                }

                // Modified Gram-Schmidt
                for (label j = 0; j <= k; ++j)
                {
                    H(j, k) = dot(w, Vk[j], sRD, sRp, sRb);
                    axpy(w, -H(j, k), Vk[j]);
                }
                H(k + 1, k) = norm(w, sRD, sRp, sRb);
                const scalar hSub = H(k + 1, k);
                if (hSub > VSMALL)
                {
                    Vk[k + 1].D = w.D/hSub;
                    Vk[k + 1].p = w.p/hSub;
                    Vk[k + 1].b = w.b/hSub;
                }

                // Givens rotations
                for (label j = 0; j < k; ++j)
                {
                    const scalar t = cs[j]*H(j, k) + sn[j]*H(j + 1, k);
                    H(j + 1, k) = -sn[j]*H(j, k) + cs[j]*H(j + 1, k);
                    H(j, k) = t;
                }
                const scalar hk = Foam::sqrt(sqr(H(k, k)) + sqr(H(k + 1, k)));
                cs[k] = H(k, k)/max(hk, VSMALL);
                sn[k] = H(k + 1, k)/max(hk, VSMALL);
                H(k, k) = hk;
                H(k + 1, k) = 0;
                g[k + 1] = -sn[k]*g[k];
                g[k] = cs[k]*g[k];

                // Converged, or an invariant subspace (the rotation has
                // zeroed the sub-diagonal entry: test its value before)
                if (mag(g[k + 1]) < 0.01*newtonTol_*r0 || hSub < VSMALL)
                {
                    ++k;
                    break;
                }
            }

            // Back substitution and the update dx += sum y_j z_j
            scalarField y(k, Zero);
            for (label i = k - 1; i >= 0; --i)
            {
                scalar sum = g[i];
                for (label j = i + 1; j < k; ++j)
                {
                    sum -= H(i, j)*y[j];
                }
                y[i] = sum/H(i, i);
            }
            for (label j = 0; j < k; ++j)
            {
                axpy(dx, y[j], Zk[j]);
            }
            gmresNorm = mag(g[k]);

            if (debug > 1)
            {
                Info<< "        GMRES restart " << restart << ": " << k
                    << " vectors, residual " << beta/r0 << " -> "
                    << gmresNorm/r0 << endl;
            }

            if (gmresNorm < 0.01*newtonTol_*r0)
            {
                break;
            }
        }

        // Newton step with backtracking on the residual norm; a step
        // that does not decrease it is rejected (the linear model is no
        // guide across a kink of the law: the caller falls back)
        scalar lambda = 1;
        nkVector xNew(n, nb), rNew(n, nb);
        scalar rNew2 = GREAT;
        bool decreased = false;
        for (label ls = 0; ls < 6; ++ls)
        {
            xNew = x;
            axpy(xNew, lambda, dx);
            residual(xNew, rNew, false);
            ++nResidual;
            rNew2 = norm(rNew, sRD, sRp, sRb);
            if (rNew2 <= rNorm*(1 - 1e-4*lambda))
            {
                decreased = true;
                break;
            }
            lambda *= 0.5;
        }

        if (!decreased && rNew2 >= rNorm)
        {
            if (debug)
            {
                Info<< "    Newton " << iNewton + 1 << ": no decrease ("
                    << rNew2/r0 << " from " << rNorm/r0 << "), rejected"
                    << endl;
            }
            lineSearchFailed = true;
            break;
        }

        // A step that reduces the residual by less than 10 % means the
        // linear model is no guide (perfectly plastic law at the limit
        // load, where Newton stalls at a fixed fraction of the initial
        // residual): give the step to the segregated correctors at once
        // rather than after maxNewton such steps
        if (rNew2 > 0.9*rNorm)
        {
            if (debug)
            {
                Info<< "    Newton " << iNewton + 1 << ": residual only "
                    << rNew2/rNorm << " of the previous one, giving up"
                    << endl;
            }
            x = xNew;
            r = rNew;
            rNorm = rNew2;
            lineSearchFailed = true;
            ++iNewton;
            break;
        }

        if (debug)
        {
            Info<< "    Newton " << iNewton + 1 << ": |R| " << rNew2/r0
                << " (GMRES " << gmresNorm/r0 << ", step " << lambda
                << ", " << nResidual << " residuals)" << endl;
        }

        x = xNew;
        r = rNew;
        rNorm = rNew2;
    }

    // Leave the state consistent with the final iterate, the traction
    // condition in its incremental form (unchanged at convergence)
    setState(x);
    residual(x, r, false);
    ++nResidual;
    mismatch = tractionMismatch(pScale);
    updateTraction();
    D_.correctBoundaryConditions();

    changeD =
        gMax(mag(D_.primitiveField() - D0))
       /max(gMax(mag(D_.primitiveField())), SMALL);
    changeP = gMax(mag(p_.primitiveField() - p0))/pScale;

    nNewton = iNewton;
    const bool converged = (rNorm <= newtonTol_*r0);

    if (debug || !converged)
    {
        Info<< "    Biot Newton-Krylov: " << iNewton << " Newton steps, "
            << nResidual << " residual evaluations, |R|/|R0| " << rNorm/r0
            << ", traction mismatch " << mismatch << ", changes D "
            << changeD << " p " << changeP
            << (converged ? "" : lineSearchFailed
              ? " (line search failed)" : " (not converged)")
            << endl;
    }

    return converged;
}


// ************************************************************************* //
