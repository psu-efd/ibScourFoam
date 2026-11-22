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
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Vector of the block-coupled linear solve: the displacement increment
// (3 per cell) and the pressure increment, with block scales for the
// norms and dot products
struct bcVector
{
    Foam::vectorField D;
    Foam::scalarField p;

    bcVector() = default;

    explicit bcVector(const Foam::label n)
    :
        D(n, Foam::Zero), p(n, Foam::Zero)
    {}
};


inline Foam::scalar bcDot
(
    const bcVector& a,
    const bcVector& b,
    const Foam::scalar sD,
    const Foam::scalar sp
)
{
    return
        Foam::gSum(a.D & b.D)/Foam::sqr(sD)
      + Foam::gSum(a.p*b.p)/Foam::sqr(sp);
}


inline Foam::scalar bcNorm
(
    const bcVector& a,
    const Foam::scalar sD,
    const Foam::scalar sp
)
{
    return Foam::sqrt(bcDot(a, a, sD, sp));
}


inline void bcAxpy(bcVector& y, const Foam::scalar a, const bcVector& x)
{
    y.D += a*x.D;
    y.p += a*x.p;
}

} // End namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::biot::solveBlockCoupled
(
    const volScalarField& impK,
    volScalarField& invM,
    const dimensionedScalar& invMsrc,
    const volScalarField::Internal& fGen,
    const label nSteps,
    const scalar pScale,
    label& nOuter,
    scalar& changeP,
    scalar& changeD
)
{
    const label n = mesh_.nCells();
    const scalar dt = mesh_.time().deltaTValue();
    const scalarField& V = mesh_.V();

    // The linear operator takes the volumetric strain rate as
    // div(dD)/dt: Euler time derivative, no volumetric smoothing
    if (volumetricSmoothing_ > 0)
    {
        FatalIOErrorInFunction(dict_)
            << "blockCoupled needs volumetricSmoothing 0"
            << exit(FatalIOError);
    }

    const auto& DpTop =
        refCast<const directionMixedFvPatchVectorField>
        (
            D_.boundaryField()[topPatchi_]
        );

    // Increment fields with homogeneous versions of the D and p conditions
    // (the traction gradient of the top faces is lagged over a linear
    // solve: their increment follows the cell)
    volVectorField dD
    (
        IOobject("dD", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedVector(D_.dimensions(), Zero),
        D_.boundaryField().types()
    );
    {
        auto& top =
            refCast<directionMixedFvPatchVectorField>
            (
                dD.boundaryFieldRef()[topPatchi_]
            );
        top.valueFraction() = DpTop.valueFraction();
        top.refValue() = Zero;
        top.refGrad() = Zero;
        dD.correctBoundaryConditions();
    }
    volScalarField dp
    (
        IOobject("dp", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(p_.dimensions(), Zero),
        p_.boundaryField().types()
    );
    dp.oldTime();

    const surfaceScalarField mobFace(nSteps*mobilityFace(mobility_));
    const dimensionedScalar dtDim("dt", dimTime, dt);

    // Homogeneous operators of the increment (their sources vanish):
    // A_D = laplacian(impK), A_p = storage/dt - laplacian(mobility)
    fvVectorMatrix AD(fvm::laplacian(impK, dD, "laplacian(DD,D)"));
    const volScalarField::Internal storageCoef
    (
        IOobject("bcStorage", mesh_.time().timeName(), mesh_),
        mesh_,
        storage_/dtDim
    );
    fvScalarMatrix Ap
    (
        fvm::Sp(storageCoef, dp)
      - fvm::laplacian(mobFace, dp)
    );

    // Traction directions of the top faces
    const symmTensorField freeDirTop
    (
        symmTensor::I - DpTop.valueFraction()
    );

    boolList masked(n, false);
    for (const label celli : maskedCells_)
    {
        masked[celli] = true;
    }

    // L x = [div(sigma(dD)) - grad dp ; A_p dp + div(dD)/dt] (matrix
    // units; the stress operator as the implicit Laplacian A_D plus its
    // explicit correction), the structure cells keep zero displacement
    // increments
    auto applyL = [&](const bcVector& x, bcVector& y)
    {
        dD.primitiveFieldRef() = x.D;
        dD.correctBoundaryConditions();
        dp.primitiveFieldRef() = x.p;
        dp.correctBoundaryConditions();

        // Full linearised stress: the implicit Laplacian plus the explicit
        // correction div(sigma(dD) - impK grad(dD)) with the law's tangent
        // moduli, formed as updateStress forms it
        const volTensorField gradDd(fvc::grad(dD));
        const volSymmTensorField sigmaD
        (
            law_->mu()*twoSymm(gradDd) + law_->lambda()*(I*tr(gradDd))
        );
        const volScalarField impKExp("impKExp", impKScale_*law_->impK());
        volVectorField expl
        (
            IOobject("bcExpl", mesh_.time().timeName(), mesh_),
            fvc::div(sigmaD, "div(sigmaD)")
        );
        if (stabilisationScale_ < 1)
        {
            expl -=
                stabilisationScale_
               *fvc::div(impKExp*gradDd, "div(sigmaD)")
              + (1 - stabilisationScale_)
               *fvc::laplacian(impKExp, dD, "laplacian(DD,D)");
        }
        else
        {
            expl -= fvc::div(impKExp*gradDd, "div(sigmaD)");
        }

        y.D =
            (AD & dD)().primitiveField()*V
          + expl.primitiveField()*V
          - fvc::grad(dp)().primitiveField()*V;

        // Linearised traction condition of the top faces: the effective
        // stress flux through a traction face is the prescribed traction,
        // its increment vanishes in the traction directions (the face
        // displacement adjusts); remove the explicit face flux there
        {
            const fvPatch& tp = mesh_.boundary()[topPatchi_];
            const vectorField& Sf = tp.Sf();
            const labelUList& fc = tp.faceCells();
            const symmTensorField& sb = sigmaD.boundaryField()[topPatchi_];
            const tensorField& gb = gradDd.boundaryField()[topPatchi_];
            const scalarField& kb = impKExp.boundaryField()[topPatchi_];
            const scalar sStab = min(stabilisationScale_, scalar(1));
            forAll(fc, i)
            {
                const vector F = (Sf[i] & sb[i]) - sStab*kb[i]*(Sf[i] & gb[i]);
                y.D[fc[i]] -= (freeDirTop[i] & F);
            }
        }
        y.p =
            (Ap & dp)().primitiveField()*V
          + fvc::div(dD, "div(D)")().primitiveField()*V/dt;

        forAll(masked, celli)
        {
            if (masked[celli])
            {
                y.D[celli] = x.D[celli];
            }
        }
    };

    // Inner solves of the preconditioner to blockInnerRelTol: the GMRES
    // keeps the preconditioned vectors (flexible GMRES), so a
    // preconditioner that varies between Krylov vectors is admissible
    dictionary solverDictD(mesh_.solverDict(D_.name()));
    dictionary solverDictP(mesh_.solverDict(p_.name()));
    for (dictionary* d : {&solverDictD, &solverDictP})
    {
        d->set("relTol", blockInnerRelTol_);
        d->set("tolerance", 1e-14);
        d->set("maxIter", 2000);
    }

    // Preconditioner: one fixed-stress sweep, A_D zD = r_D, then
    // (A_p + 1/(K dt)) zp = r_p - div(zD)/dt
    auto precondition = [&](const bcVector& r, bcVector& z)
    {
        dD.primitiveFieldRef() = Zero;
        dD.correctBoundaryConditions();
        {
            // fvMatrix convention: A psi = source (volume-integrated)
            fvVectorMatrix M(AD);
            M.source() = r.D;
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

        dp.primitiveFieldRef() = Zero;
        dp.correctBoundaryConditions();
        {
            fvScalarMatrix M(Ap);
            M.diag() += invM.primitiveField()*V/dt;
            M.source() =
                r.p - fvc::div(dD, "div(D)")().primitiveField()*V/dt;
            M.solve(solverDictP);
        }
        z.p = dp.primitiveField();
    };

    vectorField rD;
    scalarField rp;
    vectorField rb;

    // Block scales of the residual (the blocks have different units):
    // a pressure gradient pScale over a cell for D, the storage of pScale
    // over the step for p
    const scalar hCell = Foam::cbrt(gAverage(V));
    const scalar sqrtN =
        Foam::sqrt(scalar(returnReduce(n, sumOp<label>())));
    const scalar sRD = max(pScale*sqr(hCell)*sqrtN, VSMALL);
    const scalar sRp =
        max
        (
            pScale*(storage_.value() + gAverage(invM.primitiveField()))
           *gAverage(V)/dt*sqrtN,
            VSMALL
        );

    bool converged = false;
    changeP = GREAT;
    changeD = GREAT;
    scalar rNorm0 = 0;

    for (nOuter = 0; nOuter < nCorrectors_; ++nOuter)
    {
        // Traction condition of the top faces converged for the current
        // cells (a single refresh contracts slowly in a soft soil, by
        // about lambda/(lambda + 2 mu) per refresh), so that the residual
        // is a function of the cell state
        for (label k = 0; k < 50; ++k)
        {
            updateStress();
            updateTraction();
            D_.correctBoundaryConditions();
            if (tractionMismatch(pScale) < 1e-3*tolerance_)
            {
                break;
            }
        }

        // Residual of the current state: the law's stress, the traction
        // condition and the explicit stress terms refreshed
        coupledResidual
        (
            impK, invM, invMsrc, fGen, nSteps, nullptr, rD, rp, rb
        );
        bcVector r(n);
        r.D = rD;
        r.p = rp;
        const scalar rNorm = bcNorm(r, sRD, sRp);
        if (nOuter == 0)
        {
            rNorm0 = max(rNorm, VSMALL);
        }

        if
        (
            nOuter > 0
         && (max(changeP, changeD) < tolerance_ || rNorm < 1e-12*rNorm0)
        )
        {
            converged = true;
            break;
        }

        // Flexible (right-preconditioned) restarted GMRES on L x = r
        bcVector x(n);
        const scalar linTol = 1e-3;
        for (label restart = 0; restart < krylovRestarts_; ++restart)
        {
            bcVector b(r);
            if (restart > 0)
            {
                bcVector Lx(n);
                applyL(x, Lx);
                bcAxpy(b, -1, Lx);
            }
            const scalar beta = bcNorm(b, sRD, sRp);
            if (beta < linTol*rNorm)
            {
                break;
            }

            List<bcVector> Vk(krylovDim_ + 1, bcVector(n));
            List<bcVector> Zk(krylovDim_, bcVector(n));
            scalarSquareMatrix H(krylovDim_ + 1, Zero);
            scalarField cs(krylovDim_, Zero);
            scalarField sn(krylovDim_, Zero);
            scalarField g(krylovDim_ + 1, Zero);
            g[0] = beta;
            Vk[0].D = b.D/beta;
            Vk[0].p = b.p/beta;

            label k = 0;
            for (; k < krylovDim_; ++k)
            {
                precondition(Vk[k], Zk[k]);
                bcVector w(n);
                applyL(Zk[k], w);

                // Modified Gram-Schmidt
                for (label j = 0; j <= k; ++j)
                {
                    H(j, k) = bcDot(w, Vk[j], sRD, sRp);
                    bcAxpy(w, -H(j, k), Vk[j]);
                }
                H(k + 1, k) = bcNorm(w, sRD, sRp);
                const scalar hSub = H(k + 1, k);
                if (hSub > VSMALL)
                {
                    Vk[k + 1].D = w.D/hSub;
                    Vk[k + 1].p = w.p/hSub;
                }

                // Givens rotations
                for (label j = 0; j < k; ++j)
                {
                    const scalar t = cs[j]*H(j, k) + sn[j]*H(j + 1, k);
                    H(j + 1, k) = -sn[j]*H(j, k) + cs[j]*H(j + 1, k);
                    H(j, k) = t;
                }
                const scalar hk =
                    Foam::sqrt(sqr(H(k, k)) + sqr(H(k + 1, k)));
                cs[k] = H(k, k)/max(hk, VSMALL);
                sn[k] = H(k + 1, k)/max(hk, VSMALL);
                H(k, k) = hk;
                H(k + 1, k) = 0;
                g[k + 1] = -sn[k]*g[k];
                g[k] = cs[k]*g[k];

                if (mag(g[k + 1]) < linTol*rNorm || hSub < VSMALL)
                {
                    ++k;
                    break;
                }
            }

            // Back substitution and x += sum y_j z_j
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
                bcAxpy(x, y[j], Zk[j]);
            }

            if (debug > 1)
            {
                Info<< "        block GMRES restart " << restart << ": " << k
                    << " vectors, residual " << beta/rNorm << " -> "
                    << mag(g[k])/rNorm << endl;
            }

            if (mag(g[k]) < linTol*rNorm)
            {
                break;
            }
        }

        // Update of the state
        vectorField& D = D_.primitiveFieldRef();
        scalarField& p = p_.primitiveFieldRef();
        D += x.D;
        p += x.p;
        D_.correctBoundaryConditions();
        p_.correctBoundaryConditions();

        changeD = gMax(mag(x.D))/max(gMax(mag(D)), SMALL);
        changeP = gMax(mag(x.p))/pScale;

        if (debug)
        {
            Info<< "    Biot block-coupled outer " << nOuter << ": residual "
                << rNorm/rNorm0 << ", changes D " << changeD << " p "
                << changeP << endl;
        }
    }

    // The stress and the explicit terms of the final state
    updateStress();
    invM =
        dimensionedScalar(storage_.dimensions()*dimPressure, 1)
       /law_->stabilisationModulusField();

    return converged;
}


// ************************************************************************* //
