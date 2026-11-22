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

#include "pastorZienkiewicz.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace soilConstitutiveModels
{
    defineTypeNameAndDebug(pastorZienkiewicz, 0);
    addToRunTimeSelectionTable
    (
        soilConstitutiveModel,
        pastorZienkiewicz,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::soilConstitutiveModels::pastorZienkiewicz::floored
(
    const scalar pRaw
) const
{
    const scalar x = pRaw - pMin_;
    const scalar w = floorWidth_;
    if (w <= 0)
    {
        return max(pRaw, pMin_);
    }
    if (x >= w)
    {
        return pRaw;
    }
    if (x <= -w)
    {
        return pMin_;
    }
    return pMin_ + sqr(x + w)/(4*w);
}


Foam::scalar Foam::soilConstitutiveModels::pastorZienkiewicz::zeta
(
    const scalar p,
    const scalar eta
) const
{
    const scalar base = 1 - alphaF_*eta/((1 + alphaF_)*Mf_);
    return p*Foam::pow(max(base, SMALL), -1/alphaF_);
}


Foam::label Foam::soilConstitutiveModels::pastorZienkiewicz::integrate
(
    const state& s0,
    const symmTensor& epsTotal,
    state& s,
    label& dir
) const
{
    s = s0;

    const symmTensor dEpsTotal(epsTotal - s0.eps);
    const scalar magDEps = mag(dEpsTotal);

    if (magDEps < SMALL)
    {
        s.eps = epsTotal;
        return 0;
    }

    // Sub-steps: the elastic stress increment of a sub-step at most
    // substepFraction x p'
    scalar p = floored(-tr(s.sigma)/3 + pT_);
    const scalar Kc = K0_*min(p, s.pCap)/pRef_;
    const scalar Gc = G0_*min(p, s.pCap)/pRef_;
    const scalar magDSigmaE =
        mag(Kc*tr(dEpsTotal)*I + 2*Gc*dev(dEpsTotal));
    const label nSub =
        min
        (
            max
            (
                max
                (
                    label(std::ceil(magDSigmaE/(substepFraction_*p))),
                    label(std::ceil(magDEps/substepStrain_))
                ),
                1
            ),
            maxSubsteps_
        );
    const symmTensor dEps(dEpsTotal/nSub);

    label nPlastic = 0;

    for (label sub = 0; sub < nSub; ++sub)
    {
        p = floored(-tr(s.sigma)/3 + pT_);
        const symmTensor sDev(dev(s.sigma));
        const scalar q = Foam::sqrt(1.5*(sDev && sDev));
        const scalar eta = q/p;

        // Stress ratio at which the loading modulus vanishes (Hf = 0) and
        // the mobilised-stress function is singular
        const scalar etaF = (1 + 1/alphaF_)*Mf_;

        const scalar K = K0_*min(p, s.pCap)/pRef_;
        const scalar G = G0_*min(p, s.pCap)/pRef_;

        // Elastic trial
        const symmTensor dSigmaE(K*tr(dEps)*I + 2*G*dev(dEps));

        // Deviatoric unit direction (3/(2q)) s; zero at q = 0
        const symmTensor nDev(q > SMALL ? (1.5/q)*sDev : symmTensor(Zero));

        // Loading direction in (p', q): (d_f, 1) normalised; the tensor
        // form (compression positive for the volumetric part)
        const scalar df = (1 + alphaF_)*(Mf_ - eta);
        const scalar nfv = df/Foam::sqrt(1 + sqr(df));
        const scalar nfs = 1/Foam::sqrt(1 + sqr(df));
        const symmTensor nf(-nfv/3*I + nfs*nDev);

        const scalar trial = (nf && dSigmaE);

        if (mag(trial) < SMALL)
        {
            s.sigma += dSigmaE;
            s.eps += dEps;
            continue;
        }

        // Loading direction: the frozen one of the step, else decided
        // here (and frozen for the correctors to come)
        if (dir == 0 || !freezeDirection_)
        {
            dir = (trial > 0 ? 1 : -1);
        }
        const bool loading = (dir > 0);

        // Flow direction: (d_g, 1) normalised, contractive in unloading
        const scalar dg = (1 + alphaG_)*(Mg_ - eta);
        const scalar ngv =
            (loading ? dg : -mag(dg))/Foam::sqrt(1 + sqr(dg));
        const scalar ngs = 1/Foam::sqrt(1 + sqr(dg));
        const symmTensor ng(-ngv/3*I + ngs*nDev);

        // Reversal bookkeeping
        if (loading && s.loading <= 0)
        {
            s.loading = 1;
        }
        else if (!loading && s.loading >= 0)
        {
            s.loading = -1;
            s.etaU = max(eta, SMALL);
        }

        // Plastic modulus
        scalar H;
        if (loading)
        {
            const scalar Hf = Foam::pow(max(1 - eta/etaF, scalar(0)), 4);
            const scalar Hv = 1 - eta/Mg_;
            const scalar Hs = beta0_*beta1_*Foam::exp(-beta0_*s.xi);
            // zeta is singular at eta = etaF, where Hf = 0 anyway: keep
            // the ratio below it so that zetaMax stays finite (an
            // overflowing HDM times Hf = 0 gave NaN at the free surface)
            const scalar z = zeta(p, min(eta, 0.99*etaF));
            s.zetaMax = max(s.zetaMax, z);
            const scalar HDM =
                Foam::pow(max(s.zetaMax, SMALL)/max(z, SMALL), gammaDM_);
            // The softening branch (Hv + Hs < 0 beyond the peak) is
            // floored at zero: the explicit integration from the
            // committed state sees the transient strains of the
            // segregated iterations, which a negative modulus turns into
            // a singular multiplier
            H = H0_*p*Hf*max(Hv + Hs, scalar(0))*HDM;
        }
        else
        {
            const scalar r = Mg_/max(s.etaU, SMALL);
            H = (mag(r) > 1 ? Hu0_*Foam::pow(mag(r), gammaU_) : Hu0_);
        }

        // Plastic multiplier (strain driven)
        const symmTensor DeNg(K*tr(ng)*I + 2*G*dev(ng));
        const scalar denom = H + (nf && DeNg);
        const scalar dLambda = (denom > SMALL ? trial/denom : 0);

        // The multiplier is positive in loading and negative in
        // unloading (n_f : dsigma^e < 0), where the contractive unloading
        // direction turns it into contraction and reversed shear
        symmTensor dEpsP(Zero);
        if ((loading && dLambda > 0) || (!loading && dLambda < 0))
        {
            dEpsP = dLambda*ng;
            ++nPlastic;
        }

        s.sigma += K*tr(dEps - dEpsP)*I + 2*G*dev(dEps - dEpsP);
        s.eps += dEps;
        s.epsP += dEpsP;

        // Accumulated deviatoric plastic strain
        const symmTensor eP(dev(dEpsP));
        s.xi += Foam::sqrt(2.0/3.0*(eP && eP));

        // Tension cut-off: keep p' at the floor (C1 with floorWidth)
        const scalar pNew = -tr(s.sigma)/3 + pT_;
        const scalar pFloored = floored(pNew);
        if (pFloored > pNew)
        {
            s.sigma += (pNew - pFloored)*symmTensor::I;
        }

        // Stress-ratio cut-off at the failure line eta_F (the law has no
        // strength beyond it; reached only at the floor, where p' is
        // raised but q is not)
        {
            const scalar pS = -tr(s.sigma)/3 + pT_;
            const symmTensor sD(dev(s.sigma));
            const scalar qS = Foam::sqrt(1.5*(sD && sD));
            if (qS > etaF*pS && qS > SMALL)
            {
                s.sigma += (etaF*pS/qS - 1)*sD;
            }
        }

        if (!std::isfinite(mag(s.sigma)) || !std::isfinite(s.xi))
        {
            FatalErrorInFunction
                << "Non-finite state after sub-step " << sub << " of " << nSub
                << ": p " << p << " q " << q << " eta " << eta << " H " << H
                << " denom " << denom << " dLambda " << dLambda
                << " trial " << trial << " dEps " << dEps << " dEpsTotal "
                << dEpsTotal << " sigma0 " << s0.sigma << " loading "
                << s.loading << " etaU " << s.etaU << " xi " << s.xi
                << " zetaMax " << s.zetaMax << exit(FatalError);
        }
    }

    return nPlastic;
}


void Foam::soilConstitutiveModels::pastorZienkiewicz::initialise
(
    const volSymmTensorField& sigma0
)
{
    const symmTensorField& s0 = sigma0.primitiveField();

    cellOld_.resize(s0.size());
    forAll(s0, celli)
    {
        state& s = cellOld_[celli];
        s = state();
        s.sigma = s0[celli];
        const scalar p = max(-tr(s0[celli])/3 + pT_, pMin_);
        const symmTensor sDev(dev(s0[celli]));
        const scalar q = Foam::sqrt(1.5*(sDev && sDev));
        s.zetaMax = zeta(p, q/p);
        s.pCap = modulusCap_*p;
    }
    cell_ = cellOld_;

    patchOld_.resize(sigma0.boundaryField().size());
    forAll(sigma0.boundaryField(), patchi)
    {
        const symmTensorField& s0B = sigma0.boundaryField()[patchi];
        patchOld_[patchi].resize(s0B.size());
        forAll(s0B, facei)
        {
            state& s = patchOld_[patchi][facei];
            s = state();
            s.sigma = s0B[facei];
            const scalar p = max(-tr(s0B[facei])/3 + pT_, pMin_);
            const symmTensor sDev(dev(s0B[facei]));
            const scalar q = Foam::sqrt(1.5*(sDev && sDev));
            s.zetaMax = zeta(p, q/p);
            s.pCap = modulusCap_*p;
        }
    }
    patch_ = patchOld_;

    // Implicit operator: the elastic moduli of the geostatic p'
    scalarField& mu = mu_.primitiveFieldRef();
    scalarField& lambda = lambda_.primitiveFieldRef();
    scalarField& impK = impK_.primitiveFieldRef();
    forAll(mu, celli)
    {
        // Implicit operator at the cap of the moduli
        const scalar p = modulusCap_*max(-tr(s0[celli])/3 + pT_, pMin_);
        const scalar K = K0_*p/pRef_;
        const scalar G = G0_*p/pRef_;
        mu[celli] = G;
        lambda[celli] = K - 2*G/3;
        impK[celli] = K + 4*G/3;
    }
    mu_.correctBoundaryConditions();
    lambda_.correctBoundaryConditions();
    impK_.correctBoundaryConditions();

    initialised_ = true;

    Info<< "    Pastor-Zienkiewicz: states initialised, p' "
        << gMin(-tr(s0)/3) << " to " << gMax(-tr(s0)/3) << " Pa" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::soilConstitutiveModels::pastorZienkiewicz::pastorZienkiewicz
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    soilConstitutiveModel(mesh, dict),
    K0_(dict.getOrDefault<scalar>("K0", 35e6)),
    G0_(dict.getOrDefault<scalar>("G0", 52.5e6)),
    pRef_(dict.getOrDefault<scalar>("pRef", 100e3)),
    Mg_(dict.getOrDefault<scalar>("Mg", 1.5)),
    Mf_(dict.getOrDefault<scalar>("Mf", 1.1)),
    alphaF_(dict.getOrDefault<scalar>("alphaF", 0.45)),
    alphaG_(dict.getOrDefault<scalar>("alphaG", 0.45)),
    H0_(dict.getOrDefault<scalar>("H0", 350)),
    beta0_(dict.getOrDefault<scalar>("beta0", 4.2)),
    beta1_(dict.getOrDefault<scalar>("beta1", 0.2)),
    gammaDM_(dict.getOrDefault<scalar>("gammaDM", 2)),
    Hu0_(dict.getOrDefault<scalar>("Hu0", 2e6)),
    gammaU_(dict.getOrDefault<scalar>("gammaU", 2)),
    pMin_(dict.getOrDefault<scalar>("minPressure", 1000)),
    pT_(dict.getOrDefault<scalar>("tensileStrength", 0)),
    floorWidth_(dict.getOrDefault<scalar>("floorWidth", 0)),
    modulusCap_(dict.getOrDefault<scalar>("modulusCap", 1.5)),
    substepFraction_(dict.getOrDefault<scalar>("substepFraction", 0.02)),
    substepStrain_(dict.getOrDefault<scalar>("substepStrain", 1e-4)),
    maxSubsteps_(dict.getOrDefault<label>("maxSubsteps", 500)),
    freezeDirection_(dict.getOrDefault<bool>("freezeDirection", true)),
    mu_
    (
        IOobject("soilMu", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(dimPressure, G0_*pMin_/pRef_),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    lambda_
    (
        IOobject("soilLambda", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(dimPressure, (K0_ - 2*G0_/3)*pMin_/pRef_),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    impK_
    (
        IOobject("soilImpK", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(dimPressure, (K0_ + 4*G0_/3)*pMin_/pRef_),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    cellOld_(),
    cell_(),
    patchOld_(),
    patch_(),
    cellDir_(),
    patchDir_(),
    pPrime_
    (
        IOobject
        (
            "pPrime",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            dict.getOrDefault<bool>("writeState", true)
          ? IOobject::AUTO_WRITE
          : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPressure, Zero)
    ),
    qDev_
    (
        IOobject
        (
            "qDev",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            dict.getOrDefault<bool>("writeState", true)
          ? IOobject::AUTO_WRITE
          : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPressure, Zero)
    ),
    xiField_
    (
        IOobject
        (
            "xiPlastic",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            dict.getOrDefault<bool>("writeState", true)
          ? IOobject::AUTO_WRITE
          : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    ),
    nPlastic_(0),
    initialised_(false),
    geostaticSeen_(false)
{
    Info<< "    Pastor-Zienkiewicz soil: K0 " << K0_ << ", G0 " << G0_
        << " Pa at p_ref " << pRef_ << " Pa; Mg " << Mg_ << ", Mf " << Mf_
        << ", alpha " << alphaF_ << "/" << alphaG_ << ", H0 " << H0_
        << ", beta0 " << beta0_ << ", beta1 " << beta1_ << ", gammaDM "
        << gammaDM_ << ", Hu0 " << Hu0_ << ", gammaU " << gammaU_
        << ", p' floor " << pMin_ << " Pa (C1 width " << floorWidth_
        << " Pa), tensile strength " << pT_ << " Pa" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::soilConstitutiveModels::pastorZienkiewicz::geostatic
(
    const volSymmTensorField& sigma0
)
{
    // Once: the committed states carry the loading history afterwards
    if (!geostaticSeen_)
    {
        initialise(sigma0);
        geostaticSeen_ = true;
    }
}


void Foam::soilConstitutiveModels::pastorZienkiewicz::correct
(
    const volTensorField& gradD,
    const volSymmTensorField& sigma0,
    volSymmTensorField& sigma
)
{
    if (!initialised_)
    {
        initialise(sigma0);
    }

    const symmTensorField eps(symm(gradD.primitiveField()));
    const symmTensorField& s0 = sigma0.primitiveField();
    symmTensorField& sig = sigma.primitiveFieldRef();

    scalarField& pP = pPrime_.primitiveFieldRef();
    scalarField& qD = qDev_.primitiveFieldRef();
    scalarField& xi = xiField_.primitiveFieldRef();

    nPlastic_ = 0;

    if (cellDir_.size() != sig.size())
    {
        cellDir_.resize(sig.size(), 0);
    }

    forAll(sig, celli)
    {
        if (integrate(cellOld_[celli], eps[celli], cell_[celli], cellDir_[celli]) > 0)
        {
            ++nPlastic_;
        }
        sig[celli] = cell_[celli].sigma - s0[celli];

        pP[celli] = -tr(cell_[celli].sigma)/3;
        const symmTensor sDev(dev(cell_[celli].sigma));
        qD[celli] = Foam::sqrt(1.5*(sDev && sDev));
        xi[celli] = cell_[celli].xi;
    }

    // Boundary faces: their own states along the patch strain
    auto& sigBf = sigma.boundaryFieldRef();
    if (patchDir_.size() != sigBf.size())
    {
        patchDir_.resize(sigBf.size());
    }
    forAll(sigBf, patchi)
    {
        const symmTensorField epsB(symm(gradD.boundaryField()[patchi]));
        const symmTensorField& s0B = sigma0.boundaryField()[patchi];
        symmTensorField& sigB = sigBf[patchi];
        if (patchDir_[patchi].size() != sigB.size())
        {
            patchDir_[patchi].resize(sigB.size(), 0);
        }

        forAll(sigB, facei)
        {
            integrate
            (
                patchOld_[patchi][facei],
                epsB[facei],
                patch_[patchi][facei],
                patchDir_[patchi][facei]
            );
            sigB[facei] = patch_[patchi][facei].sigma - s0B[facei];
        }
    }

    pPrime_.correctBoundaryConditions();
    qDev_.correctBoundaryConditions();
    xiField_.correctBoundaryConditions();
}


Foam::tmp<Foam::volScalarField>
Foam::soilConstitutiveModels::pastorZienkiewicz::stabilisationModulusField() const
{
    auto tM = tmp<volScalarField>::New
    (
        IOobject("stabilisationModulus", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimPressure, K0_*pMin_/pRef_),
        zeroGradientFvPatchField<scalar>::typeName
    );
    scalarField& M = tM.ref().primitiveFieldRef();

    if (cell_.size() == M.size())
    {
        forAll(M, celli)
        {
            const state& s = cell_[celli];
            const scalar p =
                min(max(-tr(s.sigma)/3 + pT_, pMin_), s.pCap);
            M[celli] = K0_*p/pRef_;
        }
    }
    tM.ref().correctBoundaryConditions();

    return tM;
}


void Foam::soilConstitutiveModels::pastorZienkiewicz::unfreeze()
{
    cellDir_ = 0;
    forAll(patchDir_, patchi)
    {
        patchDir_[patchi] = 0;
    }
}


void Foam::soilConstitutiveModels::pastorZienkiewicz::beginStep()
{
    unfreeze();

    if (!geostaticSeen_)
    {
        return;
    }

    // Implicit operator and modulus cap of the step from the committed
    // p' of each point: the law's elastic stiffness stays at 1/modulusCap
    // of the implicit one, so the segregated correction keeps the same
    // contraction as the skeleton softens (with the geostatic cap the
    // ratio falls with p' and the correctors stall)
    scalarField& mu = mu_.primitiveFieldRef();
    scalarField& lambda = lambda_.primitiveFieldRef();
    scalarField& impK = impK_.primitiveFieldRef();

    forAll(cellOld_, celli)
    {
        state& s = cellOld_[celli];
        const scalar p = max(-tr(s.sigma)/3 + pT_, pMin_);
        s.pCap = modulusCap_*p;
        const scalar K = K0_*s.pCap/pRef_;
        const scalar G = G0_*s.pCap/pRef_;
        mu[celli] = G;
        lambda[celli] = K - 2*G/3;
        impK[celli] = K + 4*G/3;
    }
    forAll(patchOld_, patchi)
    {
        forAll(patchOld_[patchi], facei)
        {
            state& s = patchOld_[patchi][facei];
            s.pCap = modulusCap_*max(-tr(s.sigma)/3 + pT_, pMin_);
        }
    }
    cell_ = cellOld_;
    patch_ = patchOld_;

    mu_.correctBoundaryConditions();
    lambda_.correctBoundaryConditions();
    impK_.correctBoundaryConditions();
}


void Foam::soilConstitutiveModels::pastorZienkiewicz::updateState()
{
    cellOld_ = cell_;
    patchOld_ = patch_;

    const label n = returnReduce(nPlastic_, sumOp<label>());
    if (n > 0)
    {
        Info<< "    Pastor-Zienkiewicz: " << n << " plastic cells, p' "
            << gMin(pPrime_.primitiveField()) << " to "
            << gMax(pPrime_.primitiveField()) << " Pa, max xi "
            << gMax(xiField_.primitiveField()) << endl;
    }
}


void Foam::soilConstitutiveModels::pastorZienkiewicz::remap
(
    const labelList& c0,
    const labelList& c1,
    const scalarField& w
)
{
    // Cell i takes (1 - w) of old cell c0[i] and w of c1[i] in
    // soilModel::remapField; the state (path dependent) is taken from the
    // heavier source
    List<state> old(cellOld_);
    forAll(cellOld_, i)
    {
        cellOld_[i] = old[w[i] < 0.5 ? c0[i] : c1[i]];
    }
    cell_ = cellOld_;
}


// ************************************************************************* //
