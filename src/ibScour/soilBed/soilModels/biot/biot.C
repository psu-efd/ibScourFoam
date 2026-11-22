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
#include "syncTools.H"
#include "poreGenerationModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvc.H"
#include "fixedValueFvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "slipFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(biot, 0);
    addToRunTimeSelectionTable(soilModel, biot, dictionary);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::biot::readBulkModulus(const dictionary& dict)
{
    if (dict.found("bulkModulus"))
    {
        return dict.get<scalar>("bulkModulus");
    }

    const scalar Sr = dict.getOrDefault<scalar>("saturation", 1);
    const scalar pAbs = dict.getOrDefault<scalar>("absolutePressure", 1.2e5);
    const scalar Kw = dict.getOrDefault<scalar>("waterBulkModulus", 2.2e9);

    return 1/(1/Kw + (1 - Sr)/pAbs);
}


Foam::wordList Foam::biot::pTypes() const
{
    const fvBoundaryMesh& bm = mesh_.boundary();

    wordList types(bm.size(), zeroGradientFvPatchScalarField::typeName);

    forAll(bm, patchi)
    {
        if (bm[patchi].coupled())
        {
            types[patchi] = bm[patchi].type();
        }
    }

    types[topPatchi_] =
        sealedTop_
      ? zeroGradientFvPatchScalarField::typeName
      : fixedValueFvPatchScalarField::typeName;

    if (bottomPressure_)
    {
        types[bottomPatchi_] = fixedValueFvPatchScalarField::typeName;
    }

    return types;
}


Foam::wordList Foam::biot::DTypes(const dictionary& dict) const
{
    const fvBoundaryMesh& bm = mesh_.boundary();
    const wordList fixed
    (
        dict.getOrDefault<wordList>("fixedDisplacementPatches", wordList())
    );

    wordList types(bm.size(), slipFvPatchVectorField::typeName);

    forAll(bm, patchi)
    {
        if (bm[patchi].coupled())
        {
            types[patchi] = bm[patchi].type();
        }
        else if (fixed.found(bm[patchi].name()))
        {
            types[patchi] = fixedValueFvPatchVectorField::typeName;
        }
    }

    types[topPatchi_] = directionMixedFvPatchVectorField::typeName;
    types[bottomPatchi_] = fixedValueFvPatchVectorField::typeName;

    return types;
}


void Foam::biot::updateTraction()
{
    // Traction condition (as tractionDisplacementFvPatchVectorField): the
    // gradient reproducing the effective traction t with the implicit
    // stiffness, corrected by the explicit stress
    // Faces with a prescribed displacement (value fraction 1) are fixed
    // instead
    auto& Dp =
        refCast<directionMixedFvPatchVectorField>
        (
            D_.boundaryFieldRef()[topPatchi_]
        );

    const fvPatch& patch = mesh_.boundary()[topPatchi_];
    const vectorField n(patch.nf());
    const scalarField impK
    (
        impKScale_*law_->impK().boundaryField()[topPatchi_]
    );
    const symmTensorField sigmaP(sigma_.boundaryField()[topPatchi_]);

    vectorField t(patch.size(), Zero);
    if (topTraction_.size() == patch.size())
    {
        t = topTraction_;
    }

    if (sealedTop_ && pTop_.size() == patch.size())
    {
        // Constant total traction: the effective traction increment is
        // the excess pore pressure of the top cells along the normal
        t += n*(p_.boundaryField()[topPatchi_].patchInternalField() - pTop_);
    }

    if (topFraction_.size() == patch.size())
    {
        // All components, or the tangential ones only (the normal
        // component keeps the traction condition)
        symmTensorField vf(patch.size(), symmTensor::I);
        if (topTangential_)
        {
            vf = symmTensor::I - sqr(n);
        }
        Dp.valueFraction() = topFraction_*vf;
        Dp.refValue() = topDisplacement_;
    }
    else
    {
        Dp.valueFraction() = Zero;
        Dp.refValue() = Zero;
    }

    // Gradient of the traction faces (the mixed snGrad is refGrad there)
    const vectorField snGradD(Dp.snGrad());
    Dp.refGrad() = (t + impK*snGradD - (n & sigmaP))/impK;
}


Foam::tmp<Foam::vectorField> Foam::biot::patchEffectiveTraction
(
    const label patchi
) const
{
    const fvPatch& patch = mesh_.boundary()[patchi];
    const fvPatchVectorField& Dp = D_.boundaryField()[patchi];

    const symmTensorField sigmaFace
    (
        sigma0_.boundaryField()[patchi] + sigma_.boundaryField()[patchi]
    );
    const symmTensorField sigmaCell
    (
        sigma0_.boundaryField()[patchi].patchInternalField()
      + sigma_.boundaryField()[patchi].patchInternalField()
    );

    // Fraction of the faces with a prescribed displacement
    scalarField f(patch.size(), Dp.fixesValue() ? 1 : 0);
    if (isA<directionMixedFvPatchVectorField>(Dp))
    {
        f = tr
        (
            refCast<const directionMixedFvPatchVectorField>(Dp)
           .valueFraction()
        )/3;
    }

    return patch.nf() & (f*sigmaCell + (1 - f)*sigmaFace);
}


Foam::vector Foam::biot::structureForce() const
{
    vector F(Zero);

    // All ranks take part in the processor exchange and the reduction,
    // also those without structure cells
    if (cellMask_.size() != mesh_.nCells())
    {
        reduce(F, sumOp<vector>());
        return F;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const vectorField& Sf = mesh_.faceAreas();
    const symmTensorField sigmaCell
    (
        sigma0_.primitiveField() + sigma_.primitiveField()
    );

    forAll(nei, facei)
    {
        const bool mo = cellMask_[own[facei]] > 0.5;
        const bool mn = cellMask_[nei[facei]] > 0.5;

        if (mo == mn)
        {
            continue;
        }

        // Sf points from the owner to the neighbour: out of the structure
        // when the owner is the structure
        if (mo)
        {
            F += sigmaCell[nei[facei]] & Sf[facei];
        }
        else
        {
            F -= sigmaCell[own[facei]] & Sf[facei];
        }
    }

    // Processor faces: the structure on this side, the soil across
    labelList maskLabels(cellMask_.size());
    forAll(maskLabels, celli)
    {
        maskLabels[celli] = (cellMask_[celli] > 0.5 ? 1 : 0);
    }
    labelList nbrMask;
    syncTools::swapBoundaryCellList(mesh_, maskLabels, nbrMask);

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fp = mesh_.boundary()[patchi];

        if (!fp.coupled())
        {
            continue;
        }

        const labelUList& fc = fp.faceCells();
        const vectorField& pSf = fp.Sf();
        const symmTensorField sigmaNbr
        (
            sigma0_.boundaryField()[patchi].patchNeighbourField()
          + sigma_.boundaryField()[patchi].patchNeighbourField()
        );
        const label start = fp.start() - mesh_.nInternalFaces();

        forAll(fc, i)
        {
            if (maskLabels[fc[i]] == 1 && nbrMask[start + i] == 0)
            {
                F += sigmaNbr[i] & pSf[i];
            }
        }
    }

    reduce(F, sumOp<vector>());

    return F;
}


void Foam::biot::structureLoads
(
    vector& Feff,
    vector& Fp,
    scalar& V
) const
{
    vector Meff;
    vector Mp;
    point xV;
    structureLoads(Feff, Fp, V, Meff, Mp, xV, Zero);
}


void Foam::biot::structureLoads
(
    vector& Feff,
    vector& Fp,
    scalar& V,
    vector& Meff,
    vector& Mp,
    point& xV,
    const point& xRef
) const
{
    Feff = Zero;
    Fp = Zero;
    V = 0;
    Meff = Zero;
    Mp = Zero;
    xV = Zero;

    // All ranks take part in the processor exchange and the reduction
    if (cellMask_.size() == mesh_.nCells())
    {
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const vectorField& Sf = mesh_.faceAreas();
        const scalarField& w = mesh_.weights().primitiveField();
        const symmTensorField& s0 = sigma0_.primitiveField();
        const symmTensorField& s = sigma_.primitiveField();
        const scalarField& pc = p_.primitiveField();
        const scalarField& Vc = mesh_.V();
        const pointField& C = mesh_.cellCentres();
        const vectorField& Cf = mesh_.faceCentres();

        forAll(cellMask_, celli)
        {
            if (cellMask_[celli] > 0.5)
            {
                V += Vc[celli];
                xV += Vc[celli]*C[celli];
            }
        }

        forAll(nei, facei)
        {
            const bool mo = cellMask_[own[facei]] > 0.5;
            const bool mn = cellMask_[nei[facei]] > 0.5;

            if (mo == mn)
            {
                continue;
            }

            // Face area vector out of the structure and the soil cell
            const vector S(mo ? Sf[facei] : -Sf[facei]);
            const label soilCell = mo ? nei[facei] : own[facei];

            // The geostatic stress is linear in depth: interpolated to the
            // face it closes the buoyancy of the stepped surface exactly
            const symmTensor s0f
            (
                w[facei]*s0[own[facei]] + (1 - w[facei])*s0[nei[facei]]
            );

            vector fEff((s0f + s[soilCell]) & S);
            if
            (
                faceContactIndex_.size()
             && faceContactIndex_[facei] >= 0
             && contactStates_[faceContactIndex_[facei]] == 1
            )
            {
                fEff = Zero;    // open face: no effective traction
            }
            const vector fP(-pc[soilCell]*S);
            Feff += fEff;
            Fp += fP;
            Meff += (Cf[facei] - xRef) ^ fEff;
            Mp += (Cf[facei] - xRef) ^ fP;
        }

        labelList maskLabels(cellMask_.size());
        forAll(maskLabels, celli)
        {
            maskLabels[celli] = (cellMask_[celli] > 0.5 ? 1 : 0);
        }
        labelList nbrMask;
        syncTools::swapBoundaryCellList(mesh_, maskLabels, nbrMask);

        // Contact states across the processor patches: one global swap
        // before the patch loop (a swap per coupled patch deadlocked the
        // ranks with different numbers of processor patches, 2026-10-28)
        const labelField procStates
        (
            contactFrictional()
          ? procContactStates()()
          : labelField(mesh_.nBoundaryFaces(), 0)
        );

        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fp = mesh_.boundary()[patchi];

            if (!fp.coupled())
            {
                continue;
            }

            const labelUList& fc = fp.faceCells();
            const vectorField& pSf = fp.Sf();
            const scalarField& pw = mesh_.weights().boundaryField()[patchi];
            const vectorField& pCf = fp.Cf();
            const symmTensorField s0Nbr
            (
                sigma0_.boundaryField()[patchi].patchNeighbourField()
            );
            const symmTensorField sNbr
            (
                sigma_.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField pNbr
            (
                p_.boundaryField()[patchi].patchNeighbourField()
            );
            const label start = fp.start() - mesh_.nInternalFaces();

            // Structure on this side, soil across (the other side counts
            // the opposite case)
            forAll(fc, i)
            {
                if (maskLabels[fc[i]] == 1 && nbrMask[start + i] == 0)
                {
                    const symmTensor s0f
                    (
                        pw[i]*s0[fc[i]] + (1 - pw[i])*s0Nbr[i]
                    );
                    vector fEff((s0f + sNbr[i]) & pSf[i]);
                    if (procStates.size() && procStates[start + i] == 1)
                    {
                        fEff = Zero;    // open face: no effective traction
                    }
                    const vector fP(-pNbr[i]*pSf[i]);
                    Feff += fEff;
                    Fp += fP;
                    Meff += (pCf[i] - xRef) ^ fEff;
                    Mp += (pCf[i] - xRef) ^ fP;
                }
            }
        }
    }

    reduce(Feff, sumOp<vector>());
    reduce(Fp, sumOp<vector>());
    reduce(V, sumOp<scalar>());
    reduce(Meff, sumOp<vector>());
    reduce(Mp, sumOp<vector>());
    reduce(xV, sumOp<vector>());
    if (V > VSMALL)
    {
        xV /= V;
    }
}


void Foam::biot::ensureStructureRefPoint() const
{
    if (structureRefSet_ || !structure_)
    {
        return;
    }
    vector Feff;
    vector Fp;
    scalar V;
    vector Meff;
    vector Mp;
    point xV;
    structureLoads(Feff, Fp, V, Meff, Mp, xV, Zero);
    structure_->setRefPoint(xV);
    structureRefSet_ = true;
}


Foam::vector Foam::biot::structureResultant() const
{
    ensureStructureRefPoint();
    vector Feff;
    vector Fp;
    scalar V;
    vector Meff;
    vector Mp;
    point xV;
    const point xRef(structure_ ? structure_->refPoint() : point(Zero));
    structureLoads(Feff, Fp, V, Meff, Mp, xV, xRef);

    // p is the excess over hydrostatic: the hydrostatic part is the
    // buoyancy rho_w g V; the weight and external force of a free
    // structure
    vector Fg(0, 0, rhoW_*g_*V);
    vector Fext(Zero);
    if (structureFree())
    {
        vector Mg;
        structure_->gravityLoads(V, xV, rhoW_, g_, Fg, Mg);
        Fext = structure_->externalForce(mesh_.time().value());
    }
    return Feff + Fp + Fg + Fext;
}


Foam::vector Foam::biot::structureResultantMoment() const
{
    ensureStructureRefPoint();
    vector Feff;
    vector Fp;
    scalar V;
    vector Meff;
    vector Mp;
    point xV;
    const point xRef(structure_ ? structure_->refPoint() : point(Zero));
    structureLoads(Feff, Fp, V, Meff, Mp, xV, xRef);

    vector Mg((xV - xRef) ^ vector(0, 0, rhoW_*g_*V));
    vector Mext(Zero);
    if (structureFree())
    {
        vector Fg;
        structure_->gravityLoads(V, xV, rhoW_, g_, Fg, Mg);
        Mext = structure_->externalMoment(mesh_.time().value());
    }
    return Meff + Mp + Mg + Mext;
}


void Foam::biot::updateStress()
{
    gradD_ = fvc::grad(D_);
    divD_ = fvc::div(D_);

    if (contactFrictional())
    {
        correctOpenFaceStrain();
    }

    if (volumetricSmoothing_ > 0)
    {
        // B-bar: replace the volumetric part of the cell strain by the
        // neighbour-averaged divergence, so that the volumetric
        // constraint of an isochoric (plastic) flow is not imposed cell by
        // cell (locking) and a cell squeezed by an imposed face
        // displacement shares its volume change with its neighbours
        volScalarField divDBar("divDBar", divD_);
        for (label i = 0; i < volumetricSmoothing_; ++i)
        {
            divDBar = fvc::average(fvc::interpolate(divDBar));
        }
        volTensorField gradDBar
        (
            "gradDBar",
            gradD_ + (divDBar - divD_)/3*tensor::I
        );
        law_->correct(gradDBar, sigma0_, sigma_);
    }
    else
    {
        law_->correct(gradD_, sigma0_, sigma_);
    }

    const volScalarField impK("impKExp", impKScale_*law_->impK());

    if (stabilisationScale_ < 1)
    {
        // The explicit counterpart of the implicit compact Laplacian is
        // taken partly with the compact scheme itself, which cancels the
        // stabilisation at convergence for that fraction
        divSigmaExp_ =
            fvc::div(sigma_, "div(sigmaD)")
          - stabilisationScale_*fvc::div(impK*gradD_, "div(sigmaD)")
          - (1 - stabilisationScale_)
           *fvc::laplacian(impK, D_, "laplacian(DD,D)");
    }
    else
    {
        divSigmaExp_ = fvc::div(sigma_ - impK*gradD_, "div(sigmaD)");
    }
}


void Foam::biot::correctOpenFaceStrain()
{
    // Gauss grad/div took the face value w D_own + (1 - w) D_nei across an
    // open face, which includes the jump across the gap: the soil cell
    // then saw a spurious strain of the size jump/cell, and with a low
    // conductivity its pore pressure ran away (pile case, 2026-10-28). On
    // an open face each side uses its own displacement (zero gradient:
    // the face is a free surface of that cell)
    if (contactFaces_.empty() && !contactProcAny_)
    {
        return;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const scalarField& w = mesh_.weights().primitiveField();
    const vectorField& Sf = mesh_.Sf().primitiveField();
    const scalarField& V = mesh_.V();
    const vectorField& D = D_.primitiveField();
    tensorField& gD = gradD_.primitiveFieldRef();
    scalarField& dD = divD_.primitiveFieldRef();

    forAll(contactFaces_, i)
    {
        if (contactStates_[i] != 1)
        {
            continue;
        }
        const label facei = contactFaces_[i];
        const label o = own[facei];
        const label n = nei[facei];
        const vector Df(w[facei]*D[o] + (1 - w[facei])*D[n]);
        // Sf points out of the owner
        gD[o] += Sf[facei]*(D[o] - Df)/V[o];
        dD[o] += (Sf[facei] & (D[o] - Df))/V[o];
        gD[n] -= Sf[facei]*(D[n] - Df)/V[n];
        dD[n] -= (Sf[facei] & (D[n] - Df))/V[n];
    }

    if (contactProcAny_)
    {
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fp = mesh_.boundary()[patchi];
            if (!fp.coupled())
            {
                continue;
            }
            // Collective: every rank evaluates every coupled patch
            const vectorField DNbr
            (
                D_.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField& pw = mesh_.weights().boundaryField()[patchi];
            const vectorField& pSf = fp.Sf();
            forAll(contactProcFaces_, i)
            {
                if (contactProcPatches_[i] != patchi || contactProcStates_[i] != 1)
                {
                    continue;
                }
                const label fi = contactProcFaces_[i];
                const label P = contactProcSoilCells_[i];
                const vector Df(pw[fi]*D[P] + (1 - pw[fi])*DNbr[fi]);
                gD[P] += pSf[fi]*(D[P] - Df)/V[P];
                dD[P] += (pSf[fi] & (D[P] - Df))/V[P];
            }
        }
    }

    gradD_.correctBoundaryConditions();
    divD_.correctBoundaryConditions();
}


Foam::tmp<Foam::scalarField> Foam::biot::maskedPressureImage() const
{
    // Mean pressure of the unmasked face neighbours of each masked cell
    // (processor neighbours included); masked cells without one take the
    // mean of their masked neighbours from the first sweep, and the rest
    // keep their value
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const scalarField& p = p_.primitiveField();

    scalarField image(p);
    bitSet isMasked(mesh_.nCells());
    for (const label celli : maskedCells_)
    {
        isMasked.set(celli);
    }

    labelList maskLabels(mesh_.nCells(), 0);
    for (const label celli : maskedCells_)
    {
        maskLabels[celli] = 1;
    }
    labelList nbrMask;
    syncTools::swapBoundaryCellList(mesh_, maskLabels, nbrMask);
    scalarField nbrP;
    syncTools::swapBoundaryCellList(mesh_, p, nbrP);

    for (label sweep = 0; sweep < 2; ++sweep)
    {
        scalarField sum(mesh_.nCells(), Zero);
        scalarField cnt(mesh_.nCells(), Zero);
        // Sweep 0: unmasked neighbours only; sweep 1: masked neighbours
        // that found a value in sweep 0 (image of the previous sweep)
        const scalarField& src = (sweep == 0 ? p : image);
        forAll(nei, facei)
        {
            const label o = own[facei];
            const label n = nei[facei];
            const bool mo = isMasked.test(o);
            const bool mn = isMasked.test(n);
            if (mo && (sweep == 0 ? !mn : mn))
            {
                sum[o] += src[n];
                cnt[o] += 1;
            }
            if (mn && (sweep == 0 ? !mo : mo))
            {
                sum[n] += src[o];
                cnt[n] += 1;
            }
        }
        if (sweep == 0)
        {
            forAll(mesh_.boundaryMesh(), patchi)
            {
                const polyPatch& pp = mesh_.boundaryMesh()[patchi];
                if (!pp.coupled())
                {
                    continue;
                }
                const labelUList& fc = pp.faceCells();
                forAll(fc, i)
                {
                    const label bfacei = pp.start() + i - mesh_.nInternalFaces();
                    if (isMasked.test(fc[i]) && nbrMask[bfacei] == 0)
                    {
                        sum[fc[i]] += nbrP[bfacei];
                        cnt[fc[i]] += 1;
                    }
                }
            }
        }
        scalarField next(image);
        for (const label celli : maskedCells_)
        {
            if (cnt[celli] > 0)
            {
                next[celli] = sum[celli]/cnt[celli];
            }
        }
        image = next;
    }

    auto tv = tmp<scalarField>::New(maskedCells_.size());
    scalarField& v = tv.ref();
    forAll(maskedCells_, i)
    {
        v[i] = image[maskedCells_[i]];
    }
    return tv;
}


void Foam::biot::setStructureMask(const scalarField& mask)
{
    soilModel::setStructureMask(mask);
    if (contactFrictional())
    {
        buildContactFaces();
    }
}


void Foam::biot::buildContactFaces()
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const vectorField& Sf = mesh_.faceAreas();
    const scalarField& magSf = mesh_.magSf().primitiveField();

    DynamicList<label> faces;
    DynamicList<label> soilCells;
    DynamicList<label> structCells;
    DynamicList<vector> normals;
    faceContactIndex_ = labelList(mesh_.nInternalFaces(), -1);

    if (cellMask_.size() == mesh_.nCells())
    {
        forAll(nei, facei)
        {
            const bool mo = cellMask_[own[facei]] > 0.5;
            const bool mn = cellMask_[nei[facei]] > 0.5;
            if (mo == mn)
            {
                continue;
            }
            faceContactIndex_[facei] = faces.size();
            faces.append(facei);
            soilCells.append(mo ? nei[facei] : own[facei]);
            structCells.append(mo ? own[facei] : nei[facei]);
            normals.append((mo ? Sf[facei] : -Sf[facei])/magSf[facei]);
        }
    }
    contactFaces_ = labelList(std::move(faces));
    contactSoilCells_ = labelList(std::move(soilCells));
    contactStructureCells_ = labelList(std::move(structCells));
    contactNormals_ = vectorField(std::move(normals));
    contactStates_ = labelList(contactFaces_.size(), 0);
    contactSlip_ = vectorField(contactFaces_.size(), Zero);
    contactSlip0_ = contactSlip_;

    // Processor faces with the soil on this side and the structure across
    DynamicList<label> pPatches;
    DynamicList<label> pFaces;
    DynamicList<label> pCells;
    DynamicList<vector> pNormals;
    if (cellMask_.size() == mesh_.nCells())
    {
        labelList maskLabels(cellMask_.size());
        forAll(maskLabels, celli)
        {
            maskLabels[celli] = (cellMask_[celli] > 0.5 ? 1 : 0);
        }
        labelList nbrMask;
        syncTools::swapBoundaryCellList(mesh_, maskLabels, nbrMask);
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fp = mesh_.boundary()[patchi];
            if (!fp.coupled())
            {
                continue;
            }
            const labelUList& fc = fp.faceCells();
            const vectorField& pSf = fp.Sf();
            const label start = fp.start() - mesh_.nInternalFaces();
            forAll(fc, i)
            {
                if (maskLabels[fc[i]] == 0 && nbrMask[start + i] == 1)
                {
                    pPatches.append(patchi);
                    pFaces.append(i);
                    pCells.append(fc[i]);
                    // Sf points out of the domain, i.e. into the structure
                    pNormals.append(-pSf[i]/mag(pSf[i]));
                }
            }
        }
    }
    contactProcPatches_ = labelList(std::move(pPatches));
    contactProcFaces_ = labelList(std::move(pFaces));
    contactProcSoilCells_ = labelList(std::move(pCells));
    contactProcNormals_ = vectorField(std::move(pNormals));
    contactProcStates_ = labelList(contactProcFaces_.size(), 0);
    contactProcSlip_ = vectorField(contactProcFaces_.size(), Zero);
    contactProcSlip0_ = contactProcSlip_;
    contactProcAny_ =
        returnReduce(contactProcFaces_.size(), sumOp<label>()) > 0;

    Info<< "    Structure contact: "
        << returnReduce(contactFaces_.size(), sumOp<label>())
        << " internal faces, "
        << returnReduce(contactProcFaces_.size(), sumOp<label>())
        << " processor faces" << endl;
}


Foam::tmp<Foam::labelField> Foam::biot::procContactStates() const
{
    // States of the soil-side faces to the structure side across the
    // processor patches; boundary-face numbering
    labelList states(mesh_.nBoundaryFaces(), 0);
    forAll(contactProcFaces_, i)
    {
        const fvPatch& fp = mesh_.boundary()[contactProcPatches_[i]];
        states[fp.start() - mesh_.nInternalFaces() + contactProcFaces_[i]] =
            contactProcStates_[i];
    }
    labelList nbr(states);
    syncTools::swapBoundaryFaceList(mesh_, nbr);
    return tmp<labelField>::New(nbr);
}


void Foam::biot::contactCounts(label& nOpen, label& nSlip) const
{
    nOpen = 0;
    nSlip = 0;
    for (const label s : contactStates_)
    {
        nOpen += (s == 1);
        nSlip += (s == 2);
    }
    for (const label s : contactProcStates_)
    {
        nOpen += (s == 1);
        nSlip += (s == 2);
    }
    reduce(nOpen, sumOp<label>());
    reduce(nSlip, sumOp<label>());
}


void Foam::biot::updateContact
(
    const volScalarField& impK,
    const bool updateStates
)
{
    // No early return for a rank without contact faces: the boundary
    // evaluation of the correction field and the neighbour exchange of
    // the processor-face branch are collective (a rank that skipped them
    // deadlocked the 16-rank pile case, 2026-10-28)
    contactCorr_.primitiveFieldRef() = Zero;
    label nChanged = 0;

    const labelUList& own = mesh_.owner();
    const scalarField& magSf = mesh_.magSf().primitiveField();
    const scalarField& w = mesh_.weights().primitiveField();
    const scalarField& dc = mesh_.nonOrthDeltaCoeffs().primitiveField();
    const scalarField& V = mesh_.V();
    const vectorField& D = D_.primitiveField();
    const tensorField& gD = gradD_.primitiveField();
    const symmTensorField& s0 = sigma0_.primitiveField();
    const symmTensorField& sg = sigma_.primitiveField();
    const scalarField& k = impK.primitiveField();
    const scalar s = stabilisationScale_;
    vectorField& corr = contactCorr_.primitiveFieldRef();

    // The force on the soil cell P across a contact face as the solver
    // represents it: implicit a_f (D_S - D_P) (the Laplacian) plus the
    // explicit flux of divSigmaExp (div(sigma) minus the explicit part of
    // the Laplacian, in the stabilised form); Sf_out points out of the
    // structure into the soil, so the outward normal of P at the face is
    // -n
    auto contactFace = [&]
    (
        const label P,
        const vector& n,
        const scalar A,
        const scalar af,
        const vector& DS,
        const symmTensor& sigF,
        const tensor& gradDF,
        const scalar kF,
        label& state,
        vector& slip,
        vector& slip0,
        vector& corrP
    )
    {
        const vector SfOut(n*A);
        const vector dRel(D[P] - DS);
        const vector implicitF(af*(DS - D[P]));
        const vector explicitF
        (
          - (sigF & SfOut)
          + s*((kF*gradDF) & SfOut)
          - (1 - s)*implicitF
        );
        // Effective traction of the soil cell on the face (tension > 0)
        const vector t((s0[P] + sg[P]) & n);
        const scalar tn = t & n;

        // Tangential force on the soil of a closed face [N]: the force the
        // bonded coupling transmits (implicit plus explicit flux, the
        // continuum traction at convergence) minus the accumulated slip
        // times the coupling coefficient. Coulomb: when it exceeds the
        // limit, the excess becomes slip (return mapping along the force)
        // and the transmitted force is at the limit, in the direction the
        // structure drags the soil. The slip is frozen between the state
        // updates, so the iterations in between see the bonded (linear)
        // coupling with a constant offset and converge; a force that
        // follows the iterate or a coupling cut when slipping cycles (a
        // cut face never re-sticks: its relative displacement grows freely).
        // The mapping starts from the converged slip of the previous step
        // (slip0) at every update: the slip of an overshooting iterate is
        // undone once the iterate corrects, so the converged state does
        // not depend on the iteration path
        const tensor T(tensor::I - sqr(n));
        const scalar limitF =
            (contactFriction_*max(-tn, scalar(0)) + contactAdhesion_)*A;
        const vector bondedT(T & (implicitF + explicitF));

        if (updateStates)
        {
            const label state0 = state;
            if (state != 1 && tn > contactTension_)
            {
                state = 1;
            }
            else if (state == 1 && (dRel & n) < 0)
            {
                // Closing: no tangential force stored across the opening
                state = 0;
                slip = bondedT/max(af, SMALL);
                slip0 = slip;
            }
            if (state != 1)
            {
                const vector trial(bondedT - af*(T & slip0));
                const scalar magT = mag(trial);
                vector target(slip0);
                if (magT > limitF && magT > SMALL)
                {
                    state = 2;
                    target += (magT - limitF)/max(af, SMALL)*(trial/magT);
                }
                else
                {
                    state = 0;
                }
                slip += contactFrictionRelax_*(target - slip);
            }
            nChanged += (state != state0);
        }

        if (state == 1)
        {
            // Traction-free: cancel the explicit flux (the implicit
            // coupling is removed from the matrix)
            corrP = -explicitF;
        }
        else
        {
            // Bonded coupling (implicit and explicit) minus the slip
            corrP = -af*(T & slip);
        }
    };

    forAll(contactFaces_, i)
    {
        const label facei = contactFaces_[i];
        const label P = contactSoilCells_[i];
        const label S = contactStructureCells_[i];
        const scalar wf = w[facei];
        const label o = own[facei];
        // Linear interpolation to the face (owner weight wf)
        const scalar wP = (o == P ? wf : 1 - wf);
        const symmTensor sigF
        (
            wP*(s0[P] + sg[P]) + (1 - wP)*(s0[S] + sg[S])
        );
        const tensor gradDF(wP*gD[P] + (1 - wP)*gD[S]);
        const scalar kF = wP*k[P] + (1 - wP)*k[S];
        const scalar af = kF*magSf[facei]*dc[facei];
        vector corrP;
        contactFace
        (
            P, contactNormals_[i], magSf[facei], af, D[S],
            sigF, gradDF, kF, contactStates_[i], contactSlip_[i],
            contactSlip0_[i], corrP
        );
        corr[P] += corrP/V[P];
    }

    // Processor faces: the structure values from across the patch
    if (contactProcAny_)
    {
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fp = mesh_.boundary()[patchi];
            if (!fp.coupled())
            {
                continue;
            }
            const vectorField DNbr(D_.boundaryField()[patchi].patchNeighbourField());
            const tensorField gNbr(gradD_.boundaryField()[patchi].patchNeighbourField());
            const symmTensorField s0Nbr(sigma0_.boundaryField()[patchi].patchNeighbourField());
            const symmTensorField sNbr(sigma_.boundaryField()[patchi].patchNeighbourField());
            const scalarField kNbr(impK.boundaryField()[patchi].patchNeighbourField());
            const scalarField& pw = mesh_.weights().boundaryField()[patchi];
            const scalarField& pdc = mesh_.nonOrthDeltaCoeffs().boundaryField()[patchi];
            const scalarField& pMagSf = fp.magSf();
            forAll(contactProcFaces_, i)
            {
                if (contactProcPatches_[i] != patchi)
                {
                    continue;
                }
                const label fi = contactProcFaces_[i];
                const label P = contactProcSoilCells_[i];
                const symmTensor sigF
                (
                    pw[fi]*(s0[P] + sg[P]) + (1 - pw[fi])*(s0Nbr[fi] + sNbr[fi])
                );
                const tensor gradDF(pw[fi]*gD[P] + (1 - pw[fi])*gNbr[fi]);
                const scalar kF = pw[fi]*k[P] + (1 - pw[fi])*kNbr[fi];
                const scalar af = kF*pMagSf[fi]*pdc[fi];
                vector corrP;
                contactFace
                (
                    P, contactProcNormals_[i], pMagSf[fi], af, DNbr[fi],
                    sigF, gradDF, kF, contactProcStates_[i],
                    contactProcSlip_[i], contactProcSlip0_[i], corrP
                );
                corr[P] += corrP/V[P];
            }
        }
    }
    contactCorr_.correctBoundaryConditions();

    if (updateStates)
    {
        reduce(nChanged, sumOp<label>());
        contactChanged_ = nChanged;
    }
}


void Foam::biot::applyContact(fvVectorMatrix& DEqn) const
{
    // Open faces: no implicit coupling between the soil and the structure
    // cell (the Laplacian coefficient leaves both diagonals)
    scalarField& upper = DEqn.upper();
    scalarField& diag = DEqn.diag();
    const bool hasLower = DEqn.hasLower();
    forAll(contactFaces_, i)
    {
        if (contactStates_[i] == 1)
        {
            const label facei = contactFaces_[i];
            const scalar a = upper[facei];
            upper[facei] = 0;
            if (hasLower)
            {
                DEqn.lower()[facei] = 0;
            }
            diag[contactSoilCells_[i]] += a;
            diag[contactStructureCells_[i]] += a;
        }
    }
    forAll(contactProcFaces_, i)
    {
        if (contactProcStates_[i] == 1)
        {
            const label patchi = contactProcPatches_[i];
            const label fi = contactProcFaces_[i];
            DEqn.internalCoeffs()[patchi][fi] = Zero;
            DEqn.boundaryCoeffs()[patchi][fi] = Zero;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::biot::biot(const fvMesh& mesh, const dictionary& dict)
:
    soilModel(mesh, dict),
    porosity_(dict.getOrDefault<scalar>("porosity", 0.4)),
    conductivity_(dict.get<scalar>("conductivity")),
    bulkModulus_(readBulkModulus(dict)),
    storage_
    (
        "storage",
        dimensionSet(-1, 1, 2, 0, 0),
        porosity_/bulkModulus_
    ),
    mobility_
    (
        "mobility",
        dimensionSet(-1, 3, 1, 0, 0),
        conductivity_/gammaW()
    ),
    nCorrectors_(dict.getOrDefault<label>("nCorrectors", 30)),
    tolerance_(dict.getOrDefault<scalar>("tolerance", 1e-6)),
    aitken_(dict.getOrDefault<bool>("aitken", false)),
    aitkenPressure_(dict.getOrDefault<bool>("aitkenPressure", false)),
    solverType_(dict.getOrDefault<word>("solver", "segregated")),
    krylovDim_(dict.getOrDefault<label>("krylovDim", 20)),
    krylovRestarts_(dict.getOrDefault<label>("krylovRestarts", 2)),
    maxNewton_(dict.getOrDefault<label>("maxNewton", 10)),
    newtonTol_(dict.getOrDefault<scalar>("newtonTolerance", 1e-6)),
    fdEpsilon_(dict.getOrDefault<scalar>("fdEpsilon", 1e-7)),
    blockInnerRelTol_(dict.getOrDefault<scalar>("blockInnerRelTol", 1e-3)),
    volumetricSmoothing_
    (
        dict.getOrDefault<label>("volumetricSmoothing", 0)
    ),
    impKScale_(dict.getOrDefault<scalar>("impKScale", 1)),
    stabilisationScale_
    (
        dict.getOrDefault<scalar>("stabilisationScale", 1)
    ),
    bottomPressure_
    (
        Function1<scalar>::NewIfPresent
        (
            "bottomPressure", dict, word::null, &mesh
        )
    ),
    sealedTop_(dict.getOrDefault<bool>("sealedTop", false)),
    ownStructure_
    (
        dict.found("structureMotion")
      ? new rigidStructure
        (
            "structure",
            dict.subDict("structureMotion"),
            vector(0, 0, 1)
        )
      : nullptr
    ),
    structureRefSet_(false),
    restored_(false),
    restoredStressPending_(false),
    pTop_(),
    law_(soilConstitutiveModel::New(mesh, dict)),
    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPressure, dict.getOrDefault<scalar>("p0", 0)),
        pTypes()
    ),
    D_
    (
        IOobject
        (
            "D",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimLength, Zero),
        DTypes(dict)
    ),
    gradD_
    (
        IOobject("gradD", mesh.time().timeName(), mesh),
        mesh,
        dimensionedTensor(dimless, Zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigmaEff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor(dimPressure, Zero)
    ),
    sigma0_
    (
        IOobject("sigma0", mesh.time().timeName(), mesh),
        mesh,
        dimensionedSymmTensor(dimPressure, Zero),
        zeroGradientFvPatchField<symmTensor>::typeName
    ),
    divD_
    (
        IOobject("divD", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(dimless, Zero)
    ),
    divSigmaExp_
    (
        IOobject("divSigmaExp", mesh.time().timeName(), mesh),
        mesh,
        dimensionedVector(dimPressure/dimLength, Zero)
    ),
    contactType_
    (
        dict.subOrEmptyDict("structureContact")
           .getOrDefault<word>("type", "bonded")
    ),
    contactFriction_
    (
        dict.subOrEmptyDict("structureContact")
           .getOrDefault<scalar>("friction", 0.5)
    ),
    contactAdhesion_
    (
        dict.subOrEmptyDict("structureContact")
           .getOrDefault<scalar>("adhesion", 0)
    ),
    contactTension_
    (
        dict.subOrEmptyDict("structureContact")
           .getOrDefault<scalar>("tensionTolerance", 0)
    ),
    contactGate_
    (
        dict.subOrEmptyDict("structureContact")
           .getOrDefault<scalar>("gate", 1e-3)
    ),
    contactMaxUpdates_
    (
        dict.subOrEmptyDict("structureContact")
           .getOrDefault<label>("maxUpdates", 50)
    ),
    contactFrictionRelax_
    (
        dict.subOrEmptyDict("structureContact")
           .getOrDefault<scalar>("slipRelaxation", 1)
    ),
    contactChanged_(0),
    contactProcAny_(false),
    contactFaces_(),
    contactSoilCells_(),
    contactStructureCells_(),
    contactNormals_(),
    contactStates_(),
    faceContactIndex_(),
    contactSlip_(),
    contactSlip0_(),
    contactProcPatches_(),
    contactProcFaces_(),
    contactProcSoilCells_(),
    contactProcNormals_(),
    contactProcStates_(),
    contactProcSlip_(),
    contactProcSlip0_(),
    contactCorr_
    (
        IOobject("contactCorr", mesh.time().timeName(), mesh),
        mesh,
        dimensionedVector(dimPressure/dimLength, Zero)
    ),
    topTraction_(),
    topFraction_(),
    topDisplacement_(),
    topTangential_(false)
{
    // The mixed top condition is constructed with unset reference fields:
    // traction-free until the first update
    {
        auto& Dp =
            refCast<directionMixedFvPatchVectorField>
            (
                D_.boundaryFieldRef()[topPatchi_]
            );
        Dp.refValue() = Zero;
        Dp.refGrad() = Zero;
        Dp.valueFraction() = Zero;
        Dp.evaluate();
    }

    // Keep the old-time volumetric strain for the storage term
    divD_.oldTime();

    Info<< "    Biot soil: porosity " << porosity_
        << ", conductivity " << conductivity_ << " m/s, bulk modulus "
        << bulkModulus_ << " Pa, diffusivity " << diffusivity() << " m2/s, "
        << nCorrectors_ << " correctors, tolerance " << tolerance_
        << (aitken_ ? ", Aitken" : "")
        << (volumetricSmoothing_ > 0 ? ", B-bar" : "")
        << (stabilisationScale_ < 1 ? ", stabilisation scaled" : "")
        << (bottomPressure_ ? ", fixed bottom pressure" : "") << endl;

    if (contactType_ != "bonded" && contactType_ != "frictional")
    {
        FatalIOErrorInFunction(dict)
            << "structureContact type " << contactType_
            << ": bonded or frictional" << exit(FatalIOError);
    }
    if (contactFrictional())
    {
        if (solverType_ != "segregated")
        {
            FatalIOErrorInFunction(dict)
                << "the frictional structure contact needs the segregated"
                << " solver" << exit(FatalIOError);
        }
        Info<< "    Structure contact: frictional, friction "
            << contactFriction_ << ", adhesion " << contactAdhesion_
            << " Pa, tension tolerance " << contactTension_ << " Pa" << endl;
    }

    // Restart: the displacement and pore pressure written at the
    // restart time are read when present (the region positions its
    // columns on the written bed before the model is constructed, so
    // that the fields are read on the mesh they were written on)
    restored_ =
        IOobject
        (
            "D",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ).typeHeaderOk<volVectorField>(true);

    if (restored_)
    {
        D_.oldTime() = D_;
        p_.oldTime() = p_;
        gradD_ = fvc::grad(D_);
        divD_ = fvc::div(D_);
        divD_.oldTime() = divD_;
        restoredStressPending_ = true;
        Info<< "    Biot soil: displacement and pore pressure restored"
            << " from time " << mesh.time().timeName()
            << (contactFrictional() ? " (contact slip starts from stick)" : "")
            << endl;
    }

    if (ownStructure_)
    {
        attachStructure(ownStructure_.get());
        ownStructure_->readState(mesh.time());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::biot::diffusivity() const
{
    const scalar M = law_->constrainedModulus();
    return conductivity_/(gammaW()*(1/M + porosity_/bulkModulus_));
}


void Foam::biot::writeState() const
{
    if (ownStructure_)
    {
        ownStructure_->writeState(mesh_.time());
    }
}


void Foam::biot::initialise(const scalarField& pCells)
{
    p_.primitiveFieldRef() = pCells;
    p_.correctBoundaryConditions();
    p_.oldTime() = p_;
}


void Foam::biot::remap
(
    const labelList& c0,
    const labelList& c1,
    const scalarField& w
)
{
    initialise(remapField(p_.primitiveField(), c0, c1, w));

    D_.primitiveFieldRef() = remapField(D_.primitiveField(), c0, c1, w);
    D_.correctBoundaryConditions();
    D_.oldTime() = D_;

    law_->remap(c0, c1, w);
    remapGeneration(c0, c1, w);

    updateStress();
    divD_.oldTime() = divD_;
}


void Foam::biot::geostatic(const scalarField& sigmaV0, const scalar k0)
{
    symmTensorField& s0 = sigma0_.primitiveFieldRef();

    forAll(s0, celli)
    {
        // Compression negative (tension positive)
        s0[celli] = symmTensor
        (
            -k0*sigmaV0[celli], 0, 0,
            -k0*sigmaV0[celli], 0,
            -sigmaV0[celli]
        );
    }

    sigma0_.correctBoundaryConditions();

    law_->geostatic(sigma0_);
}


void Foam::biot::correct(const scalarField& pTop, const label nSteps)
{
    pTop_ = pTop;
    if (!sealedTop_)
    {
        p_.boundaryFieldRef()[topPatchi_] == pTop;
    }

    if (bottomPressure_)
    {
        p_.boundaryFieldRef()[bottomPatchi_] ==
            bottomPressure_->value(mesh_.time().value());
    }

    law_->beginStep();

    // Restart: the stress and its explicit divergence of the restored
    // displacement, which the first corrector takes from the previous
    // iterate (zero at construction; the geostatic stress is set by now)
    if (restoredStressPending_)
    {
        updateStress();
        restoredStressPending_ = false;
    }

    const volScalarField impK("impKEqn", impKScale_*law_->impK());

    // Stabilised segregated coupling (Tang et al. 2015): the skeleton
    // compressibility 1/M is added to the storage implicitly and removed
    // explicitly with the latest iterate, so that the fixed-point
    // iteration contracts for any pore-fluid stiffness (the plain
    // explicit coupling has a gain K'/(n M), > 1 for saturated soil)
    // Per cell, 1/K of the law's current drained bulk modulus (the
    // fixed-stress split is stable when the estimated compliance is not
    // below the true one: the bulk modulus, not the constrained one, and
    // a softening skeleton reports its softer modulus); refreshed every
    // corrector
    volScalarField invM
    (
        IOobject("invM", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(storage_.dimensions(), 1/law_->constrainedModulus())
    );
    const dimensionedScalar invMsrc
    (
        "invMsrc",
        storage_.dimensions(),
        1/law_->constrainedModulus()
    );
    auto updateInvM = [&]()
    {
        invM =
            dimensionedScalar(storage_.dimensions()*dimPressure, 1)
           /law_->stabilisationModulusField();
    };
    updateInvM();

    // Scale of the relative pressure change: the bed pressure or the pore
    // pressure, floored by 1 % of the geostatic stress so that a drained
    // run (p ~ 0) converges on the displacement
    const scalar pScale =
        max
        (
            max(gMax(mag(pTop)), gMax(mag(p_.primitiveField()))),
            max(0.01*gMax(mag(sigma0_.primitiveField())), SMALL)
        );

    // Generation rate of this interval (evaluated once: it advances the
    // generated pressure)
    volScalarField::Internal fGen
    (
        IOobject("poreGenerationRate", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimPressure/dimTime, Zero)
    );
    {
        const tmp<scalarField> tf(generationRate(nSteps, pTop));
        if (tf().size())
        {
            fGen.field() = tf();
        }
        else
        {
            fGen.field().clear();
        }
    }

    // Balance of the free degrees of freedom of the structure (force and
    // moment about its reference point; quasi-static or dynamic with the
    // structure's inertia) from the loads of the current soil iterate;
    // returns the structure residual
    auto updateStructure = [&](const label iCorr) -> scalar
    {
        vector Feff;
        vector Fp;
        scalar V;
        vector Meff;
        vector Mp;
        point xV;
        structureLoads
        (
            Feff, Fp, V, Meff, Mp, xV, structure_->refPoint()
        );
        const scalar t = mesh_.time().value();
        vector Fg;
        vector Mg;
        structure_->gravityLoads(V, xV, rhoW_, g_, Fg, Mg);
        const vector Fext(structure_->externalForce(t));
        const vector Mext(structure_->externalMoment(t));
        const vector F(Feff + Fp + Fg + Fext);
        const vector M(Meff + Mp + Mg + Mext);

        // Residual scales: the weight or the external load, and that
        // force over the structure size
        const scalar fScale =
            max(max(structure_->mass(V)*g_, mag(Fext)), SMALL);
        const scalar mScale =
            max(max(fScale*cbrt(max(V, VSMALL)), mag(Mext)), SMALL);

        const scalar changeS =
            structure_->correct(F, M, structure_->mass(V), fScale, mScale);
        setStructureMotion
        (
            structure_->displacement(),
            structure_->rotation(),
            structure_->refPoint()
        );
        if (debug)
        {
            Info<< "    structure corrector " << iCorr << ": d "
                << structure_->displacement() << ", theta "
                << structure_->rotation() << ", F " << F << ", M " << M
                << endl;
        }
        return changeS;
    };

    label iCorr = 0;
    // Pressure of the structure cells for this solve: the image of the
    // soil pressure next to them, taken once per solve so that it does
    // not move with the correctors (refreshed every corrector it kept the
    // change test from converging, lateralInclusion 2026-10-30).
    // Collective (processor exchange inside): every rank calls it
    const scalarField pImage(maskedPressureImage());

    scalar residualD = GREAT;
    scalar residualP = GREAT;
    scalar changeP = GREAT;
    scalar changeD = GREAT;
    scalar changeS = 0;

    label nContactUpdates = 0;

    // Free structure: not converged before its first update of the step
    const bool structureFree_ = structureFree();
    if (structureFree_)
    {
        if (solverType_ != "segregated")
        {
            FatalErrorInFunction
                << "a free structure needs the segregated solver"
                << exit(FatalError);
        }
        changeS = structure_->outerUpdate() ? 0 : GREAT;
        structure_->newStep
        (
            mesh_.time().deltaTValue(),
            mesh_.time().value()
        );
        ensureStructureRefPoint();

        // The structure cells follow the structure's state (its last
        // update, or the state restored at a restart: the mask was set
        // before the state was read, with the cells at rest)
        setStructureMotion
        (
            structure_->displacement(),
            structure_->rotation(),
            structure_->refPoint()
        );
    }

    bool fallback = false;
    if (solverType_ == "newtonKrylov")
    {
        const bool converged =
            solveNewtonKrylov
            (
                impK, invM, invMsrc, fGen, nSteps, pScale,
                iCorr, changeP, changeD
            );

        if (!converged)
        {
            // The segregated correctors finish the step from the last
            // Newton iterate (a Picard fixed point does not need the
            // smoothness that Newton does across the kinks of a law)
            Info<< "    Biot Newton-Krylov: falling back to the segregated"
                << " correctors" << endl;
            fallback = true;
            iCorr = 0;
        }
    }

    if (solverType_ == "blockCoupled")
    {
        const bool converged =
            solveBlockCoupled
            (
                impK, invM, invMsrc, fGen, nSteps, pScale,
                iCorr, changeP, changeD
            );

        if (!converged || debug)
        {
            Info<< "    Biot block-coupled: " << iCorr << " outer iterations,"
                << " relative changes p " << changeP << " D " << changeD
                << (converged ? "" : " (not converged)") << endl;
        }

        law_->updateState();

        if (generation_)
        {
            const symmTensorField sigmaCell
            (
                sigma0_.primitiveField() + sigma_.primitiveField()
            );
            generation_->updateEnvelope
            (
                sigmaCell.component(symmTensor::XZ),
                mesh_.time().value()
            );
        }

        return;
    }

    if (solverType_ == "newtonKrylov" && !fallback)
    {
        law_->updateState();

        if (generation_)
        {
            const symmTensorField sigmaCell
            (
                sigma0_.primitiveField() + sigma_.primitiveField()
            );
            generation_->updateEnvelope
            (
                sigmaCell.component(symmTensor::XZ),
                mesh_.time().value()
            );
        }

        return;
    }
    else if (solverType_ != "segregated" && !fallback)
    {
        FatalIOErrorInFunction(dict_)
            << "Unknown solver " << solverType_
            << ": segregated, newtonKrylov or blockCoupled"
            << exit(FatalIOError);
    }

    // Aitken (delta-squared) relaxation of the displacement iterates, and
    // of the pressure iterates (the fixed-stress split contracts slowly
    // when the skeleton is much softer than the pore fluid, near
    // liquefaction)
    vectorField rPrev;
    scalar omega = 1;
    scalarField rPrevP;
    scalar omegaP = 1;

    do
    {
        updateTraction();

        if (contactFrictional())
        {
            // The active set changes once the soil iterates have settled
            // for the current set, a bounded number of times per step
            const bool updateStates =
                iCorr == 0
             || (max(changeP, changeD) < contactGate_
              && nContactUpdates < contactMaxUpdates_);
            if (updateStates)
            {
                ++nContactUpdates;
            }
            updateContact(impK, updateStates);
        }

        fvVectorMatrix DEqn
        (
            fvm::laplacian(impK, D_, "laplacian(DD,D)")
          + divSigmaExp_
          + contactCorr_
         ==
            fvc::grad(p_)
        );

        if (contactFrictional())
        {
            applyContact(DEqn);
        }

        if (maskedCells_.size())
        {
            DEqn.setValues(maskedCells_, maskedDisplacement_);
        }

        const vectorField DPrev(D_.primitiveField());

        residualD = cmptMax(DEqn.solve().initialResidual());

        if (aitken_)
        {
            const vectorField r(D_.primitiveField() - DPrev);

            if (iCorr > 0)
            {
                const vectorField dr(r - rPrev);
                const scalar denom = gSum(dr & dr);
                if (denom > VSMALL)
                {
                    omega = -omega*gSum(rPrev & dr)/denom;
                    omega = min(max(omega, scalar(0.1)), scalar(10));
                    D_.primitiveFieldRef() = DPrev + omega*r;
                    D_.correctBoundaryConditions();
                }
            }

            rPrev = r;
        }

        changeD =
            gMax(mag(D_.primitiveField() - DPrev))
           /max(gMax(mag(D_.primitiveField())), SMALL);

        updateStress();
        updateInvM();

        const scalarField pPrev(p_.primitiveField());

        // Over nSteps flow time steps (dt-based derivatives, mobility
        // scaled by nSteps)
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
            // Residual generation: the undrained rate f is
            // defined for the laterally confined column under constant
            // total vertical stress, where the skeleton takes p/M
            // elastically, so the storage seen by f is n/K' + 1/M (the
            // constrained modulus, independent of the bulk-modulus
            // compliance of the stabilisation)
            pEqn -= nSteps*(storage_ + invMsrc)*fGen;
        }

        if (maskedCells_.size())
        {
            // Cells inside a structure are impermeable and decoupled (zero
            // mobility on their faces). Left free, their pressure drifted
            // to tens of kPa, which inflated the residual normalisation of
            // the linear solver so that the cells beside the structure were
            // never solved (a 2 kPa ring that did not drain, pile case
            // 2026-10-29). Pinned to zero, the soil cells beside the
            // structure saw a spurious pressure-gradient force towards it
            // (the lateral inclusion stiffened by 10 %). They are pinned to
            // the image of the soil pressure (zero normal gradient at the
            // impermeable face), which is also what the Gauss gradient of
            // the neighbouring soil cells needs.
            pEqn.setValues(maskedCells_, pImage);
        }

        residualP = pEqn.solve().initialResidual();

        if (aitkenPressure_)
        {
            const scalarField r(p_.primitiveField() - pPrev);

            if (iCorr > 0)
            {
                const scalarField dr(r - rPrevP);
                const scalar denom = gSum(dr*dr);
                if (denom > VSMALL)
                {
                    omegaP = -omegaP*gSum(rPrevP*dr)/denom;
                    omegaP = min(max(omegaP, scalar(0.1)), scalar(10));
                    p_.primitiveFieldRef() = pPrev + omegaP*r;
                    p_.correctBoundaryConditions();
                }
            }

            rPrevP = r;
        }

        changeP = gMax(mag(p_.primitiveField() - pPrev))/pScale;

        // With a frictional contact the structure moves once the active
        // set has settled (or stopped updating): the loads of a changing
        // set are not a function of the structure's state. (Updating the
        // structure jointly with the set made the block rocking test
        // diverge, 2026-10-29; the pile placement, where the set cannot
        // settle at the first-guess position, is an open item.)
        const bool contactSettled =
            !contactFrictional()
         || contactChanged_ == 0
         || nContactUpdates >= contactMaxUpdates_;

        if
        (
            structureFree_
         && !structure_->outerUpdate()
         && contactSettled
         && max(changeP, changeD) < structure_->gate()
        )
        {
            // Updated once the soil iterates have settled for the current
            // state (the loads of unsettled iterates are not a function
            // of the state and the joint fixed point diverges)
            changeS = updateStructure(iCorr);
        }

    } while
    (
        max(max(changeP, changeD), changeS) > tolerance_
     && ++iCorr < nCorrectors_
    );

    // Outer update: once per solve, from the converged reaction (a pure
    // function of the structure's state: a clean secant across the
    // solves; the balance converges over the solves)
    if (structureFree_ && structure_->outerUpdate())
    {
        changeS = updateStructure(iCorr);
    }

    if (contactFrictional())
    {
        // Converged slip: the reference of the return mapping of the
        // next step
        contactSlip0_ = contactSlip_;
        contactProcSlip0_ = contactProcSlip_;
    }

    law_->updateState();

    if (generation_)
    {
        // Envelope of the solved shear stress for the computed amplitude
        const symmTensorField sigmaCell
        (
            sigma0_.primitiveField() + sigma_.primitiveField()
        );
        generation_->updateEnvelope
        (
            sigmaCell.component(symmTensor::XZ),
            mesh_.time().value()
        );
    }

    if (debug || iCorr >= nCorrectors_ - 1)
    {
        Info<< "    Biot: " << iCorr + 1 << " correctors, residuals D "
            << residualD << " p " << residualP << ", relative changes p "
            << changeP << " D " << changeD;
        if (structureFree_)
        {
            Info<< ", structure residual " << changeS;
        }
        if (contactFrictional())
        {
            label nOpen, nSlip;
            contactCounts(nOpen, nSlip);
            Info<< ", contact faces open " << nOpen << " slip " << nSlip
                << " (" << nContactUpdates << " state updates)";
        }
        Info<< endl;
    }
    if (structureFree_ && debug)
    {
        structure_->report();
    }
}


// ************************************************************************* //
