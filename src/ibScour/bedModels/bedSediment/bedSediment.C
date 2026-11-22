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

#include "bedSediment.H"
#include "faCFD.H"
#include "indexedOctree.H"
#include "PatchFunction1.H"
#include "triSurfaceMesh.H"
#include "treeDataPoint.H"
#include "syncTools.H"
#include "surfaceWriter.H"
#include "IOdictionary.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bedSediment, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::bedSediment::makeFaMesh(const word& patchName)
{
    const polyBoundaryMesh& pbm = bedMesh_.boundaryMesh();

    if (pbm.findPatchID(patchName) < 0)
    {
        FatalIOErrorInFunction(dict_)
            << "Bed patch " << patchName << " not found"
            << exit(FatalIOError);
    }

    const polyPatch& bedPatch = pbm[patchName];

    // The area mesh written by makeFaMesh (constant/finite-area/faMesh),
    // if present: decomposePar and reconstructPar then handle the bed
    // fields like any other field
    aMeshPtr_ = faMesh::TryNew(bedMesh_);
    if (aMeshPtr_)
    {
        label nOff = 0;
        for (const label facei : aMeshPtr_->faceLabels())
        {
            const label patchFacei = facei - bedPatch.start();
            if (patchFacei < 0 || patchFacei >= bedPatch.size())
            {
                ++nOff;
            }
        }
        const label nArea =
            returnReduce(aMeshPtr_->nFaces(), sumOp<label>());
        const label nPatch = returnReduce(bedPatch.size(), sumOp<label>());
        reduce(nOff, sumOp<label>());

        if (nOff || nArea != nPatch)
        {
            FatalIOErrorInFunction(dict_)
                << "The finite-area mesh in constant/finite-area has "
                << nArea << " faces, " << nOff << " of them not on the bed "
                << "patch " << patchName << " (" << nPatch << " faces): "
                << "polyMeshPatches of system/finite-area/faMeshDefinition "
                << "must be the bed patch alone" << exit(FatalIOError);
        }

        Info<< "    Height-field bed on patch " << patchName << ": "
            << nArea << " bed faces (finite-area mesh read from constant)"
            << endl;
        return;
    }

    // Otherwise built here: one faPatch per neighbouring (non-coupled)
    // polyPatch; edges between processors become processor faPatches
    dictionary definition;
    definition.add("polyMeshPatches", wordList(1, patchName));

    dictionary boundary;
    for (const polyPatch& pp : pbm)
    {
        if (pp.coupled() || pp.name() == patchName)
        {
            continue;
        }

        dictionary patchDict;
        patchDict.add("type", "patch");
        patchDict.add("neighbourPolyPatch", pp.name());
        boundary.add(pp.name(), patchDict);
    }
    definition.add("boundary", boundary);

    dictionary defaultPatch;
    defaultPatch.add("name", "bedBoundary");
    defaultPatch.add("type", "patch");
    definition.add("defaultPatch", defaultPatch);

    aMeshPtr_.reset(new faMesh(bedMesh_, definition));

    Info<< "    Height-field bed on patch " << patchName << ": "
        << returnReduce(aMeshPtr_->nFaces(), sumOp<label>())
        << " bed faces" << endl;
}


bool Foam::bedSediment::isOpenBoundary(const faPatch& p) const
{
    // By the name of the faPatch or of its neighbouring polyPatch (a
    // faMeshDefinition may name the edge patches differently)
    if (openBoundaries_.found(p.name()))
    {
        return true;
    }
    const label ngb = p.ngbPolyPatchIndex();
    return
        ngb >= 0
     && openBoundaries_.found(bedMesh_.boundaryMesh()[ngb].name());
}


Foam::tmp<Foam::scalarField> Foam::bedSediment::bedLevelInput
(
    const word& patchName,
    const word& valueKey,
    const word& surfaceKey,
    const scalar defaultValue
) const
{
    const faMesh& aMesh = *aMeshPtr_;
    const labelList& faceLabels = aMesh.faceLabels();
    const polyPatch& pp = bedMesh_.boundaryMesh()[patchName];

    auto tz = tmp<scalarField>::New(faceLabels.size(), defaultValue);
    scalarField& z = tz.ref();

    if (dict_.found(surfaceKey))
    {
        const fileName surfName(dict_.get<fileName>(surfaceKey));
        const triSurfaceMesh surf
        (
            IOobject
            (
                surfName,
                bedMesh_.time().constant(),
                "triSurface",
                bedMesh_.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        // Vertical lines through the face centres, from above the surface
        // to below it: the first hit from the top is the bed
        const boundBox& bb = surf.bounds();
        const scalar margin = 1 + bb.span().z();
        pointField start(aMesh.areaCentres().primitiveField());
        pointField end(start);
        start.replace(vector::Z, bb.max().z() + margin);
        end.replace(vector::Z, bb.min().z() - margin);

        List<pointIndexHit> hits;
        surf.findLine(start, end, hits);

        label nMiss = 0;
        forAll(hits, i)
        {
            if (hits[i].hit())
            {
                z[i] = hits[i].point().z();
            }
            else
            {
                ++nMiss;
            }
        }
        reduce(nMiss, sumOp<label>());

        if (nMiss)
        {
            FatalIOErrorInFunction(dict_)
                << surfaceKey << " " << surfName << ": " << nMiss
                << " bed faces lie outside the plan view of the surface"
                << exit(FatalIOError);
        }

        Info<< "    Bed: " << surfaceKey << " from constant/triSurface/"
            << surfName << ", " << gMin(z) << " to " << gMax(z) << " m"
            << endl;
    }
    else if (dict_.found(valueKey))
    {
        const autoPtr<PatchFunction1<scalar>> func
        (
            PatchFunction1<scalar>::New(pp, valueKey, dict_)
        );
        const scalarField onPatch(func->value(bedMesh_.time().value()));

        forAll(faceLabels, i)
        {
            z[i] = onPatch[faceLabels[i] - pp.start()];
        }

        if (!func->constant())
        {
            Info<< "    Bed: " << valueKey << " (" << func->type() << "), "
                << gMin(z) << " to " << gMax(z) << " m" << endl;
        }
    }

    return tz;
}


void Foam::bedSediment::makeFixedMask()
{
    scalarField& fixed = fixedPtr_->primitiveFieldRef();
    fixed = 0;

    const faMesh& aMesh = *aMeshPtr_;
    const vectorField& Cf = aMesh.areaCentres().primitiveField();

    for (const entry& e : fixedBed_)
    {
        const word& patchName = e.keyword();
        const scalar width = fixedBed_.get<scalar>(patchName);

        const label patchi = mesh_.boundaryMesh().findPatchID(patchName);

        if (patchi < 0)
        {
            FatalIOErrorInFunction(fixedBed_)
                << "Unknown patch " << patchName << " in fixedBed"
                << exit(FatalIOError);
        }

        // Face centres of the patch on all processors, projected on z = 0
        List<pointField> allCentres(Pstream::nProcs());
        allCentres[Pstream::myProcNo()] =
            pointField(mesh_.boundaryMesh()[patchi].faceCentres());
        allCentres[Pstream::myProcNo()].replace(vector::Z, 0);
        Pstream::allGatherList(allCentres);

        pointField centres(ListListOps::combine<pointField>
        (
            allCentres, accessOp<pointField>()
        ));

        if (centres.empty())
        {
            continue;
        }

        treeBoundBox bb(centres);
        bb.inflate(1e-4);

        const indexedOctree<treeDataPoint> tree
        (
            treeDataPoint(centres),
            bb,
            10,
            3.0,
            3.0
        );

        // Horizontal distance from the floor-face centres to the patch
        forAll(Cf, facei)
        {
            const point p(Cf[facei].x(), Cf[facei].y(), 0);
            const pointIndexHit hit = tree.findNearest(p, sqr(width));

            if (hit.hit())
            {
                fixed[facei] = 1;
            }
        }
    }

    fixedPtr_->correctBoundaryConditions();

    Info<< "    Bed: " << returnReduce(label(sum(fixed)), sumOp<label>())
        << " fixed faces " << fixedBed_.toc() << endl;
}


void Foam::bedSediment::smoothShearStress()
{
    if (smoothShear_ < 1)
    {
        return;
    }

    const faMesh& aMesh = *aMeshPtr_;
    const labelUList& own = aMesh.edgeOwner();
    const labelUList& nei = aMesh.edgeNeighbour();
    const scalarField& active = activePtr_->primitiveField();
    areaVectorField& tauField = *tauPtr_;

    for (label pass = 0; pass < smoothShear_; ++pass)
    {
        const vectorField tau0(tauField.primitiveField());
        vectorField sum(tau0.size(), Zero);
        scalarField n(tau0.size(), Zero);

        // Active edge neighbours: internal edges, then coupled patches
        forAll(nei, edgei)
        {
            const label o = own[edgei];
            const label nb = nei[edgei];

            if (active[o] > 0.5 && active[nb] > 0.5)
            {
                sum[o] += tau0[nb];
                n[o] += 1;
                sum[nb] += tau0[o];
                n[nb] += 1;
            }
        }

        forAll(aMesh.boundary(), patchi)
        {
            const faPatch& p = aMesh.boundary()[patchi];

            if (!p.coupled())
            {
                continue;
            }

            const labelUList& ec = p.edgeFaces();
            const vectorField tn
            (
                tauField.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField an
            (
                activePtr_->boundaryField()[patchi].patchNeighbourField()
            );

            forAll(ec, i)
            {
                if (active[ec[i]] > 0.5 && an[i] > 0.5)
                {
                    sum[ec[i]] += tn[i];
                    n[ec[i]] += 1;
                }
            }
        }

        vectorField& tau = tauField.primitiveFieldRef();

        forAll(tau, facei)
        {
            if (n[facei] > 0)
            {
                tau[facei] = 0.5*tau0[facei] + 0.5*sum[facei]/n[facei];
            }
        }

        tauField.correctBoundaryConditions();
    }
}


void Foam::bedSediment::calcBedLoad()
{
    const areaScalarField& zb = *zbPtr_;
    const vectorField& tau = tauPtr_->primitiveField();
    vectorField& q = qPtr_->primitiveFieldRef();

    const areaVectorField gradZb(fac::grad(zb));

    const scalar sg = (s_ - 1)*g_*d50_;
    const scalar qScale = Foam::sqrt(sg)*d50_;

    // van Rijn particle parameter D* = d50 ((s - 1) g/nu^2)^(1/3)
    const scalar Dstar = d50_*Foam::cbrt((s_ - 1)*g_/sqr(nu_));

    const scalarField& iSeep = seepagePtr_->primitiveField();
    const vectorField& seepSlope = seepageSlopePtr_->primitiveField();

    slopeCoeff_.resize(q.size());
    slopeCoeff_ = 0;

    forAll(q, facei)
    {
        const scalar magTau = mag(tau[facei]);

        // Generalised Shields number (Francalanci et al. 2008: unity for
        // the other seepage models)
        scalar theta =
            magTau/sg*seepage_->thetaFactor(iSeep[facei]);
        if (maxShields_ > 0 && theta > maxShields_)
        {
            theta = maxShields_;
            ++windowClipped_;
        }

        if (magTau < VSMALL)
        {
            q[facei] = Zero;
            continue;
        }

        const vector tHat(tau[facei]/magTau);

        // Seepage: reduced effective weight (critical Shields number) and
        // friction coefficient
        const scalar mus =
            max
            (
                mus_*seepage_->reposeFactor(iSeep[facei])
               *liquefiedFactor(facei, liqReposeFactor_),
                scalar(1e-3)
            );
        const scalar thetaCFlat =
            thetaC0_*seepage_->thetaCFactor(iSeep[facei])
           *liquefiedFactor(facei, liqThetaCFactor_);
        const scalar betaMax = Foam::atan(mus);

        // Bed slope angle and angle between the flow and the downslope
        // direction (Roulund et al. 2005, Eq. 25); the horizontal seepage
        // force adds an apparent slope
        vector grad(gradZb[facei] + seepSlope[facei]);
        grad.z() = 0;
        const scalar magGrad = mag(grad);
        const scalar beta = min(Foam::atan(magGrad), betaMax);

        scalar cosPhi = 0;
        if (magGrad > VSMALL)
        {
            cosPhi = (tHat & (-grad/magGrad));
        }
        const scalar sin2Phi = max(1 - sqr(cosPhi), scalar(0));

        const scalar root =
            max(1 - sin2Phi*sqr(Foam::tan(beta)/mus), scalar(0));

        const scalar thetaC =
            max
            (
                thetaCFlat
               *(
                    Foam::cos(beta)*Foam::sqrt(root)
                  - cosPhi*Foam::sin(beta)/mus
                ),
                scalar(0)
            );

        scalar qStar = 0;

        if (theta > thetaC)
        {
            if (bedLoad_ == "meyerPeterMuller")
            {
                qStar = 8*pow(theta - thetaC, 1.5);
            }
            else if (bedLoad_ == "vanRijn")
            {
                // van Rijn (1984): transport stage T against the
                // slope-corrected threshold, bounded below by a tenth of the
                // flat-bed value so that T stays finite on steep down-slopes;
                // T >= 0 also where theta lies between the slope-corrected
                // threshold and that bound
                const scalar T =
                    max(theta/max(thetaC, 0.1*thetaCFlat) - 1, scalar(0));
                qStar = 0.053*pow(T, 2.1)/pow(Dstar, 0.3);
            }
            else
            {
                qStar =
                    18.74*(theta - thetaC)
                   *(Foam::sqrt(theta) - 0.7*Foam::sqrt(thetaC));
            }
        }

        q[facei] = qStar*qScale*tHat;

        if (qStar > 0 && slopeModel_ != "none")
        {
            slopeCoeff_[facei] = slopeGamma(theta)*qStar*qScale;
        }
    }

    if (!fixedBedLoad_)
    {
        // Fixed strips without sand supply (a rigid apron): no bed load
        const scalarField& fixed = fixedPtr_->primitiveField();
        forAll(q, facei)
        {
            if (fixed[facei] > 0.5)
            {
                q[facei] = Zero;
                slopeCoeff_[facei] = 0;
            }
        }
    }

    if (adaptLength_ > 0)
    {
        adaptBedLoad();
    }

    qPtr_->correctBoundaryConditions();
}


void Foam::bedSediment::adaptBedLoad()
{
    const faMesh& aMesh = *aMeshPtr_;
    const labelUList& own = aMesh.edgeOwner();
    const labelUList& nei = aMesh.edgeNeighbour();
    const vectorField& Le = aMesh.Le().primitiveField();
    const scalarField& S = aMesh.S().field();
    const scalarField& active = activePtr_->primitiveField();
    const vectorField& tau = tauPtr_->primitiveField();
    vectorField& q = qPtr_->primitiveFieldRef();

    // Transport direction: that of the bed shear stress
    areaVectorField eHat
    (
        IOobject
        (
            "bedLoadDirection",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        aMesh,
        dimensionedVector(dimless, Zero),
        "zeroGradient"
    );
    vectorField& e = eHat.primitiveFieldRef();
    forAll(e, facei)
    {
        const scalar magTau = mag(tau[facei]);
        if (magTau > VSMALL)
        {
            e[facei] = tau[facei]/magTau;
        }
    }
    eHat.correctBoundaryConditions();

    // Equilibrium magnitude
    const scalarField qEq(mag(q));

    // Discretisation of L (e . grad) q on face P with upwind inflow edges:
    //     (L/S_P) sum_in |F_e| (q_P - q_N),  F_e = e_e . Le_e
    // sumIn: sum_in |F_e| per face; bSource: sum_in |F_e| q_b over the
    // non-coupled boundary inflow edges (q_b = q_eq at open boundaries,
    // 0 at closed ones); no exchange with inactive faces
    scalarField sumIn(aMesh.nFaces(), Zero);
    scalarField bSource(aMesh.nFaces(), Zero);

    scalarField F(nei.size(), Zero);
    forAll(nei, edgei)
    {
        const label o = own[edgei];
        const label n = nei[edgei];

        if (active[o] < 0.5 || active[n] < 0.5)
        {
            continue;
        }

        F[edgei] = (0.5*(e[o] + e[n]) & Le[edgei]);

        if (F[edgei] > 0)
        {
            sumIn[n] += F[edgei];
        }
        else
        {
            sumIn[o] -= F[edgei];
        }
    }

    // Inflow edges of the coupled patches: -F > 0 per edge, else 0
    PtrList<scalarField> coupledIn(aMesh.boundary().size());

    forAll(aMesh.boundary(), patchi)
    {
        const faPatch& p = aMesh.boundary()[patchi];
        const labelUList& ec = p.edgeFaces();
        const vectorField& pLe = aMesh.Le().boundaryField()[patchi];

        if (p.coupled())
        {
            const vectorField en
            (
                eHat.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField an
            (
                activePtr_->boundaryField()[patchi].patchNeighbourField()
            );

            coupledIn.set(patchi, new scalarField(ec.size(), Zero));
            scalarField& cin = coupledIn[patchi];

            forAll(ec, i)
            {
                if (active[ec[i]] < 0.5 || an[i] < 0.5)
                {
                    continue;
                }

                const scalar Fb = (0.5*(e[ec[i]] + en[i]) & pLe[i]);

                if (Fb < 0)
                {
                    cin[i] = -Fb;
                    sumIn[ec[i]] -= Fb;
                }
            }
        }
        else
        {
            const bool open = isOpenBoundary(p);

            forAll(ec, i)
            {
                if (active[ec[i]] < 0.5)
                {
                    continue;
                }

                const scalar Fb = (e[ec[i]] & pLe[i]);

                if (Fb < 0)
                {
                    sumIn[ec[i]] -= Fb;

                    if (open)
                    {
                        bSource[ec[i]] -= Fb*qEq[ec[i]];
                    }
                }
            }
        }
    }

    // Jacobi sweeps of
    //     q_P (1 + L sumIn_P/S_P) - (L/S_P) sum_in |F_e| q_N = q_eq,P
    // (an M-matrix: the iteration converges for any L)
    const scalar L = adaptLength_;
    const scalarField diag(1 + L*sumIn/S);
    const scalarField b(qEq + L*bSource/S);

    areaScalarField qa
    (
        IOobject
        (
            "bedLoadAdapted",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        aMesh,
        dimensionedScalar(sqr(dimLength)/dimTime, Zero),
        "zeroGradient"
    );
    scalarField& qaf = qa.primitiveFieldRef();
    qaf = qEq;
    qa.correctBoundaryConditions();

    const scalar qRef = max(gMax(qEq), VSMALL);

    adaptIter_ = 0;
    while (adaptIter_ < adaptMaxIter_)
    {
        scalarField rhs(b);

        forAll(nei, edgei)
        {
            if (F[edgei] > 0)
            {
                const label n = nei[edgei];
                rhs[n] += L*F[edgei]*qaf[own[edgei]]/S[n];
            }
            else if (F[edgei] < 0)
            {
                const label o = own[edgei];
                rhs[o] -= L*F[edgei]*qaf[nei[edgei]]/S[o];
            }
        }

        forAll(aMesh.boundary(), patchi)
        {
            if (!coupledIn.set(patchi))
            {
                continue;
            }

            const labelUList& ec = aMesh.boundary()[patchi].edgeFaces();
            const scalarField& cin = coupledIn[patchi];
            const scalarField qn
            (
                qa.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(ec, i)
            {
                rhs[ec[i]] += L*cin[i]*qn[i]/S[ec[i]];
            }
        }

        const scalarField qNew(rhs/diag);
        const scalar change = gMax(mag(qNew - qaf)());

        qaf = qNew;
        qa.correctBoundaryConditions();
        ++adaptIter_;

        if (change < adaptTolerance_*qRef)
        {
            break;
        }
    }

    forAll(q, facei)
    {
        q[facei] = qaf[facei]*e[facei];
    }
}


Foam::scalarField Foam::bedSediment::netOutflow
(
    const areaVectorField& q,
    scalar& boundaryOutflow,
    const areaScalarField* donorScale,
    scalarField* outgoing
) const
{
    const faMesh& aMesh = *aMeshPtr_;
    const labelUList& own = aMesh.edgeOwner();
    const labelUList& nei = aMesh.edgeNeighbour();
    const vectorField& Le = aMesh.Le().primitiveField();

    scalarField out(aMesh.nFaces(), Zero);
    boundaryOutflow = 0;

    if (outgoing)
    {
        *outgoing = scalarField(aMesh.nFaces(), Zero);
    }
    const scalarField* scaleP =
        (donorScale ? &donorScale->primitiveField() : nullptr);

    const scalarField& active = activePtr_->primitiveField();

    // Down-slope transport: compact two-point gradient on the edges with
    // the edge-average coefficient gamma |q0| (zero on inactive faces)
    const bool slope = (slopeModel_ != "none");
    const areaScalarField& zb = *zbPtr_;
    const scalarField& z = zb.primitiveField();
    const scalarField& dc = aMesh.deltaCoeffs().primitiveField();
    const scalarField& magLe = aMesh.magLe().primitiveField();

    const tmp<areaScalarField> tkappa(slopeCoeffField());
    const areaScalarField& kappa = tkappa();

    // Internal edges: upwind in the direction of the edge flux; no flux
    // into or out of inactive (structure-covered) faces, none into a
    // face that is full up to the underside of a solid, none out of a
    // fixed face without sand supply (a rigid apron donates nothing, not
    // even down its edge into the scour hole next to it)
    const scalarField& capped = cappedPtr_->primitiveField();
    const scalarField& fixed = fixedPtr_->primitiveField();
    const bool noFixedDonor = !fixedBedLoad_;

    forAll(nei, edgei)
    {
        if (active[own[edgei]] < 0.5 || active[nei[edgei]] < 0.5)
        {
            continue;
        }

        const scalar fo = (q[own[edgei]] & Le[edgei]);
        const scalar fn = (q[nei[edgei]] & Le[edgei]);
        scalar F = (fo + fn > 0 ? fo : fn);

        if (slope)
        {
            const label o = own[edgei];
            const label n = nei[edgei];
            F -=
                0.5*(slopeCoeff_[o] + slopeCoeff_[n])
               *(z[n] - z[o])*dc[edgei]*magLe[edgei];
        }

        if ((F > 0 && capped[nei[edgei]] > 0.5) || (F < 0 && capped[own[edgei]] > 0.5))
        {
            continue;
        }

        if
        (
            noFixedDonor
         && (
                (F > 0 && fixed[own[edgei]] > 0.5)
             || (F < 0 && fixed[nei[edgei]] > 0.5)
            )
        )
        {
            continue;
        }

        if (scaleP)
        {
            F *= (*scaleP)[F > 0 ? own[edgei] : nei[edgei]];
        }
        if (outgoing)
        {
            (*outgoing)[F > 0 ? own[edgei] : nei[edgei]] += mag(F);
        }

        out[own[edgei]] += F;
        out[nei[edgei]] -= F;
    }

    // Boundary edges
    forAll(aMesh.boundary(), patchi)
    {
        const faPatch& p = aMesh.boundary()[patchi];
        const labelUList& ec = p.edgeFaces();
        const vectorField& pLe = aMesh.Le().boundaryField()[patchi];

        if (p.coupled())
        {
            const vectorField qn
            (
                q.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField an
            (
                activePtr_->boundaryField()[patchi].patchNeighbourField()
            );

            const scalarField zn
            (
                zb.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField kn
            (
                kappa.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField cn
            (
                cappedPtr_->boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField fn_
            (
                fixedPtr_->boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField& pdc =
                aMesh.deltaCoeffs().boundaryField()[patchi];
            const scalarField& pMagLe = aMesh.magLe().boundaryField()[patchi];
            scalarField sn(ec.size(), 1.0);
            if (donorScale)
            {
                sn = donorScale->boundaryField()[patchi].patchNeighbourField();
            }

            forAll(ec, i)
            {
                if (active[ec[i]] < 0.5 || an[i] < 0.5)
                {
                    continue;
                }

                const scalar fo = (q[ec[i]] & pLe[i]);
                const scalar fn = (qn[i] & pLe[i]);
                scalar F = (fo + fn > 0 ? fo : fn);

                if (slope)
                {
                    F -=
                        0.5*(slopeCoeff_[ec[i]] + kn[i])
                       *(zn[i] - z[ec[i]])*pdc[i]*pMagLe[i];
                }

                if ((F > 0 && cn[i] > 0.5) || (F < 0 && capped[ec[i]] > 0.5))
                {
                    continue;
                }

                if
                (
                    noFixedDonor
                 && ((F > 0 && fixed[ec[i]] > 0.5) || (F < 0 && fn_[i] > 0.5))
                )
                {
                    continue;
                }

                if (scaleP)
                {
                    F *= (F > 0 ? (*scaleP)[ec[i]] : sn[i]);
                }
                if (outgoing && F > 0)
                {
                    (*outgoing)[ec[i]] += F;
                }

                out[ec[i]] += F;
            }
        }
        else if (isOpenBoundary(p))
        {
            // Zero-gradient flux: equilibrium inflow, free outflow
            forAll(ec, i)
            {
                scalar F = (q[ec[i]] & pLe[i]);
                if (scaleP && F > 0)
                {
                    F *= (*scaleP)[ec[i]];
                }
                if (outgoing && F > 0)
                {
                    (*outgoing)[ec[i]] += F;
                }
                out[ec[i]] += F;
                boundaryOutflow += F;
            }
        }
        // Closed boundaries (walls, structures): no sediment flux
    }

    return out;
}


Foam::scalar Foam::bedSediment::liquefiedFactor
(
    const label facei,
    const scalar f
) const
{
    if (f == 1 || !liquefiedPtr_)
    {
        return 1;
    }

    const scalar d = liquefiedPtr_->primitiveField()[facei] - liqMinDepth_;

    if (d <= 0)
    {
        return 1;
    }

    const scalar w =
        liqTransition_ > 0 ? min(d/liqTransition_, scalar(1)) : scalar(1);

    return 1 - w*(1 - f);
}


Foam::scalar Foam::bedSediment::slopeGamma(const scalar theta) const
{
    if (slopeModel_ == "kochFlokstra")
    {
        return slopeAlpha_*Foam::sqrt(thetaC0_/max(theta, SMALL));
    }
    else if (slopeModel_ == "talmon")
    {
        return
            1.0
           /(
                9*Foam::pow(d50_/slopeDepth_, 0.3)
               *Foam::sqrt(max(theta, SMALL))
            );
    }

    return 0;
}


Foam::tmp<Foam::areaScalarField>
Foam::bedSediment::slopeCoeffField() const
{
    auto tkappa = tmp<areaScalarField>::New
    (
        IOobject
        (
            "bedSlopeCoeff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *aMeshPtr_,
        dimensionedScalar(dimArea/dimTime, Zero)
    );

    if (slopeModel_ != "none")
    {
        tkappa.ref().primitiveFieldRef() = slopeCoeff_;
        tkappa.ref().correctBoundaryConditions();
    }

    return tkappa;
}


Foam::scalar Foam::bedSediment::slopeTimeStep() const
{
    if (slopeModel_ == "none")
    {
        return GREAT;
    }

    // Positivity of the explicit diffusion: per face
    // dt sum_e kappa_e dc_e |L_e| <= 0.5 (1 - n) S, with kappa_e bounded by
    // the larger of the two face coefficients
    const faMesh& aMesh = *aMeshPtr_;
    const labelUList& own = aMesh.edgeOwner();
    const labelUList& nei = aMesh.edgeNeighbour();
    const scalarField& dc = aMesh.deltaCoeffs().primitiveField();
    const scalarField& magLe = aMesh.magLe().primitiveField();
    const scalarField& S = aMesh.S().field();

    const tmp<areaScalarField> tkappa(slopeCoeffField());

    scalarField sumK(aMesh.nFaces(), Zero);

    forAll(nei, edgei)
    {
        const scalar k =
            max(slopeCoeff_[own[edgei]], slopeCoeff_[nei[edgei]])
           *dc[edgei]*magLe[edgei];
        sumK[own[edgei]] += k;
        sumK[nei[edgei]] += k;
    }

    forAll(aMesh.boundary(), patchi)
    {
        const faPatch& p = aMesh.boundary()[patchi];
        if (!p.coupled())
        {
            continue;
        }
        const labelUList& ec = p.edgeFaces();
        const scalarField& pdc = aMesh.deltaCoeffs().boundaryField()[patchi];
        const scalarField& pMagLe = aMesh.magLe().boundaryField()[patchi];

        const scalarField kn
        (
            tkappa().boundaryField()[patchi].patchNeighbourField()
        );

        forAll(ec, i)
        {
            sumK[ec[i]] +=
                max(slopeCoeff_[ec[i]], kn[i])*pdc[i]*pMagLe[i];
        }
    }

    scalar dt = GREAT;
    forAll(sumK, facei)
    {
        if (sumK[facei] > VSMALL)
        {
            dt = min(dt, 0.5*(1 - porosity_)*S[facei]/sumK[facei]);
        }
    }

    return returnReduce(dt, minOp<scalar>());
}


Foam::label Foam::bedSediment::sandSlide(const scalar dtMorph)
{
    if (!sandSlide_)
    {
        return 0;
    }

    const faMesh& aMesh = *aMeshPtr_;
    const labelUList& own = aMesh.edgeOwner();
    const labelUList& nei = aMesh.edgeNeighbour();
    const scalarField& magLe = aMesh.magLe().primitiveField();
    const scalarField& dc = aMesh.deltaCoeffs().primitiveField();
    const scalarField& S = aMesh.S().field();

    areaScalarField& zb = *zbPtr_;
    scalarField& z = zb.primitiveFieldRef();

    // Stable pseudo-time step of the explicit diffusion: per face
    // dtau K0 sum_e |L_e| dc_e <= S (Gershgorin). The legacy step
    // 0.2 min(S)/K0 satisfies this on square faces (sum = 4) and is kept
    // there; on stretched faces (e.g. the ring next to a pier on an O-grid,
    // aspect ratio 5: sum ~ 10) it is not, and the slide grows a grid-scale
    // checkerboard
    const scalar minS = gMin(S);
    scalar dtau = 0.2*minS/slideK0_;
    {
        scalarField sumLdc(S.size(), Zero);
        forAll(nei, edgei)
        {
            const scalar c = magLe[edgei]*dc[edgei];
            sumLdc[own[edgei]] += c;
            sumLdc[nei[edgei]] += c;
        }
        forAll(aMesh.boundary(), patchi)
        {
            if (aMesh.boundary()[patchi].coupled())
            {
                const labelUList& ec = aMesh.boundary()[patchi].edgeFaces();
                const scalarField& pdc =
                    aMesh.deltaCoeffs().boundaryField()[patchi];
                const scalarField& pLe = aMesh.magLe().boundaryField()[patchi];
                forAll(ec, i)
                {
                    sumLdc[ec[i]] += pLe[i]*pdc[i];
                }
            }
        }

        scalar minRatio = GREAT;
        forAll(S, facei)
        {
            if (sumLdc[facei] > VSMALL)
            {
                minRatio = min(minRatio, S[facei]/sumLdc[facei]);
            }
        }
        const scalar dtauStable =
            0.8*returnReduce(minRatio, minOp<scalar>())/slideK0_;

        if (dtau > dtauStable*(1 + 1e-6))
        {
            dtau = dtauStable;
        }
    }

    // Friction coefficient per face, reduced by upward seepage (coupled
    // neighbours through the boundary field)
    areaScalarField musField
    (
        IOobject
        (
            "bedFriction",
            mesh_.time().timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedScalar(dimless, mus_),
        "zeroGradient"
    );
    if (seepage_->active() || liqReposeFactor_ != 1)
    {
        const scalarField& iSeep = seepagePtr_->primitiveField();
        scalarField& musFace = musField.primitiveFieldRef();
        forAll(musFace, facei)
        {
            musFace[facei] =
                mus_*seepage_->reposeFactor(iSeep[facei])
               *liquefiedFactor(facei, liqReposeFactor_);
        }
        musField.correctBoundaryConditions();
    }
    const scalarField& musFace = musField.primitiveField();

    // Rate-limited seepage failure: a face whose friction is reduced by
    // the seepage lowers at most by velocity x dt (a flow slide moves at
    // a finite rate) and, when liquefied, by the liquefied depth (only
    // that layer flows); a face that has shed its cap reverts to the full
    // repose friction for the rest of this slide
    const bool limitLiquefied =
        liquefiedSlideVelocity_ > 0 && liquefiedPtr_ && seepage_->active();
    scalarField lowerCap(z.size(), GREAT);
    if (limitLiquefied)
    {
        const scalarField& liq = liquefiedPtr_->primitiveField();
        const scalarField& iSeep = seepagePtr_->primitiveField();
        forAll(lowerCap, facei)
        {
            if (seepage_->reposeFactor(iSeep[facei]) < 1)
            {
                lowerCap[facei] = liquefiedSlideVelocity_*dtMorph;
                if (liq[facei] > 0)
                {
                    lowerCap[facei] = min(lowerCap[facei], liq[facei]);
                }
            }
        }
    }
    const scalarField z0Slide(z);

    // Apparent slope of the horizontal seepage force along each internal
    // edge (owner -> neighbour), added to the geometric slope
    const vectorField& seepSlope = seepageSlopePtr_->primitiveField();
    const vectorField& fc = aMesh.areaCentres().primitiveField();
    scalarField edgeSeep(nei.size(), Zero);
    if (seepage_->horizontalForce())
    {
        forAll(nei, edgei)
        {
            const vector d(fc[nei[edgei]] - fc[own[edgei]]);
            const vector sMean
            (
                0.5*(seepSlope[own[edgei]] + seepSlope[nei[edgei]])
            );
            edgeSeep[edgei] = (sMean & d)/max(mag(d), SMALL);
        }
    }

    // Faces at or inside a solid span cannot receive sand (the column is
    // full up to the underside of the solid); they can still give
    const areaScalarField& cappedField = *cappedPtr_;
    const scalarField& capped = cappedField.primitiveField();

    // Room below the next solid above each face, as a volume per pseudo
    // time step: the sand delivered into a face is clipped to it (the
    // rest stays with the giver); across processor edges the neighbour's
    // room is read once per iteration, so a small overshoot is possible
    // there (reported as blocked deposition)
    areaScalarField roomField
    (
        IOobject
        (
            "bedRoom",
            mesh_.time().timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedScalar(dimless, GREAT),
        "zeroGradient"
    );

    // Non-erodible floor: the volume a face may give per pseudo step,
    // shared equally among its edges (so the sum over the edges never
    // exceeds the sand above the floor); empty donors give nothing
    areaScalarField giveField
    (
        IOobject
        (
            "slideGive",
            mesh_.time().timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedScalar(dimless, GREAT),
        "zeroGradient"
    );
    scalarField nEdges(S.size(), Zero);
    if (rigidFloor_)
    {
        forAll(nei, edgei)
        {
            nEdges[own[edgei]] += 1;
            nEdges[nei[edgei]] += 1;
        }
        forAll(aMesh.boundary(), patchi)
        {
            if (aMesh.boundary()[patchi].coupled())
            {
                for (const label f : aMesh.boundary()[patchi].edgeFaces())
                {
                    nEdges[f] += 1;
                }
            }
        }
    }
    const scalar emptyDonor = 1e-12;

    label iter = 0;

    for (; iter < slideMaxIter_; ++iter)
    {
        scalarField out(z.size(), Zero);
        scalar maxExcess = 0;

        scalarField giveNbr;
        if (rigidFloor_)
        {
            giveField.primitiveFieldRef() =
                sedimentAboveRigid()*S/(dtau*max(nEdges, scalar(1)));
            giveField.correctBoundaryConditions();
        }
        const scalarField& give = giveField.primitiveField();

        const scalarField& active = activePtr_->primitiveField();
        const scalarField& fixed = fixedPtr_->primitiveField();

        updateCapped();

        scalarField roomLeft;
        scalarField roomNbr;
        if (haveSolids_)
        {
            roomField.primitiveFieldRef() = roomBelowSolid()*S/dtau;
            roomField.correctBoundaryConditions();
            roomLeft = roomField.primitiveField();
        }

        if (limitLiquefied)
        {
            const scalarField& iSeep = seepagePtr_->primitiveField();
            scalarField& mf = musField.primitiveFieldRef();
            forAll(mf, facei)
            {
                mf[facei] =
                    (z0Slide[facei] - z[facei] >= lowerCap[facei])
                  ? mus_
                  : mus_*seepage_->reposeFactor(iSeep[facei])
                   *liquefiedFactor(facei, liqReposeFactor_);
            }
            musField.correctBoundaryConditions();
        }

        // Internal edges: flux proportional to the slope above repose,
        // between active, erodible faces only
        forAll(nei, edgei)
        {
            if
            (
                active[own[edgei]] < 0.5 || active[nei[edgei]] < 0.5
             || fixed[own[edgei]] > 0.5 || fixed[nei[edgei]] > 0.5
            )
            {
                continue;
            }

            const scalar slope =
                (z[nei[edgei]] - z[own[edgei]])*dc[edgei] + edgeSeep[edgei];
            const scalar excess =
                mag(slope) - min(musFace[own[edgei]], musFace[nei[edgei]]);

            // The receiver (the low face) must have room below a solid
            if ((slope < 0 ? capped[nei[edgei]] : capped[own[edgei]]) > 0.5)
            {
                continue;
            }

            // Donor: the high face (owner if the slope is negative)
            const scalar giveMax =
                (slope < 0 ? give[own[edgei]] : give[nei[edgei]]);
            if (excess > 0 && giveMax > emptyDonor)
            {
                maxExcess = max(maxExcess, excess);

                // Flux from high to low: owner -> neighbour if slope < 0
                scalar F = -slideK0_*sign(slope)*excess*magLe[edgei];
                F = sign(F)*min(mag(F), giveMax);

                if (roomLeft.size())
                {
                    const label r = (F > 0 ? nei[edgei] : own[edgei]);
                    if (roomLeft[r] < GREAT)
                    {
                        const scalar magF = min(mag(F), roomLeft[r]);
                        roomLeft[r] -= magF;
                        F = sign(F)*magF;
                    }
                }

                out[own[edgei]] += F;
                out[nei[edgei]] -= F;
            }
        }

        // Processor edges
        forAll(aMesh.boundary(), patchi)
        {
            const faPatch& p = aMesh.boundary()[patchi];

            if (!p.coupled())
            {
                continue;
            }

            const labelUList& ec = p.edgeFaces();
            const scalarField zn
            (
                zb.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField& pdc =
                aMesh.deltaCoeffs().boundaryField()[patchi];
            const scalarField& pLe = aMesh.magLe().boundaryField()[patchi];
            const scalarField an
            (
                activePtr_->boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField fn
            (
                fixedPtr_->boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField mn
            (
                musField.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField cn
            (
                cappedField.boundaryField()[patchi].patchNeighbourField()
            );
            if (roomLeft.size())
            {
                roomNbr = roomField.boundaryField()[patchi].patchNeighbourField();
            }
            giveNbr = giveField.boundaryField()[patchi].patchNeighbourField();

            forAll(ec, i)
            {
                if
                (
                    active[ec[i]] < 0.5 || an[i] < 0.5
                 || fixed[ec[i]] > 0.5 || fn[i] > 0.5
                )
                {
                    continue;
                }

                const scalar slope = (zn[i] - z[ec[i]])*pdc[i];
                const scalar excess = mag(slope) - min(musFace[ec[i]], mn[i]);

                if ((slope < 0 ? cn[i] : capped[ec[i]]) > 0.5)
                {
                    continue;
                }

                // Donor: this face if the slope to the neighbour is negative
                const scalar giveMax = (slope < 0 ? give[ec[i]] : giveNbr[i]);
                if (excess > 0 && giveMax > emptyDonor)
                {
                    maxExcess = max(maxExcess, excess);

                    // Positive: out of this face, into the neighbour
                    scalar F = -slideK0_*sign(slope)*excess*pLe[i];
                    F = sign(F)*min(mag(F), giveMax);

                    if (roomLeft.size())
                    {
                        if (F > 0 && roomNbr[i] < GREAT)
                        {
                            F = min(F, roomNbr[i]);
                        }
                        else if (F < 0 && roomLeft[ec[i]] < GREAT)
                        {
                            const scalar magF = min(-F, roomLeft[ec[i]]);
                            roomLeft[ec[i]] -= magF;
                            F = -magF;
                        }
                    }

                    out[ec[i]] += F;
                }
            }
        }

        reduce(maxExcess, maxOp<scalar>());

        if (maxExcess < slideTolerance_*mus_)
        {
            break;
        }

        forAll(z, facei)
        {
            z[facei] -= dtau*out[facei]/S[facei];
        }

        zb.correctBoundaryConditions();
    }

    if (iter == slideMaxIter_)
    {
        WarningInFunction
            << "Sand slide not converged in " << slideMaxIter_
            << " iterations" << endl;
    }

    return iter;
}


Foam::tmp<Foam::scalarField> Foam::bedSediment::pointElevation() const
{
    const polyPatch& pp = bedMesh_.boundaryMesh()[patchID_];
    const labelList& faceLabels = aMeshPtr_->faceLabels();

    // z_b in patch-face order
    scalarField zbFace(pp.size(), Zero);
    forAll(faceLabels, i)
    {
        zbFace[faceLabels[i] - pp.start()] = zbPtr_->primitiveField()[i];
    }

    // Inverse-distance average of the surrounding faces, summed over all
    // processors so that points on processor boundaries get the same
    // elevation on each side (otherwise the merged surface has cracks
    // along them)
    const pointField& localPoints = pp.localPoints();
    const labelListList& pointFaces = pp.pointFaces();
    const vectorField faceCentres(pp.faceCentres());

    scalarField sumW(pp.nPoints(), Zero);
    scalarField sumWz(pp.nPoints(), Zero);

    forAll(pointFaces, pointi)
    {
        for (const label facei : pointFaces[pointi])
        {
            const scalar w =
                1/max(mag(faceCentres[facei] - localPoints[pointi]), VSMALL);
            sumW[pointi] += w;
            sumWz[pointi] += w*zbFace[facei];
        }
    }

    syncTools::syncPointList
    (
        bedMesh_, pp.meshPoints(), sumW, plusEqOp<scalar>(), scalar(0)
    );
    syncTools::syncPointList
    (
        bedMesh_, pp.meshPoints(), sumWz, plusEqOp<scalar>(), scalar(0)
    );

    auto tz = tmp<scalarField>::New(pp.nPoints());
    scalarField& z = tz.ref();
    forAll(z, pointi)
    {
        z[pointi] = sumWz[pointi]/max(sumW[pointi], VSMALL);
    }

    return tz;
}


void Foam::bedSediment::writeSurface() const
{
    const polyPatch& pp = bedMesh_.boundaryMesh()[patchID_];
    const labelList& faceLabels = aMeshPtr_->faceLabels();

    // faMesh face -> patch-local face
    labelList patchFace(faceLabels.size());
    forAll(faceLabels, i)
    {
        patchFace[i] = faceLabels[i] - pp.start();
    }

    // Fields in patch-face order
    scalarField zbFace(pp.size(), Zero);
    vectorField tauFace(pp.size(), Zero);
    vectorField qFace(pp.size(), Zero);
    scalarField activeFace(pp.size(), Zero);
    scalarField seepFace(pp.size(), Zero);
    scalarField liqFace(pp.size(), Zero);
    scalarField exFace(pp.size(), Zero);
    scalarField capFace(pp.size(), Zero);
    vectorField slopeFace(pp.size(), Zero);
    scalarField fixedFace(pp.size(), Zero);

    forAll(patchFace, i)
    {
        zbFace[patchFace[i]] = zbPtr_->primitiveField()[i];
        tauFace[patchFace[i]] = tauPtr_->primitiveField()[i];
        qFace[patchFace[i]] = qPtr_->primitiveField()[i];
        activeFace[patchFace[i]] = activePtr_->primitiveField()[i];
        seepFace[patchFace[i]] = seepagePtr_->primitiveField()[i];
        liqFace[patchFace[i]] = liquefiedPtr_->primitiveField()[i];
        if (exchangePtr_)
        {
            exFace[patchFace[i]] = exchangePtr_->primitiveField()[i];
        }
        capFace[patchFace[i]] = capacityPtr_->primitiveField()[i];
        slopeFace[patchFace[i]] = seepageSlopePtr_->primitiveField()[i];
        fixedFace[patchFace[i]] = fixedPtr_->primitiveField()[i];
    }

    // Points of the patch lifted to the bed elevation
    const scalarField zPoint(pointElevation());
    pointField points(pp.localPoints());
    forAll(points, pointi)
    {
        points[pointi].z() = zPoint[pointi];
    }

    const fileName outputDir
    (
        mesh_.time().globalPath()/"postProcessing"/"bed"
       /mesh_.time().timeName()
    );

    autoPtr<surfaceWriter> writer(surfaceWriter::New("vtk"));
    writer->isPointData(false);
    writer->open(points, pp.localFaces(), outputDir/"bed", Pstream::parRun());
    writer->beginTime(mesh_.time());
    writer->write("zb", zbFace);
    writer->write("tauBed", tauFace);
    writer->write("qBed", qFace);
    writer->write("bedActive", activeFace);
    writer->write("bedFixed", fixedFace);
    writer->write("bedSeepage", seepFace);
    writer->write("bedLiquefied", liqFace);
    writer->write("bedExchange", exFace);
    writer->write("bedCapacity", capFace);
    writer->write("bedSeepageSlope", slopeFace);
    writer->endTime();
    writer->close();
}


Foam::tmp<Foam::scalarField> Foam::bedSediment::sedimentAboveRigid() const
{
    const scalarField& z = zbPtr_->primitiveField();
    if (!rigidFloor_)
    {
        return tmp<scalarField>::New(z.size(), GREAT);
    }
    return max(z - zRigidPtr_->primitiveField(), scalar(0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bedSediment::bedSediment
(
    const fvMesh& mesh,
    const polyMesh& bedMesh,
    const dictionary& dict,
    const word& patchName
)
:
    mesh_(mesh),
    bedMesh_(bedMesh),
    dict_(dict),
    fixedBed_(dict.subOrEmptyDict("fixedBed")),
    openBoundaries_(dict.getOrDefault<wordList>("openBoundaries", wordList())),
    patchID_(bedMesh.boundaryMesh().findPatchID(patchName)),
    d50_(dict.get<scalar>("grainDiameter")),
    s_(dict.getOrDefault<scalar>("specificGravity", 2.65)),
    porosity_(dict.getOrDefault<scalar>("porosity", 0.4)),
    g_(dict.getOrDefault<scalar>("gravity", 9.81)),
    thetaC0_(dict.getOrDefault<scalar>("criticalShields", 0.05)),
    mus_(dict.getOrDefault<scalar>("frictionCoefficient", 0.63)),
    bedLoad_(dict.getOrDefault<word>("bedLoad", "engelundFredsoe")),
    nu_(dict.getOrDefault<scalar>("kinematicViscosity", 1e-6)),
    morphFactor_
    (
        Function1<scalar>::NewIfPresent("morphFactor", dict, word::null, &mesh)
    ),
    smoothShear_(dict.getOrDefault<label>("smoothShear", 1)),
    morphStartTime_(dict.getOrDefault<scalar>("morphStartTime", 0)),
    maxBedChange_(dict.getOrDefault<scalar>("maxBedChange", 1e-4)),
    bedCourant_(dict.getOrDefault<scalar>("bedCourant", 0)),
    maxStepChange_(dict.getOrDefault<scalar>("maxBedChangePerStep", 0)),
    morphFactorGrowth_(dict.getOrDefault<scalar>("morphFactorGrowth", 1.2)),
    morphFactorInForce_(0),
    morphTime_(0),
    lastAdvanceTime_(-GREAT),
    slopeModel_
    (
        dict.subOrEmptyDict("bedSlope").getOrDefault<word>("model", "none")
    ),
    slopeAlpha_
    (
        dict.subOrEmptyDict("bedSlope").getOrDefault<scalar>("alpha", 1.5)
    ),
    slopeDepth_
    (
        dict.subOrEmptyDict("bedSlope").getOrDefault<scalar>("flowDepth", 0)
    ),
    sandSlide_(dict.getOrDefault<bool>("sandSlide", true)),
    slideK0_(dict.getOrDefault<scalar>("slideDiffusivity", 1e-3)),
    slideTolerance_(dict.getOrDefault<scalar>("slideTolerance", 0.02)),
    slideMaxIter_(dict.getOrDefault<label>("slideMaxIter", 500)),
    liquefiedSlideVelocity_
    (
        dict.subOrEmptyDict("liquefiedSlide").getOrDefault<scalar>
        (
            "velocity", 0
        )
    ),
    liqMinDepth_
    (
        dict.subOrEmptyDict("liquefiedBed").getOrDefault<scalar>
        (
            "minDepth", 0
        )
    ),
    liqTransition_
    (
        dict.subOrEmptyDict("liquefiedBed").getOrDefault<scalar>
        (
            "transitionDepth", 1e-3
        )
    ),
    liqThetaCFactor_
    (
        dict.subOrEmptyDict("liquefiedBed").getOrDefault<scalar>
        (
            "criticalShieldsFactor", 1
        )
    ),
    liqReposeFactor_
    (
        dict.subOrEmptyDict("liquefiedBed").getOrDefault<scalar>
        (
            "reposeFactor", 1
        )
    ),
    adaptLength_(dict.getOrDefault<scalar>("adaptationLength", 0)),
    fixedBedLoad_
    (
        dict.getOrDefault<word>("fixedBedLoad", "capacity") != "none"
    ),
    adaptTolerance_
    (
        dict.getOrDefault<scalar>("adaptationTolerance", 1e-6)
    ),
    adaptMaxIter_(dict.getOrDefault<label>("adaptationMaxIter", 500)),
    fluxAveragingPeriod_(dict.getOrDefault<scalar>("fluxAveragingPeriod", 0)),
    fluxWindowTime_(0),
    haveMeanFlux_(false),
    windowMaxQ_(0),
    windowMaxQPoint_(Zero),
    maxShields_(dict.getOrDefault<scalar>("maxShields", 0)),
    windowClipped_(0),
    adaptIter_(0),
    zbRead_(false),
    cumulativeMassError_(0),
    cumulativeOverlap_(0),
    cumulativeBlocked_(0),
    seepage_(seepageModel::New(dict, (s_ - 1)*(1 - porosity_))),
    seepagePtr_(nullptr),
    liquefiedPtr_(nullptr),
    capacityPtr_(nullptr),
    seepageSlopePtr_(nullptr)
{
    if
    (
        bedLoad_ != "engelundFredsoe"
     && bedLoad_ != "meyerPeterMuller"
     && bedLoad_ != "vanRijn"
    )
    {
        FatalIOErrorInFunction(dict)
            << "Unknown bedLoad " << bedLoad_
            << "; use engelundFredsoe, meyerPeterMuller or vanRijn"
            << exit(FatalIOError);
    }

    if
    (
        slopeModel_ != "none"
     && slopeModel_ != "kochFlokstra"
     && slopeModel_ != "talmon"
    )
    {
        FatalIOErrorInFunction(dict)
            << "Unknown bedSlope model " << slopeModel_
            << "; use none, kochFlokstra or talmon" << exit(FatalIOError);
    }
    if (slopeModel_ == "talmon" && slopeDepth_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "bedSlope model talmon needs a positive flowDepth"
            << exit(FatalIOError);
    }

    makeFaMesh(patchName);
    const faMesh& aMesh = *aMeshPtr_;

    zbPtr_.reset
    (
        new areaScalarField
        (
            IOobject
            (
                "zb",
                mesh.time().timeName(),
                aMesh.thisDb(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedScalar(dimLength, Zero),
            "zeroGradient"
        )
    );
    zbRead_ = zbPtr_->headerOk();

    // A zb file of the undecomposed case that did not reach the processors:
    // decomposePar skips area fields without an area mesh on disk
    // (all ranks take the same branches: the collectives stay matched)
    if (Pstream::parRun() && returnReduceOr(!zbRead_))
    {
        const fileName global
        (
            mesh.time().globalPath()/mesh.time().timeName()
           /faMesh::prefix()/"zb"
        );
        if (returnReduceOr(isFile(global)))
        {
            FatalIOErrorInFunction(dict_)
                << global << " exists but the processor directories have no "
                << "zb: write the area mesh with makeFaMesh "
                << "(system/finite-area/faMeshDefinition) so that "
                << "decomposePar or makeFaMesh -parallel decomposes it"
                << exit(FatalIOError);
        }
    }

    // Initial bed (no zb field file): initialSurface, initialElevation
    // (constant, expression or mappedFile), else z_b = 0
    if (!zbRead_)
    {
        zbPtr_->primitiveFieldRef() = bedLevelInput
        (
            patchName, "initialElevation", "initialSurface", 0
        );
        zbPtr_->correctBoundaryConditions();
    }

    // Non-erodible floor: the field file zRigid, else the uniform
    // rigidElevation, else none
    {
        IOobject zrIO
        (
            "zRigid",
            mesh.time().timeName(),
            aMesh.thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE    // restarts read it back
        );
        const bool zrFile = zrIO.typeHeaderOk<areaScalarField>(true);
        rigidFloor_ =
            zrFile
         || dict.found("rigidElevation")
         || dict.found("rigidSurface");
        cumulativeRigidClamp_ = 0;

        if (rigidFloor_)
        {
            zRigidPtr_.reset
            (
                new areaScalarField
                (
                    zrIO,
                    aMesh,
                    dimensionedScalar(dimLength, Zero),
                    "zeroGradient"
                )
            );
            if (!zrFile)
            {
                zRigidPtr_->primitiveFieldRef() = bedLevelInput
                (
                    patchName, "rigidElevation", "rigidSurface", 0
                );
                zRigidPtr_->correctBoundaryConditions();
            }

            const scalarField& zr = zRigidPtr_->primitiveField();
            const scalarField& z = zbPtr_->primitiveField();
            label nBelow = 0;
            forAll(z, facei)
            {
                if (z[facei] < zr[facei] - SMALL)
                {
                    ++nBelow;
                }
            }
            reduce(nBelow, sumOp<label>());
            if (nBelow)
            {
                FatalIOErrorInFunction(dict)
                    << "Initial bed below the non-erodible floor on "
                    << nBelow << " faces" << exit(FatalIOError);
            }
            Info<< "    Bed: non-erodible floor "
                << (zrFile ? "from the field zRigid" : "from the bed dictionary")
                << ", z_r min/max " << gMin(zr) << " / " << gMax(zr)
                << ", sand above it max " << gMax((z - zr)()) << " m"
                << endl;
        }
    }

    tauPtr_.reset
    (
        new areaVectorField
        (
            IOobject
            (
                "tauBed",
                mesh.time().timeName(),
                aMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedVector(sqr(dimVelocity), Zero),
            "zeroGradient"
        )
    );

    qPtr_.reset
    (
        new areaVectorField
        (
            IOobject
            (
                "qBed",
                mesh.time().timeName(),
                aMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedVector(sqr(dimLength)/dimTime, Zero),
            "zeroGradient"
        )
    );
    if (fluxAveragingPeriod_ > 0)
    {
        // Window sum and mean (read back on a restart when present)
        for (const word& name : {word("qBedSum"), word("qBedMean")})
        {
            autoPtr<areaVectorField>& ptr =
                (name == "qBedSum" ? qSumPtr_ : qMeanPtr_);
            ptr.reset
            (
                new areaVectorField
                (
                    IOobject
                    (
                        name,
                        mesh.time().timeName(),
                        aMesh.thisDb(),
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    aMesh,
                    dimensionedVector(sqr(dimLength)/dimTime, Zero),
                    "zeroGradient"
                )
            );
        }
    }

    activePtr_.reset
    (
        new areaScalarField
        (
            IOobject
            (
                "bedActive",
                mesh.time().timeName(),
                aMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedScalar(dimless, 1),
            "zeroGradient"
        )
    );

    cappedPtr_.reset
    (
        new areaScalarField
        (
            IOobject
            (
                "bedCapped",
                mesh.time().timeName(),
                aMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh,
            dimensionedScalar(dimless, Zero),
            "zeroGradient"
        )
    );

    seepageSlopePtr_.reset
    (
        new areaVectorField
        (
            IOobject
            (
                "bedSeepageSlope",
                mesh.time().timeName(),
                aMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh,
            dimensionedVector(dimless, Zero),
            "zeroGradient"
        )
    );

    capacityPtr_.reset
    (
        new areaScalarField
        (
            IOobject
            (
                "bedCapacity",
                mesh.time().timeName(),
                aMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedScalar(dimPressure*dimLength, Zero),
            "zeroGradient"
        )
    );

    liquefiedPtr_.reset
    (
        new areaScalarField
        (
            IOobject
            (
                "bedLiquefied",
                mesh.time().timeName(),
                aMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedScalar(dimLength, Zero),
            "zeroGradient"
        )
    );

    seepagePtr_.reset
    (
        new areaScalarField
        (
            IOobject
            (
                "bedSeepage",
                mesh.time().timeName(),
                aMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedScalar(dimless, Zero),
            "zeroGradient"
        )
    );

    fixedPtr_.reset
    (
        new areaScalarField
        (
            IOobject
            (
                "bedFixed",
                mesh.time().timeName(),
                aMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh,
            dimensionedScalar(dimless, 0),
            "zeroGradient"
        )
    );

    makeFixedMask();

    // Morphological state of a restart
    {
        IOobject stateIO
        (
            "bedMorphology",
            mesh.time().timeName(),
            "uniform",
            mesh.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );
        if (stateIO.typeHeaderOk<IOdictionary>(true))
        {
            IOdictionary state(stateIO);
            morphTime_ = state.getOrDefault<scalar>("morphTime", 0);
            morphFactorInForce_ =
                state.getOrDefault<scalar>("morphFactorInForce", 0);
            fluxWindowTime_ = state.getOrDefault<scalar>("fluxWindowTime", 0);
            haveMeanFlux_ =
                fluxAveragingPeriod_ > 0
             && state.getOrDefault<bool>("haveMeanFlux", false)
             && qMeanPtr_->headerOk();
            Info<< "    Bed: morphological state read: t_m " << morphTime_
                << " s, factor in force " << morphFactorInForce_ << endl;
        }
    }

    if (bedCourant_ < 0 || morphFactorGrowth_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "bedCourant must be non-negative and morphFactorGrowth >= 1"
            << exit(FatalIOError);
    }

    const std::string vanRijnInfo
    (
        bedLoad_ == "vanRijn"
      ? " (D* " + Foam::name(d50_*Foam::cbrt((s_ - 1)*g_/sqr(nu_))) + ")"
      : ""
    );

    Info<< "    Bed: d50 " << d50_ << " m, bed load " << bedLoad_
        << vanRijnInfo.c_str()
        << ", morphFactor "
        << (morphFactor_ ? morphFactor_->type() : word("constant"))
        << " (" << morphFactor() << " at t = " << mesh_.time().value() << ")"
        << (bedCourant_ > 0 ? " as the maximum, bedCourant " : "")
        << (bedCourant_ > 0 ? Foam::name(bedCourant_) : word())
        << (
                fluxAveragingPeriod_ > 0
              ? ", flux averaged over " + Foam::name(fluxAveragingPeriod_) + " s"
              : word()
           )
        << (maxStepChange_ > 0 ? ", bed change per step at most " : "")
        << (maxStepChange_ > 0 ? Foam::name(maxStepChange_) + " m" : word())
        << ", open boundaries " << openBoundaries_
        << ", bed-slope transport " << slopeModel_
        << ", fixed strips " << fixedBed_.toc() << endl;

    if (adaptLength_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "adaptationLength must be non-negative" << exit(FatalIOError);
    }
    if (liqThetaCFactor_ != 1 || liqReposeFactor_ != 1)
    {
        Info<< "    Bed: liquefied bed (depth > " << liqMinDepth_
            << " m, ramp " << liqTransition_ << " m): critical Shields x "
            << liqThetaCFactor_ << ", repose friction x " << liqReposeFactor_
            << endl;
    }
    if (adaptLength_ > 0)
    {
        Info<< "    Bed: bed-load adaptation length " << adaptLength_
            << " m (" << adaptLength_/d50_ << " d50)" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::bedSediment::setElevation(const scalarField& zb)
{
    zbPtr_->primitiveFieldRef() = zb;
    zbPtr_->correctBoundaryConditions();
}


void Foam::bedSediment::setSeepageGradient(const scalarField& i)
{
    areaScalarField& seep = *seepagePtr_;

    if (i.size() != seep.size())
    {
        FatalErrorInFunction
            << "Expected " << seep.size() << " seepage gradients, got "
            << i.size() << exit(FatalError);
    }

    seep.primitiveFieldRef() = i;
    seep.correctBoundaryConditions();
}


void Foam::bedSediment::setSeepageForce
(
    const vectorField& gradHp,
    const scalar gammaS
)
{
    areaVectorField& slopeField = *seepageSlopePtr_;
    vectorField& slope = slopeField.primitiveFieldRef();

    if (gradHp.size() != slope.size())
    {
        FatalErrorInFunction
            << "Expected " << slope.size() << " seepage forces, got "
            << gradHp.size() << exit(FatalError);
    }

    if (!seepage_->horizontalForce())
    {
        slope = Zero;
        return;
    }

    const scalarField& iSeep = seepagePtr_->primitiveField();

    forAll(slope, facei)
    {
        const scalar f = max(seepage_->thetaCFactor(iSeep[facei]), 0.05);
        slope[facei] = gradHp[facei]/(gammaS*f);
        slope[facei].z() = 0;
    }

    slopeField.correctBoundaryConditions();
}


void Foam::bedSediment::setCapacityIndex(const scalarField& cap)
{
    areaScalarField& c = *capacityPtr_;

    if (cap.size() != c.size())
    {
        FatalErrorInFunction
            << "Expected " << c.size() << " capacity values, got "
            << cap.size() << exit(FatalError);
    }

    c.primitiveFieldRef() = cap;
    c.correctBoundaryConditions();
}


void Foam::bedSediment::setLiquefiedDepth(const scalarField& depth)
{
    areaScalarField& liq = *liquefiedPtr_;

    if (depth.size() != liq.size())
    {
        FatalErrorInFunction
            << "Expected " << liq.size() << " liquefied depths, got "
            << depth.size() << exit(FatalError);
    }

    liq.primitiveFieldRef() = depth;
    liq.correctBoundaryConditions();
}


void Foam::bedSediment::setSolidSpans(List<List<Pair<scalar>>>&& spans)
{
    if (spans.size() && spans.size() != aMeshPtr_->nFaces())
    {
        FatalErrorInFunction
            << "Expected " << aMeshPtr_->nFaces() << " span lists, got "
            << spans.size() << exit(FatalError);
    }

    solidSpans_ = std::move(spans);
    haveSolids_ = returnReduce(solidSpans_.size() > 0, orOp<bool>());
}


void Foam::bedSediment::updateCapped()
{
    if (!haveSolids_)
    {
        return;
    }

    const scalarField& z = zbPtr_->primitiveField();
    scalarField& capped = cappedPtr_->primitiveFieldRef();

    forAll(solidSpans_, facei)
    {
        capped[facei] = 0;
        for (const Pair<scalar>& span : solidSpans_[facei])
        {
            if (z[facei] >= span.first() - SMALL && z[facei] < span.second())
            {
                capped[facei] = 1;
                break;
            }
        }
    }

    cappedPtr_->correctBoundaryConditions();
}


Foam::tmp<Foam::scalarField> Foam::bedSediment::roomBelowSolid() const
{
    const scalarField& z = zbPtr_->primitiveField();
    auto troom = tmp<scalarField>::New(z.size(), GREAT);
    scalarField& room = troom.ref();

    forAll(solidSpans_, facei)
    {
        for (const Pair<scalar>& span : solidSpans_[facei])
        {
            if (z[facei] < span.first())
            {
                room[facei] = span.first() - z[facei];
                break;
            }
            if (z[facei] < span.second())
            {
                room[facei] = 0;
                break;
            }
        }
    }

    return troom;
}


void Foam::bedSediment::applySolidSpans
(
    const scalarField& zOld,
    scalarField& z,
    scalar& overlap,
    scalar& blocked
) const
{
    if (!haveSolids_)
    {
        return;
    }

    const scalarField& S = aMeshPtr_->S().field();

    forAll(solidSpans_, facei)
    {
        const List<Pair<scalar>>& spans = solidSpans_[facei];

        if (spans.empty() || z[facei] == zOld[facei])
        {
            continue;
        }

        if (z[facei] < zOld[facei])
        {
            // Erosion by the height zOld - z: the sediment comes from the
            // gaps between the spans, top-down from zOld; a span inside
            // the way is passed at no cost
            scalar remaining = zOld[facei] - z[facei];
            scalar zc = zOld[facei];

            forAllReverse(spans, i)
            {
                const scalar zBot = spans[i].first();
                const scalar zTop = spans[i].second();

                if (zTop <= zc - remaining)
                {
                    // The span lies below the eroded range
                    break;
                }

                if (zc > zTop)
                {
                    // Gap above the span
                    remaining -= zc - zTop;
                    zc = zTop;
                }

                if (zc > zBot)
                {
                    // Inside the span: pass through
                    overlap += (1 - porosity_)*S[facei]*(zc - zBot);
                    zc = zBot;
                }
            }

            z[facei] = zc - remaining;
        }
        else
        {
            // Deposition by the height z - zOld: fills the gap up to the
            // underside of the next span above, the rest is blocked
            scalar remaining = z[facei] - zOld[facei];
            scalar zc = zOld[facei];

            forAll(spans, i)
            {
                const scalar zBot = spans[i].first();
                const scalar zTop = spans[i].second();

                if (zc < zBot)
                {
                    // Gap below the span
                    const scalar gap = zBot - zc;
                    if (remaining <= gap)
                    {
                        zc += remaining;
                        remaining = 0;
                        break;
                    }
                    remaining -= gap;
                    zc = zBot;
                }

                if (zc >= zBot && zc < zTop)
                {
                    // At or inside the span: no room
                    break;
                }
            }

            blocked += (1 - porosity_)*S[facei]*remaining;
            z[facei] = zc;
        }
    }
}


void Foam::bedSediment::setSuspendedExchange(const scalarField& netDeposition)
{
    const faMesh& aMesh = *aMeshPtr_;

    if (!exchangePtr_)
    {
        exchangePtr_.reset
        (
            new areaScalarField
            (
                IOobject
                (
                    "bedExchange",
                    mesh_.time().timeName(),
                    aMesh.thisDb(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                aMesh,
                dimensionedScalar(dimVelocity, Zero)
            )
        );
    }

    areaScalarField& ex = *exchangePtr_;

    if (netDeposition.size() != ex.size())
    {
        FatalErrorInFunction
            << "Expected " << ex.size() << " exchange values, got "
            << netDeposition.size() << exit(FatalError);
    }

    ex.primitiveFieldRef() = netDeposition;
    ex.correctBoundaryConditions();
}


void Foam::bedSediment::advance(const scalar dt)
{
    const faMesh& aMesh = *aMeshPtr_;
    const scalarField& S = aMesh.S().field();
    areaScalarField& zb = *zbPtr_;
    scalarField& z = zb.primitiveFieldRef();

    const scalarField z0(z);

    smoothShearStress();

    const bool meanFlux = fluxAveragingPeriod_ > 0;
    if (meanFlux)
    {
        // Period-averaged bed load: the instantaneous flux of this flow
        // step enters the running mean of the window; the Exner equation
        // below takes the mean of the last completed window
        updateCapped();
        calcBedLoad();
        qSumPtr_->primitiveFieldRef() += dt*qPtr_->primitiveField();
        fluxWindowTime_ += dt;
        {
            const vectorField& qi = qPtr_->primitiveField();
            const vectorField& Cf = aMesh.areaCentres().primitiveField();
            label iMax = -1;
            scalar qMax = windowMaxQ_;
            forAll(qi, facei)
            {
                if (mag(qi[facei]) > qMax)
                {
                    qMax = mag(qi[facei]);
                    iMax = facei;
                }
            }
            const scalar gQMax = returnReduce(qMax, maxOp<scalar>());
            if (iMax >= 0 && qMax == gQMax)
            {
                windowMaxQPoint_ = Cf[iMax];
            }
            if (gQMax > windowMaxQ_)
            {
                windowMaxQ_ = gQMax;
                reduce(windowMaxQPoint_, maxMagSqrOp<point>());
            }
        }
        if (fluxWindowTime_ >= fluxAveragingPeriod_ - SMALL)
        {
            qMeanPtr_->primitiveFieldRef() =
                qSumPtr_->primitiveField()/fluxWindowTime_;
            qMeanPtr_->correctBoundaryConditions();
            qSumPtr_->primitiveFieldRef() = Zero;
            fluxWindowTime_ = 0;
            haveMeanFlux_ = true;
            Info<< "Bed: bed-load flux averaged over "
                << fluxAveragingPeriod_ << " s, max |q_mean| "
                << gMax(mag(qMeanPtr_->primitiveField())())
                << ", largest instantaneous |q| in the window "
                << windowMaxQ_ << " at " << windowMaxQPoint_;
            if (maxShields_ > 0)
            {
                Info<< ", faces clipped at theta " << maxShields_ << ": "
                    << returnReduce(windowClipped_, sumOp<label>());
            }
            Info<< endl;
            windowMaxQ_ = 0;
            windowMaxQPoint_ = Zero;
            windowClipped_ = 0;
        }
        if (!haveMeanFlux_)
        {
            // First window: the bed waits for its mean flux
            lastAdvanceTime_ = mesh_.time().value();
            return;
        }
    }
    const areaVectorField& qExner = (meanFlux ? *qMeanPtr_ : *qPtr_);

    // Exner equation, sub-cycled so that no face changes by more than
    // maxBedChange per sub-step. With bedCourant the morphological step
    // of this flow step is limited by the bed change per flow step
    // (adaptive factor), from the rate of the first sub-step.
    const scalar mfMax = morphFactor();
    scalar mf = mfMax;
    scalar dtMorph = dt*mf;
    scalar remaining = dtMorph;
    const scalarField faceSize(sqrt(S));
    scalar outflowVolume = 0;
    scalar exchangeVolume = 0;
    scalar overlapVolume = 0;
    scalar blockedVolume = 0;
    scalar blockedSlide = 0;
    scalar rigidClampVolume = 0;
    label nSub = 0;

    while (remaining > SMALL*dtMorph && nSub < 1000)
    {
        if (!meanFlux)
        {
            updateCapped();
            calcBedLoad();
        }

        scalar boundaryOutflow = 0;
        scalarField outgoing;
        scalarField out
        (
            netOutflow
            (
                qExner, boundaryOutflow, nullptr,
                (rigidFloor_ ? &outgoing : nullptr)
            )
        );

        // Rate of bed change [m/s]; fixed faces keep their elevation, and
        // their net outflow is exchanged with the outside (rigid apron).
        // Suspended load: net deposition D - E on the active, erodible
        // faces (grains settling on a fixed apron leave the balance)
        const scalarField& fixed = fixedPtr_->primitiveField();
        scalar fixedOutflow = 0;
        scalarField exApplied(S.size(), Zero);

        auto bedRate = [&](const scalarField& netOut)
        {
            scalarField r(-netOut/((1 - porosity_)*S));
            fixedOutflow = 0;
            forAll(r, facei)
            {
                if (fixed[facei] > 0.5)
                {
                    fixedOutflow += netOut[facei];
                    r[facei] = 0;
                }
            }
            if (exchangePtr_)
            {
                const scalarField& ex = exchangePtr_->primitiveField();
                const scalarField& active = activePtr_->primitiveField();
                forAll(r, facei)
                {
                    if (active[facei] > 0.5 && fixed[facei] < 0.5)
                    {
                        exApplied[facei] = ex[facei];
                        r[facei] += ex[facei]/(1 - porosity_);
                    }
                }
            }
            return r;
        };

        scalarField rate(bedRate(out));

        const scalar maxRate = gMax(mag(rate)());

        if (nSub == 0 && (bedCourant_ > 0 || maxStepChange_ > 0))
        {
            scalar mfLimit = mfMax;
            if (bedCourant_ > 0)
            {
                // Largest rate relative to the face size [1/s]
                const scalar maxRel = gMax((mag(rate)/faceSize)());
                if (maxRel > VSMALL)
                {
                    mfLimit = min(mfMax, bedCourant_/(maxRel*dt));
                }
            }
            if (maxStepChange_ > 0 && maxRate > VSMALL)
            {
                // Absolute cap on the bed change per flow step: on a
                // deforming mesh this is the mesh velocity the fluid is
                // made to follow at the bed
                mfLimit = min(mfLimit, maxStepChange_/(maxRate*dt));
            }
            if (morphFactorInForce_ > 0)
            {
                mfLimit = min(mfLimit, morphFactorGrowth_*morphFactorInForce_);
            }
            mf = max(mfLimit, SMALL);
            dtMorph = dt*mf;
            remaining = dtMorph;
        }

        scalar dts =
            (
                maxRate > VSMALL
              ? min(remaining, maxBedChange_/maxRate)
              : remaining
            );
        dts = min(dts, slopeTimeStep());

        // A face depositing under a solid stops at the underside: the
        // sub-step ends there (the face is then capped: no more inflow)
        if (haveSolids_)
        {
            const scalarField room(roomBelowSolid());
            scalar dtRoom = GREAT;
            forAll(rate, facei)
            {
                // A face already full up to the solid (no room) is capped:
                // it receives nothing and must not stall the sub-step
                if
                (
                    rate[facei] > VSMALL
                 && room[facei] > SMALL && room[facei] < GREAT
                )
                {
                    dtRoom = min(dtRoom, room[facei]/rate[facei]);
                }
            }
            dts = min(dts, returnReduce(dtRoom, minOp<scalar>()));
        }

        // Non-erodible floor: a face gives at most the sand above the
        // floor in this sub-step; its edge fluxes are scaled by the ratio
        // and the divergence is formed again (conservative per edge)
        if (rigidFloor_)
        {
            areaScalarField scale
            (
                IOobject
                (
                    "rigidDonorScale",
                    mesh_.time().timeName(),
                    aMesh.thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                aMesh,
                dimensionedScalar(dimless, 1),
                "zeroGradient"
            );
            scalarField& sc = scale.primitiveFieldRef();
            const scalarField avail
            (
                sedimentAboveRigid()*(1 - porosity_)*S
            );
            label nLimited = 0;
            forAll(sc, facei)
            {
                const scalar give = outgoing[facei]*dts;
                if (give > avail[facei] && fixed[facei] < 0.5)
                {
                    sc[facei] = max(avail[facei], scalar(0))/give;
                    ++nLimited;
                }
            }
            if (returnReduce(nLimited, sumOp<label>()))
            {
                scale.correctBoundaryConditions();
                out = netOutflow(qExner, boundaryOutflow, &scale);
                rate = bedRate(out);
            }
        }

        if (haveSolids_)
        {
            const scalarField zPrev(z);
            z += dts*rate;
            applySolidSpans(zPrev, z, overlapVolume, blockedVolume);
        }
        else
        {
            z += dts*rate;
        }

        if (rigidFloor_)
        {
            const scalarField& zr = zRigidPtr_->primitiveField();
            forAll(z, facei)
            {
                if (z[facei] < zr[facei])
                {
                    rigidClampVolume +=
                        (1 - porosity_)*(zr[facei] - z[facei])*S[facei];
                    z[facei] = zr[facei];
                }
            }
        }
        zb.correctBoundaryConditions();

        outflowVolume +=
            returnReduce(boundaryOutflow - fixedOutflow, sumOp<scalar>())*dts;
        exchangeVolume += gSum(exApplied*S)*dts;
        remaining -= dts;
        ++nSub;
    }

    if (remaining > SMALL*dtMorph)
    {
        WarningInFunction
            << "Exner sub-step limit reached: " << remaining/dtMorph*100
            << "% of the morphological time step not applied" << endl;
    }

    label nSlide = 0;
    if (haveSolids_)
    {
        const scalarField zPrev(z);
        nSlide = sandSlide(dtMorph);
        applySolidSpans(zPrev, z, overlapVolume, blockedSlide);
        zb.correctBoundaryConditions();
    }
    else
    {
        nSlide = sandSlide(dtMorph);
    }

    overlapVolume = returnReduce(overlapVolume, sumOp<scalar>());
    blockedVolume = returnReduce(blockedVolume, sumOp<scalar>());
    blockedSlide = returnReduce(blockedSlide, sumOp<scalar>());
    blockedVolume += blockedSlide;
    cumulativeOverlap_ += overlapVolume;
    cumulativeBlocked_ += blockedVolume;

    // Sediment-volume balance: change of bed volume + outflow = net
    // deposition from suspension; the passage of z_b through solid spans
    // is bed change without sediment, blocked deposition is sediment
    // without bed change
    rigidClampVolume = returnReduce(rigidClampVolume, sumOp<scalar>());
    cumulativeRigidClamp_ += rigidClampVolume;
    const scalar dVolume =
        gSum((1 - porosity_)*(z - z0)*S) + outflowVolume - exchangeVolume
      + overlapVolume + blockedVolume - rigidClampVolume;
    cumulativeMassError_ += dVolume;

    morphFactorInForce_ = mf;
    morphTime_ += dtMorph;
    lastAdvanceTime_ = mesh_.time().value();

    Info<< "Bed: MF " << mf << ", t_m " << morphTime_ << " s, "
        << nSub << " Exner sub-steps, " << nSlide
        << " slide iterations, z_b min/max " << gMin(z) << " / " << gMax(z)
        << ", max |q| " << gMax(mag(qExner.primitiveField())());
    if (seepage_->active())
    {
        Info<< ", max i " << gMax(seepagePtr_->primitiveField());
    }
    if (exchangePtr_)
    {
        Info<< ", net deposition " << exchangeVolume << " m3";
    }
    if (haveSolids_)
    {
        Info<< ", solid overlap " << overlapVolume
            << " m3 (cumulative " << cumulativeOverlap_
            << "), blocked deposition " << blockedVolume
            << " m3 (slide " << blockedSlide << ", cumulative "
            << cumulativeBlocked_ << ")";
    }
    if (rigidFloor_)
    {
        Info<< ", rigid-floor clamp " << rigidClampVolume
            << " m3 (cumulative " << cumulativeRigidClamp_ << ")";
    }
    Info<< ", volume error " << dVolume
        << " (cumulative " << cumulativeMassError_ << ")";
    if (adaptLength_ > 0)
    {
        Info<< ", adaptation iterations " << adaptIter_;
    }
    Info<< endl;
}


void Foam::bedSediment::writeMorphologyState() const
{
    const Time& runTime = mesh_.time();

    IOdictionary state
    (
        IOobject
        (
            "bedMorphology",
            runTime.timeName(),
            "uniform",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );
    state.set("morphTime", morphTime_);
    state.set("morphFactorInForce", morphFactorInForce_);
    state.set("morphFactorMax", morphFactor());
    state.set("bedCourant", bedCourant_);
    state.set("fluxWindowTime", fluxWindowTime_);
    state.set("haveMeanFlux", haveMeanFlux_);
    state.regIOobject::write();

    if (Pstream::master())
    {
        const fileName dir(runTime.globalPath()/"postProcessing"/"bed");
        mkDir(dir);
        const fileName file(dir/"morphTime.dat");
        const bool exists = isFile(file);
        OFstream os(file, IOstreamOption(), IOstreamOption::APPEND_ATE);
        if (!exists)
        {
            os  << "# flow time [s], morphological time t_m [s], "
                << "factor in force, maximum factor" << nl;
        }
        os  << runTime.timeName() << ' ' << morphTime_ << ' '
            << morphFactorInForce_ << ' ' << morphFactor() << nl;
    }
}


// ************************************************************************* //
