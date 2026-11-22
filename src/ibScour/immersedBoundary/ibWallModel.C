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
#include "ibWallModel.H"
#include "immersedBoundary.H"
#include "centredCPCCellToCellStencilObject.H"
#include "fvcGrad.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- Values of a vol field in the compact addressing of a cell stencil
//- (local cells, local boundary faces, then remote cells)
template<class Type>
static List<Type> flatField
(
    const mapDistribute& map,
    const UList<Type>& cellValues,
    const Type& boundaryValue
)
{
    List<Type> flat(map.constructSize(), boundaryValue);

    forAll(cellValues, celli)
    {
        flat[celli] = cellValues[celli];
    }

    map.distribute(flat);

    return flat;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::ibWallModel::nuIB() const
{
    const fvMesh& mesh = ib_.mesh();
    const labelList& ibc = ib_.ibCells();

    if (mesh.foundObject<volScalarField>("nu"))
    {
        const volScalarField& nu = mesh.lookupObject<volScalarField>("nu");
        return tmp<scalarField>::New(scalarField(nu.primitiveField(), ibc));
    }

    // Constant viscosity from transportProperties
    const dictionary& transportProperties =
        mesh.lookupObject<IOdictionary>("transportProperties");

    return tmp<scalarField>::New
    (
        ibc.size(),
        dimensionedScalar("nu", dimViscosity, transportProperties).value()
    );
}


bool Foam::ibWallModel::sampleNeighbours
(
    const label i,
    const UList<vector>& flatU,
    const vector& Uw,
    const scalar nu,
    scalar& uTau,
    vector& tHat,
    scalar& unRatio
) const
{
    const labelList& idx = nbrIdx_[i];
    const scalarList& w = nbrWeights_[i];
    const scalarList& yN = nbrY_[i];
    const vector& ni = nIB_[i];
    const scalar y = yIB_[i];

    uTau = 0;
    tHat = Zero;
    unRatio = 0;

    if (idx.empty())
    {
        return false;
    }

    scalar sumW = 0;
    vector utSum(Zero);

    forAll(idx, j)
    {
        const vector Ur(flatU[idx[j]] - Uw);
        const scalar un = (Ur & ni);
        const vector ut(Ur - un*ni);

        uTau += w[j]*law_->uTau(mag(ut), yN[j], nu);
        unRatio += w[j]*un*y/yN[j];
        utSum += w[j]*ut;
        sumW += w[j];
    }

    if (sumW < VSMALL)
    {
        uTau = 0;
        unRatio = 0;
        return false;
    }

    uTau /= sumW;
    unRatio /= sumW;

    const scalar magUt = mag(utSum);
    tHat = (magUt > VSMALL ? utSum/magUt : vector::zero);

    return true;
}


Foam::scalar Foam::ibWallModel::consistentNut
(
    const label i,
    const scalar uTau,
    const scalar nu,
    const UList<scalar>& flatNut
) const
{
    const labelList& idx = nbrIdx_[i];
    const scalarList& w = nbrWeights_[i];
    const scalarList& yN = nbrY_[i];
    const scalar y = yIB_[i];

    const scalar ksPlus = law_->ks()*uTau/nu;
    const scalar uPlusC = law_->uPlus(y*uTau/nu, ksPlus);

    // Effective viscosity that transmits u_tau^2 over [y_c, y_j] according
    // to the wall law, and the neighbours' eddy viscosity
    scalar sumW = 0;
    scalar nuSec = 0;
    scalar nutN = 0;

    forAll(idx, j)
    {
        const scalar dUPlus = law_->uPlus(yN[j]*uTau/nu, ksPlus) - uPlusC;

        nuSec += w[j]*uTau*(yN[j] - y)/max(dUPlus, SMALL);
        nutN += w[j]*max(flatNut[idx[j]], scalar(0));
        sumW += w[j];
    }

    nuSec /= sumW;
    nutN /= sumW;

    // Linear interpolation to the faces gives nu + (nu_t,c + nu_t,j)/2
    return max(2*(nuSec - nu) - nutN, scalar(0));
}


void Foam::ibWallModel::wallTurbulence
(
    const scalar uTau,
    const scalar y,
    const scalar nu,
    scalar& omega,
    scalar& nutLog
) const
{
    law_->wallTurbulence
    (
        uTau, y, nu, Cmu_, beta1_, Aplus_, omega, nutLog
    );
}


Foam::scalar Foam::ibWallModel::wallK
(
    const label i,
    const scalar uTau,
    const scalar nu,
    const scalar kEq,
    const UList<scalar>& flatK
) const
{
    if (flatK.empty() || nbrIdx_[i].empty())
    {
        return max(kEq, SMALL);
    }

    // Smooth (ks+ <= 2.25): k -> 0 at the wall, kept as the damped
    // equilibrium value. Transitional and rough (ks+ >= 5): dk/dn = 0
    // (Roulund et al. 2005, sec. 2.3.3), i.e. the live neighbours' k.
    // ks+ from u_tau, or (roughnessScale k) from u_k = Cmu^1/4 sqrt(k_N) of
    // the neighbours, which stays finite where u_tau vanishes (separation,
    // reattachment); cf. nutkRoughWallFunction, Xu & Liu (2021, eq. 10)
    const labelList& idx = nbrIdx_[i];
    const scalarList& w = nbrWeights_[i];

    scalar kN = 0;
    scalar sumW = 0;

    forAll(idx, j)
    {
        kN += w[j]*max(flatK[idx[j]], scalar(0));
        sumW += w[j];
    }

    kN /= max(sumW, VSMALL);

    const scalar uScale =
    (
        roughnessFromK_
      ? Foam::pow(Cmu_, 0.25)*Foam::sqrt(kN)
      : uTau
    );
    const scalar ksPlus = law_->ks()*uScale/nu;
    const scalar f =
        min
        (
            max
            (
                (ksPlus - ksPlusSmooth_)
               /max(ksPlusRough_ - ksPlusSmooth_, SMALL),
                scalar(0)
            ),
            scalar(1)
        );

    if (f <= 0)
    {
        return max(kEq, SMALL);
    }

    return max((1 - f)*kEq + f*kN, SMALL);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ibWallModel::ibWallModel
(
    const immersedBoundary& ib,
    const dictionary& dict
)
:
    ib_(ib),
    law_(ibWallLaw::New(dict.subDict("wallLaw"))),
    stressMode_(dict.getOrDefault<word>("mode", "stress") == "stress"),
    consistentNut_(dict.getOrDefault<bool>("consistentNut", true)),
    sampleBelow_(dict.getOrDefault<scalar>("sampleBelow", 0.25)),
    minLiveNeighbours_(dict.getOrDefault<label>("minLiveNeighbours", 0)),
    pinned_(),
    directionHeight_(dict.getOrDefault<scalar>("directionHeight", 0)),
    reportMaxShear_(dict.getOrDefault<bool>("reportMaxShear", false)),
    yTurbMin_(dict.getOrDefault<scalar>("yTurbMin", 0.05)),
    Cmu_(dict.getOrDefault<scalar>("Cmu", 0.09)),
    beta1_(dict.getOrDefault<scalar>("beta1", 0.075)),
    Aplus_(dict.getOrDefault<scalar>("Aplus", 26)),
    roughnessFromK_
    (
        dict.getOrDefault<word>("roughnessScale", "uTau") == "k"
    ),
    ksPlusSmooth_(dict.getOrDefault<scalar>("ksPlusSmooth", 2.25)),
    ksPlusRough_(dict.getOrDefault<scalar>("ksPlusRough", 5)),
    wallStrain_
    (
        dict.getOrDefault<word>("wallStrain", "model") == "wall"
    ),
    wallShearStress_
    (
        IOobject
        (
            "ibWallShearStress",
            ib.mesh().time().timeName(),
            ib.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ib.mesh(),
        dimensionedVector(sqr(dimVelocity), Zero)
    ),
    yPlusField_
    (
        IOobject
        (
            "ibYPlus",
            ib.mesh().time().timeName(),
            ib.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ib.mesh(),
        dimensionedScalar(dimless, Zero)
    )
{
    const word mode(dict.getOrDefault<word>("mode", "stress"));
    if (mode != "stress" && mode != "velocity")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown wall-model mode " << mode
            << "; use stress or velocity" << exit(FatalIOError);
    }
    Info<< "    Wall model mode: " << mode;
    if (stressMode_)
    {
        Info<< ", u_tau sampled at the neighbours below y = "
            << sampleBelow_ << " h";
    }
    Info<< endl;

    if (dict.found("samplingFactor"))
    {
        WarningInFunction
            << "Entry 'samplingFactor' is obsolete and ignored: the wall law"
            << " is applied at the live neighbours of each IB cell" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ibWallModel::updateGeometry()
{
    const fvMesh& mesh = ib_.mesh();
    const labelList& ibc = ib_.ibCells();
    const pointField& C = mesh.C().primitiveField();
    const scalarField& d = ib_.distance().primitiveField();
    const vectorField& n = ib_.normal().primitiveField();

    const extendedCentredCellToCellStencil& stencilObj =
        centredCPCCellToCellStencilObject::New(mesh);
    const mapDistribute& map = stencilObj.map();
    const labelListList& stencil = stencilObj.stencil();

    // Centres, status and wall distance in compact stencil addressing.
    // Boundary-face entries are excluded (status -1).
    const List<point> flatC(flatField<point>(map, C, point::max));

    labelList cellStatus(ib_.status());
    const List<label> flatS(flatField<label>(map, cellStatus, label(-1)));
    const List<scalar> flatD(flatField<scalar>(map, d, scalar(0)));

    yIB_.resize(ibc.size());
    hIB_.resize(ibc.size());
    nIB_.resize(ibc.size());
    wallPoints_.resize(ibc.size());
    nbrIdx_.resize(ibc.size());
    nbrWeights_.resize(ibc.size());
    nbrY_.resize(ibc.size());

    label nEmpty = 0;

    forAll(ibc, i)
    {
        const label celli = ibc[i];
        const vector& ni = n[celli];
        const scalar h = ib_.cellExtent(celli, ni);

        yIB_[i] = max(d[celli], SMALL);
        hIB_[i] = h;
        nIB_[i] = ni;
        wallPoints_[i] = C[celli] - d[celli]*ni;

        // Live neighbours clearly further from the wall than the IB cell
        // (a quarter of the cell extent); relax if there are none
        DynamicList<label> idx;
        DynamicList<scalar> w;
        DynamicList<scalar> y;

        for (const scalar minRise : {0.25*h, SMALL})
        {
            for (const label compacti : stencil[celli])
            {
                if
                (
                    flatS[compacti] == immersedBoundary::LIVE
                 && flatD[compacti] > yIB_[i] + minRise
                )
                {
                    const vector r(flatC[compacti] - C[celli]);
                    const scalar magSqrR = magSqr(r);
                    const scalar cosTheta = (ni & r)/Foam::sqrt(magSqrR);

                    // Directional weight with a floor: a cell under an
                    // overhanging structure may have live neighbours only
                    // beside it (cosTheta = 0), and a zero weight sum
                    // divided the sampled values by zero
                    idx.push_back(compacti);
                    w.push_back(max(sqr(cosTheta), scalar(0.01))/magSqrR);
                    y.push_back(flatD[compacti]);
                }
            }

            if (idx.size())
            {
                break;
            }
        }

        if (idx.empty())
        {
            ++nEmpty;
        }

        nbrIdx_[i].transfer(idx);
        nbrWeights_[i].transfer(w);
        nbrY_[i].transfer(y);
    }

    Info<< "    Wall model: "
        << returnReduce(ibc.size(), sumOp<label>()) << " IB cells, "
        << returnReduce(nEmpty, sumOp<label>())
        << " without live neighbours" << endl;

    updateDirectionSampling();

    // Size the results
    Uib_.resize(ibc.size(), Zero);
    uTau_.resize(ibc.size(), Zero);
    yPlus_.resize(ibc.size(), Zero);
    kIB_.resize(ibc.size(), SMALL);
    omegaIB_.resize(ibc.size(), SMALL);
    nutIB_.resize(ibc.size(), Zero);
    nutWeight_.resize(ibc.size(), Zero);
    tauW_.resize(ibc.size(), Zero);
    wallCoeff_.resize(ibc.size(), Zero);
    Uw_.resize(ibc.size(), Zero);
}


void Foam::ibWallModel::updateDirectionSampling()
{
    dirIdx_.clear();
    dirMapPtr_.clear();

    if (directionHeight_ <= 0)
    {
        return;
    }

    const fvMesh& mesh = ib_.mesh();
    const labelList& ibc = ib_.ibCells();
    const pointField& C = mesh.C().primitiveField();
    const globalIndex globalCells(mesh.nCells());

    // Sample points and local search
    labelList globalId(ibc.size(), -1);
    DynamicList<point> unresolved;
    DynamicList<label> unresolvedIdx;

    forAll(ibc, i)
    {
        if (yIB_[i] >= directionHeight_)
        {
            continue;
        }

        const point p
        (
            C[ibc[i]] + (directionHeight_ - yIB_[i])*nIB_[i]
        );
        const label celli = mesh.findCell(p, polyMesh::CELL_TETS);

        if (celli >= 0)
        {
            globalId[i] = globalCells.toGlobal(celli);
        }
        else
        {
            unresolved.push_back(p);
            unresolvedIdx.push_back(i);
        }
    }

    // Points not found locally: every processor searches all of them
    if (Pstream::parRun())
    {
        List<pointField> allPts(Pstream::nProcs());
        allPts[Pstream::myProcNo()] = pointField(unresolved);
        Pstream::allGatherList(allPts);

        labelList offsets(Pstream::nProcs() + 1, 0);
        forAll(allPts, proci)
        {
            offsets[proci + 1] = offsets[proci] + allPts[proci].size();
        }

        labelList found(offsets.last(), -1);
        forAll(allPts, proci)
        {
            if (proci == Pstream::myProcNo())
            {
                continue;
            }

            forAll(allPts[proci], j)
            {
                const label celli =
                    mesh.findCell(allPts[proci][j], polyMesh::CELL_TETS);
                if (celli >= 0)
                {
                    found[offsets[proci] + j] = globalCells.toGlobal(celli);
                }
            }
        }

        Pstream::listCombineReduce(found, maxOp<label>());

        forAll(unresolvedIdx, j)
        {
            globalId[unresolvedIdx[j]] =
                found[offsets[Pstream::myProcNo()] + j];
        }
    }

    // Map to fetch the sampled velocities; cells without a sample point
    // (above the height, or outside the mesh) keep their own velocity
    label nSampled = 0;
    labelList elements(ibc.size());
    forAll(elements, i)
    {
        if (globalId[i] >= 0)
        {
            elements[i] = globalId[i];
            ++nSampled;
        }
        else
        {
            elements[i] = globalCells.toGlobal(ibc[i]);
        }
    }

    List<Map<label>> compactMap(Pstream::nProcs());
    dirMapPtr_.reset(new mapDistribute(globalCells, elements, compactMap));
    dirIdx_ = elements;

    Info<< "    Wall model: shear direction sampled at "
        << directionHeight_ << " m for "
        << returnReduce(nSampled, sumOp<label>()) << " IB cells" << endl;
}


void Foam::ibWallModel::correct(const volVectorField& U)
{
    pinned_.clear();
    addProfiling(correctProf, "ibWallModel::correct");

    const fvMesh& mesh = ib_.mesh();
    const labelList& ibc = ib_.ibCells();

    const mapDistribute& map =
        centredCPCCellToCellStencilObject::New(mesh).map();

    const List<vector> flatU
    (
        flatField<vector>(map, U.primitiveField(), vector::zero)
    );

    // Eddy viscosity in compact stencil addressing (absent for laminar)
    const volScalarField* nutPtr =
        mesh.cfindObject<volScalarField>("nut");

    const List<scalar> flatNut
    (
        nutPtr && consistentNut_
      ? flatField<scalar>(map, nutPtr->primitiveField(), scalar(0))
      : List<scalar>()
    );

    const scalarField nu(nuIB());

    // Turbulent kinetic energy in compact stencil addressing, for the
    // zero-gradient (Neumann) k of rough walls (Roulund et al. 2005)
    const volScalarField* kPtr = mesh.cfindObject<volScalarField>("k");

    const List<scalar> flatK
    (
        kPtr && law_->ks() > 0
      ? flatField<scalar>(map, kPtr->primitiveField(), scalar(0))
      : List<scalar>()
    );

    // Wall velocity at the wall points
    const vectorField Uw(ib_.solidVelocity(ibc, wallPoints_));

    const scalar ks = law_->ks();

    Uw_ = Uw;

    // Velocities at the direction-sampling height (compact addressing)
    const List<vector> flatUs
    (
        dirMapPtr_
      ? flatField<vector>(*dirMapPtr_, U.primitiveField(), Zero)
      : List<vector>()
    );

    auto sampledDirection = [&](const label i, const vector& tHat) -> vector
    {
        if (!dirMapPtr_ || dirIdx_[i] < 0)
        {
            return tHat;
        }
        const vector Ur(flatUs[dirIdx_[i]] - Uw[i]);
        const vector ut(Ur - (Ur & nIB_[i])*nIB_[i]);
        const scalar magUt = mag(ut);
        return (magUt > VSMALL ? ut/magUt : tHat);
    };

    if (stressMode_)
    {
        const vectorField& Ui = U.primitiveField();

        forAll(ibc, i)
        {
            const vector& ni = nIB_[i];
            const scalar y = yIB_[i];

            const vector Ur(Ui[ibc[i]] - Uw[i]);
            const vector ut(Ur - (Ur & ni)*ni);
            const scalar magUt = mag(ut);

            // Friction velocity and direction of the cell itself
            scalar uTau = law_->uTau(magUt, y, nu[i]);
            vector tHat(magUt > VSMALL ? ut/magUt : vector::zero);

            // Close to the wall: blend towards the neighbour-sampled values
            const scalar yBlend = sampleBelow_*hIB_[i];
            scalar uTauN = 0;
            vector tHatN(Zero);
            scalar unRatio = 0;

            // Weight of the wall-consistent eddy viscosity (correctNut)
            scalar nutWeight = 0;
            scalar nutSample = 0;

            if (y < yBlend)
            {
                // The cell's own wall law at y_c -> 0 diverges. In the
                // viscous sublayer u_tau = sqrt(nu |u_t|/y_c), which is exact
                // for laminar flow; a rough log law has no solution once
                // y_c < Cs ks/E and grows without bound. Limit the growth
                // to the viscous scaling from the value at yBlend
                const scalar uTauBlend = law_->uTau(magUt, yBlend, nu[i]);

                if (nbrIdx_[i].empty())
                {
                    // No live neighbours to sample (e.g. a thin fluid
                    // sliver): use the value at yBlend alone
                    uTau = uTauBlend;
                }
                else
                {
                    uTau =
                        min
                        (
                            uTau,
                            uTauBlend*Foam::sqrt(yBlend/max(y, VSMALL))
                        );
                }
            }

            if
            (
                y < yBlend
             && !nbrIdx_[i].empty()
             && sampleNeighbours(i, flatU, Uw[i], nu[i], uTauN, tHatN, unRatio)
            )
            {
                const scalar lambda = y/yBlend;
                uTau = lambda*uTau + (1 - lambda)*uTauN;

                const vector t(lambda*tHat + (1 - lambda)*tHatN);
                const scalar magT = mag(t);
                tHat = (magT > VSMALL ? t/magT : tHatN);

                // Near the wall the turbulence model's eddy viscosity of the
                // IB cell (y+ small) under-transmits the wall stress to the
                // live neighbours; blend towards the wall-law value
                if (flatNut.size() && uTau > VSMALL)
                {
                    nutWeight = 1 - lambda;
                    nutSample = consistentNut(i, uTau, nu[i], flatNut);
                }
            }

            const scalar yPlus = y*uTau/nu[i];

            // Wall-law velocity at the cell centre bounds the wall coefficient
            const scalar uRef =
                uTau*law_->uPlus(yPlus, ks*uTau/nu[i]);

            // Turbulence values at a regularised distance
            scalar omega = 0;
            scalar nutLog = 0;
            wallTurbulence
            (
                uTau,
                max(y, yTurbMin_*hIB_[i]),
                nu[i],
                omega,
                nutLog
            );

            uTau_[i] = uTau;
            yPlus_[i] = yPlus;
            Uib_[i] = Ui[ibc[i]];
            omegaIB_[i] = omega;
            nutIB_[i] = nutSample;
            nutWeight_[i] = nutWeight;
            kIB_[i] = wallK(i, uTau, nu[i], nutLog*omega, flatK);
            // Momentum sink linearised with the cell's own wall law, so
            // that it stays proportional to the cell velocity and grows as
            // the bed rises to the cell centre (nu/y for the laminar law).
            // The sampled (blended) u_tau above serves the bed shear only:
            // used in the sink it gave a constant-stress force that no
            // longer held a cell being buried (10 -> 90 000 m/s, wave
            // basin 2026-11-01)
            const scalar uTauCell = law_->uTau(magUt, y, nu[i]);
            wallCoeff_[i] = sqr(uTauCell)/max(max(magUt, uRef), VSMALL);
            tauW_[i] = sqr(uTau)*sampledDirection(i, tHat);

            // Boxed-in cell: wall-law velocity instead of its own momentum
            if (label(nbrIdx_[i].size()) < minLiveNeighbours_)
            {
                scalar uTauN = 0;
                vector tHatN(Zero);
                scalar unRatioN = 0;
                if (sampleNeighbours(i, flatU, Uw[i], nu[i], uTauN, tHatN, unRatioN))
                {
                    const scalar uPlusN =
                        law_->uPlus(y*uTauN/nu[i], ks*uTauN/nu[i]);
                    Uib_[i] = Uw[i] + uTauN*uPlusN*tHatN + unRatioN*ni;
                    uTau_[i] = uTauN;
                    tauW_[i] = sqr(uTauN)*sampledDirection(i, tHatN);
                }
                else
                {
                    Uib_[i] = Uw[i];
                    uTau_[i] = 0;
                    tauW_[i] = Zero;
                }
                wallCoeff_[i] = 0;
                pinned_.append(i);
            }
        }
    }
    else forAll(ibc, i)
    {
        const vector& ni = nIB_[i];
        const scalar y = yIB_[i];

        // Wall law at each live neighbour, at its own wall distance
        scalar uTau = 0;
        vector tHat(Zero);
        scalar unRatio = 0;

        if (!sampleNeighbours(i, flatU, Uw[i], nu[i], uTau, tHat, unRatio))
        {
            Uib_[i] = Uw[i];
            uTau_[i] = 0;
            yPlus_[i] = 0;
            nutIB_[i] = 0;
            kIB_[i] = SMALL;
            omegaIB_[i] = 6*nu[i]/(beta1_*sqr(y));
            tauW_[i] = Zero;
            continue;
        }

        // Velocity at the IB cell centre
        const scalar yPlus = y*uTau/nu[i];
        const scalar ksPlus = ks*uTau/nu[i];
        const scalar uPlusC = law_->uPlus(yPlus, ksPlus);

        Uib_[i] = Uw[i] + uTau*uPlusC*tHat + unRatio*ni;

        // Turbulence quantities at the IB cell centre
        scalar omega = 0;
        scalar nutLog = 0;
        wallTurbulence(uTau, y, nu[i], omega, nutLog);

        scalar nut = nutLog;

        if (flatNut.size() && uTau > VSMALL)
        {
            nut = consistentNut(i, uTau, nu[i], flatNut);
        }

        uTau_[i] = uTau;
        yPlus_[i] = yPlus;
        nutIB_[i] = nut;
        omegaIB_[i] = omega;
        kIB_[i] = wallK(i, uTau, nu[i], nutLog*omega, flatK);
        tauW_[i] = sqr(uTau)*sampledDirection(i, tHat);
    }

    // Output fields
    vectorField& tw = wallShearStress_.primitiveFieldRef();
    scalarField& yp = yPlusField_.primitiveFieldRef();
    tw = Zero;
    yp = 0;

    forAll(ibc, i)
    {
        tw[ibc[i]] = tauW_[i];
        yp[ibc[i]] = yPlus_[i];
    }

    if (reportMaxShear_)
    {
        reportMaxShear(U, flatU, Uw, nu);
    }

    if (immersedBoundary::debug)
    {
        Info<< "IB wall model: u_tau min/max = "
            << gMin(uTau_) << " / " << gMax(uTau_)
            << ", y+ min/max = " << gMin(yPlus_) << " / " << gMax(yPlus_)
            << endl;
    }
}


void Foam::ibWallModel::wallStrainNut(volScalarField& nut) const
{
    const fvMesh& mesh = ib_.mesh();
    const labelList& ibc = ib_.ibCells();

    const volVectorField& U = mesh.lookupObject<volVectorField>("U");
    const volScalarField* kPtr = mesh.cfindObject<volScalarField>("k");
    const volScalarField* omegaPtr =
        mesh.cfindObject<volScalarField>("omega");

    // Only for the k-omega SST closure (k and omega present)
    if (!kPtr || !omegaPtr)
    {
        return;
    }

    const volScalarField& k = *kPtr;
    const volScalarField& omega = *omegaPtr;

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const surfaceScalarField& w = mesh.weights();
    const scalarField nu(nuIB());

    // Solid cells, to find the blocked faces of the IB cells
    boolList isSolid(mesh.nCells(), false);

    for (const label celli : ib_.solidCells())
    {
        isSolid[celli] = true;
    }

    // SST constants (OpenFOAM defaults)
    const scalar a1 = 0.31;
    const scalar betaStar = 0.09;

    scalarField& nuti = nut.primitiveFieldRef();

    forAll(ibc, i)
    {
        const label celli = ibc[i];
        const scalar y = yIB_[i];

        // Gauss gradient with the wall velocity on the blocked faces
        tensor G(gradU[celli]);
        const cell& c = mesh.cells()[celli];

        forAll(c, fi)
        {
            const label facei = c[fi];

            if (!mesh.isInternalFace(facei))
            {
                continue;
            }

            const bool isOwn = (own[facei] == celli);
            const label nbr = (isOwn ? nei[facei] : own[facei]);

            if (!isSolid[nbr])
            {
                continue;
            }

            // Outward normal of the IB cell. The face is blocked by the
            // solid neighbour in full, as in the momentum and diffusion
            // masks (the geometric open area of bed-plane faces is not
            // reliable: it is 1 there)
            const vector Sout = (isOwn ? 1 : -1)*mesh.Sf()[facei];
            const scalar blocked = 1;

            const vector Uf =
                w[facei]*U[own[facei]] + (1 - w[facei])*U[nei[facei]];

            // Replace the blocked share of the face flux Sout (x) Uf by
            // Sout (x) Uw
            G += blocked*(Sout*(Uw_[i] - Uf))/mesh.V()[celli];
        }

        // Strain S = sqrt(2 S:S) with S = symm(grad U) of the cell
        const symmTensor Ssym(symm(G));
        const scalar S = Foam::sqrt(2*magSqr(Ssym));

        // SST F2 at the wall distance of the IB cell
        const scalar kc = max(k[celli], SMALL);
        const scalar omc = max(omega[celli], SMALL);
        const scalar arg2 = min
        (
            max
            (
                2*Foam::sqrt(kc)/(betaStar*omc*max(y, SMALL)),
                500*nu[i]/(sqr(max(y, SMALL))*omc)
            ),
            scalar(100)
        );
        const scalar F2 = Foam::tanh(sqr(arg2));

        nuti[celli] = a1*kc/max(a1*omc, F2*S);
    }
}


void Foam::ibWallModel::correctNut(volScalarField& nut) const
{
    if (stressMode_ && wallStrain_)
    {
        wallStrainNut(nut);
    }

    scalarField& nuti = nut.primitiveFieldRef();
    const labelList& ibc = ib_.ibCells();

    // Velocity mode: the wall-model value. Stress mode: the turbulence
    // model's value, blended towards the wall-consistent value in cells
    // close to the wall
    forAll(ibc, i)
    {
        if (!stressMode_)
        {
            nuti[ibc[i]] = nutIB_[i];
        }
        else if (nutWeight_[i] > 0)
        {
            nuti[ibc[i]] =
                (1 - nutWeight_[i])*nuti[ibc[i]] + nutWeight_[i]*nutIB_[i];
        }
    }

    for (const label celli : ib_.solidCells())
    {
        nuti[celli] = 0;
    }

    nut.correctBoundaryConditions();
}


void Foam::ibWallModel::initialiseFreshCells
(
    const labelList& freshCells,
    volVectorField& U
) const
{
    const fvMesh& mesh = ib_.mesh();
    const labelList& ibc = ib_.ibCells();

    const bitSet isFresh(mesh.nCells(), freshCells);

    // Wall velocity of the cells that are not reconstructed
    const pointField& C = mesh.C().primitiveField();
    const vectorField Uw0
    (
        ib_.solidVelocity(freshCells, pointField(C, freshCells))
    );

    const mapDistribute& map =
        centredCPCCellToCellStencilObject::New(mesh).map();

    const List<vector> flatU
    (
        flatField<vector>(map, U.primitiveField(), vector::zero)
    );

    const scalarField nu(nuIB());
    const scalar ks = law_->ks();

    // Wall velocity at the wall points of the IB cells
    const vectorField Uw(ib_.solidVelocity(ibc, wallPoints_));

    vectorField& Ui = U.primitiveFieldRef();

    forAll(freshCells, i)
    {
        Ui[freshCells[i]] = Uw0[i];
    }

    label nReconstructed = 0;

    forAll(ibc, i)
    {
        if (!isFresh.test(ibc[i]))
        {
            continue;
        }

        scalar uTau = 0;
        vector tHat(Zero);
        scalar unRatio = 0;

        if (sampleNeighbours(i, flatU, Uw[i], nu[i], uTau, tHat, unRatio))
        {
            const scalar yPlus = yIB_[i]*uTau/nu[i];
            const scalar uPlus = law_->uPlus(yPlus, ks*uTau/nu[i]);

            Ui[ibc[i]] = Uw[i] + uTau*uPlus*tHat + unRatio*nIB_[i];
            ++nReconstructed;
        }
    }

    U.correctBoundaryConditions();

    Info<< "    Fresh cells: "
        << returnReduce(nReconstructed, sumOp<label>())
        << " IB cells from the wall law, "
        << returnReduce(freshCells.size() - nReconstructed, sumOp<label>())
        << " at the wall velocity" << endl;
}


void Foam::ibWallModel::reportMaxShear
(
    const volVectorField& U,
    const UList<vector>& flatU,
    const vectorField& Uw,
    const scalarField& nu
) const
{
    const labelList& ibc = ib_.ibCells();
    const pointField& C = ib_.mesh().C().primitiveField();
    const vectorField& Ui = U.primitiveField();

    // IB cell with the largest wall shear stress
    label iMax = -1;
    scalar tauMax = -1;
    forAll(ibc, i)
    {
        if (mag(tauW_[i]) > tauMax)
        {
            tauMax = mag(tauW_[i]);
            iMax = i;
        }
    }
    const scalar gTauMax = returnReduce(tauMax, maxOp<scalar>());

    // Cell with the largest velocity
    label cMax = -1;
    scalar uMax = -1;
    forAll(Ui, celli)
    {
        if (mag(Ui[celli]) > uMax)
        {
            uMax = mag(Ui[celli]);
            cMax = celli;
        }
    }
    const scalar gUMax = returnReduce(uMax, maxOp<scalar>());

    if (iMax >= 0 && tauMax == gTauMax)
    {
        const label celli = ibc[iMax];
        const vector& ni = nIB_[iMax];
        const vector Ur(Ui[celli] - Uw[iMax]);
        const vector ut(Ur - (Ur & ni)*ni);

        scalar uTauN = 0;
        vector tHatN(Zero);
        scalar unRatio = 0;
        const bool sampled =
            sampleNeighbours(iMax, flatU, Uw[iMax], nu[iMax], uTauN, tHatN,
            unRatio);

        Pout<< "IB max shear: tau " << gTauMax << " at " << C[celli]
            << " y_c/h " << yIB_[iMax]/hIB_[iMax]
            << " h " << hIB_[iMax]
            << " |U| " << mag(Ui[celli]) << " |u_t| " << mag(ut)
            << " u_tau " << uTau_[iMax]
            << " u_tau,cell " << law_->uTau(mag(ut), yIB_[iMax], nu[iMax])
            << " u_tau,nbr " << (sampled ? uTauN : -1)
            << " nNbr " << nbrIdx_[iMax].size()
            << " nutIB " << nutIB_[iMax] << " w " << nutWeight_[iMax]
            << " fresh " << ib_.freshCells().found(celli)
            << endl;
    }

    if (cMax >= 0 && uMax == gUMax)
    {
        Pout<< "IB max |U|: " << gUMax << " at " << C[cMax]
            << " status " << ib_.status()[cMax]
            << " distance " << ib_.distance()[cMax]
            << endl;
    }
}


Foam::tmp<Foam::scalarField> Foam::ibWallModel::epsilonIB() const
{
    return Cmu_*kIB_*omegaIB_;
}


// ************************************************************************* //
