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

#include "soilRegion.H"
#include "mathematicalConstants.H"
#include "poreGenerationModel.H"
#include "bedRegionMesh.H"
#include "IOField.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalarList Foam::soilRegion::levelDepths
(
    const label nLayers,
    const scalar grading
)
{
    // Layer thickness growing geometrically from the top: eta of the level
    // j layers below the top is (1 - q^j)/(1 - q^N)
    scalarList eta(nLayers + 1);

    const scalar q =
        (nLayers > 1 ? Foam::pow(grading, 1.0/(nLayers - 1)) : 1);

    for (label j = 0; j <= nLayers; ++j)
    {
        const scalar e =
        (
            mag(q - 1) < SMALL
          ? scalar(j)/nLayers
          : (1 - Foam::pow(q, j))/(1 - Foam::pow(q, nLayers))
        );

        eta[nLayers - j] = e;
    }

    return eta;
}


Foam::scalarList Foam::soilRegion::initialLevels
(
    const scalarList& eta,
    const scalar zBottom
)
{
    // Unit-height columns above the bottom until the bed is set
    scalarList z(eta.size());

    forAll(z, k)
    {
        z[k] = zBottom + (1 - eta[k]);
    }

    return z;
}


Foam::scalar Foam::soilRegion::zTop(const label facei) const
{
    const face& f = floorPatch_.localFaces()[facei];

    scalar z = 0;
    for (const label pointi : f)
    {
        z += zBedPoints_[pointi];
    }

    return z/f.size();
}


void Foam::soilRegion::updateDiagnostics()
{
    const fvMesh& mesh = *meshPtr_;
    const scalarField& p = modelPtr_->p().primitiveField();
    const vectorField& C = mesh.cellCentres();
    const label nCol = nColumns();

    scalarField& pe = pExcess_.primitiveFieldRef();
    scalarField& s0 = sigmaV0_.primitiveFieldRef();

    depth_.resize(pe.size());

    forAll(pe, celli)
    {
        const label facei = celli % nCol;
        pe[celli] = p[celli] - pTop_[facei];
        depth_[celli] = zTop(facei) - C[celli].z();
        s0[celli] = overburden_ + gammaS_*depth_[celli];
    }

    modelPtr_->setColumnState(s0, depth_, pe);

    // With a displacement-solving model the loss of vertical effective
    // stress is the (tension-positive) stress increment sigma'_zz, which
    // includes the pore-pressure effect through equilibrium; the Level 1
    // models assume a rigid skeleton under constant total stress, where it
    // is the excess pore pressure
    if (modelPtr_->solvesDisplacement())
    {
        const symmTensorField& sig = modelPtr_->sigma().primitiveField();
        forAll(pe, celli)
        {
            pe[celli] = sig[celli].zz();
        }
    }

    pExcess_.correctBoundaryConditions();
    sigmaV0_.correctBoundaryConditions();

    modelPtr_->geostatic(s0, k0_);
}


void Foam::soilRegion::positionColumns(const scalarField& zBedPoints)
{
    fvMesh& mesh = *meshPtr_;
    const label nPts = floorPatch_.nPoints();

    zBedPoints_ = zBedPoints;

    pointField newPoints(mesh.points());

    for (label k = 0; k <= nLayers_; ++k)
    {
        for (label pointi = 0; pointi < nPts; ++pointi)
        {
            newPoints[pointi + k*nPts].z() =
                zBedPoints[pointi] - eta_[k]*(zBedPoints[pointi] - zBottom_);
        }
    }

    mesh.movePoints(newPoints);

    // The column motion is a remap of the soil coordinates, not a material
    // motion: the time derivatives use the current volumes, not V0
    mesh.moving(false);
}


void Foam::soilRegion::readBedState()
{
    IOobject io
    (
        "zBedPoints",
        flowMesh_.time().timeName(),
        "uniform"/meshPtr_->name(),
        flowMesh_.time(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (!io.typeHeaderOk<IOField<scalar>>(true))
    {
        return;
    }

    const IOField<scalar> zBed(io);
    if (zBed.size() != floorPatch_.nPoints())
    {
        FatalErrorInFunction
            << "Soil region " << meshPtr_->name() << ": " << io.objectPath()
            << " has " << zBed.size() << " bed elevations, the floor patch "
            << floorPatch_.nPoints() << " points" << exit(FatalError);
    }

    positionColumns(zBed);

    Info<< "    Soil region " << meshPtr_->name()
        << ": columns positioned on the bed of time "
        << flowMesh_.time().timeName() << endl;

    // The relaxed bed pressure of the restart time (absent in time
    // directories of earlier versions: then set from the first flow step)
    IOobject ioP
    (
        "pBedRelaxed",
        flowMesh_.time().timeName(),
        "uniform"/meshPtr_->name(),
        flowMesh_.time(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (ioP.typeHeaderOk<IOField<scalar>>(true))
    {
        const IOField<scalar> pBed(ioP);
        if (pBed.size() == pBedRelaxed_.size())
        {
            pBedRelaxed_ = pBed;
            pBedRelaxedSet_ = true;
            Info<< "    Soil region " << meshPtr_->name()
                << ": relaxed bed pressure restored" << endl;
        }
    }
}


void Foam::soilRegion::writeBedState() const
{
    IOField<scalar> zBed
    (
        IOobject
        (
            "zBedPoints",
            flowMesh_.time().timeName(),
            "uniform"/meshPtr_->name(),
            flowMesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        zBedPoints_
    );
    zBed.write();

    IOField<scalar> pBed
    (
        IOobject
        (
            "pBedRelaxed",
            flowMesh_.time().timeName(),
            "uniform"/meshPtr_->name(),
            flowMesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pBedRelaxed_
    );
    pBed.write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::soilRegion::soilRegion
(
    const fvMesh& flowMesh,
    const word& floorPatch,
    const dictionary& dict
)
:
    flowMesh_(flowMesh),
    floorPatch_(flowMesh.boundaryMesh()[floorPatch]),
    dict_(dict),
    nLayers_(dict.get<label>("nLayers")),
    zBottom_(dict.get<scalar>("bottom")),
    eta_(levelDepths(nLayers_, dict.getOrDefault<scalar>("grading", 1))),
    gammaS_(0),
    k0_(dict.getOrDefault<scalar>("k0", 0.5)),
    overburden_(dict.getOrDefault<scalar>("overburden", 0)),
    meshPtr_
    (
        bedRegionMesh::NewFvMesh
        (
            flowMesh,
            floorPatch_.name(),
            dict.getOrDefault<word>("region", "soil"),
            initialLevels(eta_, zBottom_),
            soilModel::bottomPatchName,
            soilModel::topPatchName
        )
    ),
    modelPtr_(),
    unloading_
    (
        Function1<scalar>::NewIfPresent
        (
            "unloading", dict, word::null, &flowMesh
        )
    ),
    pressureRampTime_(dict.getOrDefault<scalar>("pressureRampTime", 0)),
    waveAmp_
    (
        dict.subOrEmptyDict("surfaceWave").getOrDefault<scalar>("amplitude", 0)
    ),
    waveK_
    (
        waveAmp_ != 0
      ? constant::mathematical::twoPi
       /dict.subDict("surfaceWave").get<scalar>("wavelength")
      : 0
    ),
    waveOmega_
    (
        waveAmp_ != 0
      ? constant::mathematical::twoPi
       /dict.subDict("surfaceWave").get<scalar>("period")
      : 0
    ),
    waveDir_
    (
        normalised
        (
            dict.subOrEmptyDict("surfaceWave")
               .getOrDefault<vector>("direction", vector(1, 0, 0))
        )
    ),
    zBedPoints_(floorPatch_.nPoints(), zBottom_ + 1),
    pTop_(floorPatch_.size(), Zero),
    pExcess_
    (
        IOobject
        (
            "pExcess",
            flowMesh.time().timeName(),
            *meshPtr_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *meshPtr_,
        dimensionedScalar(dimPressure, Zero)
    ),
    sigmaV0_
    (
        IOobject
        (
            "sigmaV0",
            flowMesh.time().timeName(),
            *meshPtr_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *meshPtr_,
        dimensionedScalar(dimPressure, Zero)
    ),
    meshWritten_(false),
    equilibrate_
    (
        dict.getOrDefault<bool>("equilibrateInitial", !dict.found("p0"))
    ),
    gradientDepth_(dict.getOrDefault<scalar>("gradientDepth", 0)),
    pressureRelaxation_(dict.getOrDefault<scalar>("pressureRelaxation", 0)),
    interval_(max(dict.getOrDefault<label>("interval", 1), label(1))),
    stepsSinceSolve_(0),
    geostaticSet_(false),
    pBedRelaxed_(floorPatch_.size(), Zero),
    pBedRelaxedSet_(false)
{
    // Restart: the columns stand under the bed written at the restart
    // time before the model reads its fields on them
    readBedState();
    modelPtr_ = soilModel::New(*meshPtr_, dict);

    if (modelPtr_->restored())
    {
        // The pore pressure is the written one, not the bed pressure
        equilibrate_ = false;
    }

    const scalar s = dict.getOrDefault<scalar>("specificGravity", 2.65);
    const scalar n = dict.getOrDefault<scalar>("porosity", 0.4);

    gammaS_ = (s - 1)*(1 - n)*modelPtr_->gammaW();

    Info<< "    Soil region " << meshPtr_->name() << ": "
        << returnReduce(nColumns(), sumOp<label>()) << " columns of "
        << nLayers_ << " cells, bottom at z = " << zBottom_
        << " m, gamma' " << gammaS_ << " N/m3"
        << (unloading_ ? ", prescribed unloading" : "")
        << (waveAmp_ != 0 ? ", prescribed surface wave" : "") << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::soilRegion::setBed(const scalarField& zBedPoints)
{
    fvMesh& mesh = *meshPtr_;
    const label nPts = floorPatch_.nPoints();

    if (zBedPoints.size() != nPts)
    {
        FatalErrorInFunction
            << "Expected " << nPts << " bed elevations, got "
            << zBedPoints.size() << exit(FatalError);
    }

    // Cell elevations before the motion, for the remap
    const scalarField zOld(mesh.cellCentres().component(vector::Z));

    positionColumns(zBedPoints);

    // Remap the fields of every column onto the new cell elevations
    // (linear in z; above the old top the old top value, below the old
    // bottom the old bottom value), so that a cell does not carry the
    // state of its old depth to a new one
    const scalarField zNew(mesh.cellCentres().component(vector::Z));
    const label nCol = nColumns();
    labelList c0(zNew.size());
    labelList c1(zNew.size());
    scalarField w(zNew.size(), Zero);

    for (label facei = 0; facei < nCol; ++facei)
    {
        for (label k = 0; k < nLayers_; ++k)
        {
            const label celli = k*nCol + facei;
            const scalar z = zNew[celli];

            // Old cells of this column are ordered upwards with k
            label kl = 0;
            while (kl < nLayers_ - 1 && zOld[(kl + 1)*nCol + facei] < z)
            {
                ++kl;
            }

            c0[celli] = kl*nCol + facei;
            c1[celli] = min(kl + 1, nLayers_ - 1)*nCol + facei;
            const scalar z0 = zOld[c0[celli]];
            const scalar z1 = zOld[c1[celli]];

            if (z <= z0 || c0[celli] == c1[celli])
            {
                w[celli] = 0;
            }
            else if (z >= z1)
            {
                w[celli] = 1;
            }
            else
            {
                w[celli] = (z - z0)/max(z1 - z0, SMALL);
            }
        }
    }

    modelPtr_->remap(c0, c1, w);
}


void Foam::soilRegion::evolve(const scalarField& pBed)
{
    if (pBed.size() != nColumns())
    {
        FatalErrorInFunction
            << "Expected " << nColumns() << " bed pressures, got "
            << pBed.size() << exit(FatalError);
    }

    // Flow pressure ramped in over the start-up
    scalarField pFlow(pBed);
    if (pressureRampTime_ > 0)
    {
        pFlow *= min(flowMesh_.time().value()/pressureRampTime_, scalar(1));
    }

    if (pressureRelaxation_ > 0 && !equilibrate_ && pBedRelaxedSet_)
    {
        const scalar w =
            min(flowMesh_.time().deltaTValue()/pressureRelaxation_, scalar(1));
        pBedRelaxed_ += w*(pFlow - pBedRelaxed_);
    }
    else
    {
        // Initial value (also of a restart without a written relaxed
        // pressure): the relaxation does not start from zero
        pBedRelaxed_ = pFlow;
        pBedRelaxedSet_ = true;
    }

    pTop_ = pBedRelaxed_;

    if (unloading_)
    {
        pTop_ += unloading_->value(flowMesh_.time().value());
    }

    if (waveAmp_ != 0)
    {
        const vectorField Cf(floorPatch_.faceCentres());
        const scalar t = flowMesh_.time().value();
        forAll(pTop_, facei)
        {
            pTop_[facei] +=
                waveAmp_
               *Foam::cos(waveK_*(Cf[facei] & waveDir_) - waveOmega_*t);
        }
    }

    if (equilibrate_)
    {
        // Initial equilibrium: the pressure of each column equals its
        // bed pressure (no initial transient from the flow pressure)
        const label nCol = nColumns();
        scalarField pCells(meshPtr_->nCells());
        forAll(pCells, celli)
        {
            pCells[celli] = pTop_[celli % nCol];
        }
        modelPtr_->initialise(pCells);
        equilibrate_ = false;

        Info<< "    Soil region: pore pressure equilibrated with the bed"
            << " pressure" << endl;
    }

    if (bedTraction_.size())
    {
        modelPtr_->topTraction(bedTraction_);
    }

    if (bedFraction_.size())
    {
        modelPtr_->topDisplacement(bedFraction_, bedDisplacement_, bedTangential_);
    }

    if (++stepsSinceSolve_ >= interval_)
    {
        if (!geostaticSet_)
        {
            updateDiagnostics();
            geostaticSet_ = true;
        }

        modelPtr_->correct(pTop_, stepsSinceSolve_);
        stepsSinceSolve_ = 0;
        updateDiagnostics();
    }
}


void Foam::soilRegion::setStructureMask(const scalarField& mask)
{
    if (mask.size() != meshPtr_->nCells())
    {
        FatalErrorInFunction
            << "Expected " << meshPtr_->nCells() << " mask values, got "
            << mask.size() << exit(FatalError);
    }

    modelPtr_->setStructureMask(mask);
}


void Foam::soilRegion::setBedTraction(const vectorField& t)
{
    if (t.size() != nColumns())
    {
        FatalErrorInFunction
            << "Expected " << nColumns() << " bed tractions, got "
            << t.size() << exit(FatalError);
    }

    bedTraction_ = t;
}


void Foam::soilRegion::setBedDisplacement
(
    const scalarField& fraction,
    const vectorField& D,
    const bool tangential
)
{
    bedTangential_ = tangential;

    if (fraction.size() != nColumns() || D.size() != nColumns())
    {
        FatalErrorInFunction
            << "Expected " << nColumns() << " bed displacements, got "
            << fraction.size() << " and " << D.size() << exit(FatalError);
    }

    bedFraction_ = fraction;
    bedDisplacement_ = D;
}


Foam::tmp<Foam::scalarField> Foam::soilRegion::topGradient() const
{
    const fvMesh& mesh = *meshPtr_;
    const label nCol = nColumns();
    const scalarField& p = modelPtr_->p().primitiveField();
    const vectorField& C = mesh.cellCentres();
    const scalar gammaW = modelPtr_->gammaW();

    auto tgrad = tmp<scalarField>::New(nCol, Zero);
    scalarField& grad = tgrad.ref();

    // Columns covered by a structure: the top soil cell is inside the
    // structure (rigid and impermeable, its pressure is no pore pressure)
    // and the bed above it is not exposed to the flow: no seepage
    // feedback there
    const scalarField& sMask = structureMask();

    forAll(grad, facei)
    {
        const scalar z0 = zTop(facei);
        const label topCelli = (nLayers_ - 1)*nCol + facei;

        if (sMask.size() == p.size() && sMask[topCelli] > 0.5)
        {
            grad[facei] = 0;
            continue;
        }

        if (gradientDepth_ <= 0)
        {
            const label celli = topCelli;
            const scalar dz = z0 - C[celli].z();
            grad[facei] = (p[celli] - pTop_[facei])/(gammaW*max(dz, SMALL));
            continue;
        }

        // Pressure at depth gradientDepth, linear between the cells of
        // the column (top value at the bed face)
        const scalar zd = z0 - gradientDepth_;
        scalar zUp = z0;
        scalar pUp = pTop_[facei];
        scalar pd = pUp;

        for (label k = nLayers_ - 1; k >= 0; --k)
        {
            const label celli = k*nCol + facei;
            const scalar zc = C[celli].z();

            if (zc <= zd)
            {
                const scalar w = (zUp - zd)/max(zUp - zc, SMALL);
                pd = (1 - w)*pUp + w*p[celli];
                break;
            }

            zUp = zc;
            pUp = p[celli];
            pd = pUp;
        }

        grad[facei] = (pd - pTop_[facei])/(gammaW*gradientDepth_);
    }

    // Diagnostics: a gradient far above critical points at a column whose
    // top cell is degenerate or whose bed pressure is not a flow value
    {
        label iMax = -1;
        scalar gMaxLocal = 0;
        forAll(grad, facei)
        {
            if (mag(grad[facei]) > gMaxLocal)
            {
                gMaxLocal = mag(grad[facei]);
                iMax = facei;
            }
        }
        const scalar gMax = returnReduce(gMaxLocal, maxOp<scalar>());
        if (gMax > 10 && iMax >= 0 && gMaxLocal == gMax)
        {
            const label celli = (nLayers_ - 1)*nCol + iMax;
            Pout<< "Soil region: bed gradient " << grad[iMax]
                << " at column " << iMax << " x " << C[celli].x()
                << " y " << C[celli].y() << ": z_top " << zTop(iMax)
                << " top-cell centre z " << C[celli].z()
                << " p " << p[celli] << " pTop " << pTop_[iMax] << endl;
            // The column profile (top down): is the column draining?
            Pout<< "    column " << iMax << " profile (z, p, masked):";
            for (label k = nLayers_ - 1; k >= max(nLayers_ - 6, 0); --k)
            {
                const label ck = k*nCol + iMax;
                Pout<< " (" << C[ck].z() << ", " << p[ck] << ", "
                    << (sMask.size() == p.size() ? sMask[ck] : -1) << ")";
            }
            Pout<< endl;
        }
    }

    return tgrad;
}


Foam::tmp<Foam::scalarField> Foam::soilRegion::pressureAtDepth() const
{
    // p_top + gamma_w * i * depth: consistent with topGradient()
    const scalarField i(topGradient()());
    const label nCol = nColumns();
    const vectorField& C = meshPtr_->cellCentres();
    const scalar gammaW = modelPtr_->gammaW();

    auto tp = tmp<scalarField>::New(nCol, Zero);
    scalarField& pd = tp.ref();

    forAll(pd, facei)
    {
        const scalar depth =
        (
            gradientDepth_ > 0
          ? gradientDepth_
          : zTop(facei) - C[(nLayers_ - 1)*nCol + facei].z()
        );
        pd[facei] = pTop_[facei] + gammaW*i[facei]*depth;
    }

    return tp;
}


Foam::tmp<Foam::scalarField> Foam::soilRegion::liquefiedDepth() const
{
    const fvMesh& mesh = *meshPtr_;
    const label nCol = nColumns();
    const scalarField& pe = pExcess_.primitiveField();
    const scalarField& s0 = sigmaV0_.primitiveField();
    const vectorField& C = mesh.cellCentres();

    auto tdepth = tmp<scalarField>::New(nCol, Zero);
    scalarField& depth = tdepth.ref();

    forAll(depth, facei)
    {
        const scalar z0 = zTop(facei);
        scalar zPrev = z0;
        scalar fPrev = 0;

        const scalarField& mask = modelPtr_->cellMask();

        for (label k = nLayers_ - 1; k >= 0; --k)
        {
            const label celli = k*nCol + facei;
            // Liquefied at p_e >= sigma'_v0 within a relative tolerance
            // (a cell capped at sigma'_v0 by the pore-pressure generation
            // sits a round-off below it after the drainage of a step)
            const scalar f =
                (mask.size() && mask[celli] > 0.5)
              ? -GREAT
              : pe[celli] - s0[celli]*(1 - 1e-4);

            if (f < 0)
            {
                // Crossing between the previous (liquefied) point and this cell
                if (k == nLayers_ - 1)
                {
                    break;
                }

                depth[facei] =
                    (z0 - zPrev)
                  + fPrev*(zPrev - C[celli].z())/max(fPrev - f, SMALL);
                break;
            }

            zPrev = C[celli].z();
            fPrev = f;
            depth[facei] = z0 - zPrev;
        }
    }

    return tdepth;
}


Foam::tmp<Foam::scalarField> Foam::soilRegion::capacityIndex() const
{
    const label nCol = nColumns();
    const scalarField& pe = pExcess_.primitiveField();
    const scalarField& s0 = sigmaV0_.primitiveField();

    auto tcap = tmp<scalarField>::New(nCol, Zero);
    scalarField& cap = tcap.ref();

    forAll(cap, facei)
    {
        const scalar height = zTop(facei) - zBottom_;

        for (label k = 0; k < nLayers_; ++k)
        {
            const label celli = k*nCol + facei;
            const scalar dz = height*(eta_[k] - eta_[k + 1]);
            cap[facei] += max(s0[celli] - pe[celli], scalar(0))*dz;
        }
    }

    return tcap;
}


Foam::tmp<Foam::scalarField> Foam::soilRegion::columnZ(const label facei) const
{
    const vectorField& C = meshPtr_->cellCentres();
    const label nCol = nColumns();

    auto tz = tmp<scalarField>::New(nLayers_);
    scalarField& z = tz.ref();

    forAll(z, j)
    {
        z[j] = C[(nLayers_ - 1 - j)*nCol + facei].z();
    }

    return tz;
}


Foam::tmp<Foam::scalarField> Foam::soilRegion::columnExcess
(
    const label facei
) const
{
    const scalarField& pe = pExcess_.primitiveField();
    const label nCol = nColumns();

    auto tp = tmp<scalarField>::New(nLayers_);
    scalarField& p = tp.ref();

    forAll(p, j)
    {
        p[j] = pe[(nLayers_ - 1 - j)*nCol + facei];
    }

    return tp;
}


Foam::tmp<Foam::scalarField> Foam::soilRegion::columnPorePressure
(
    const label facei
) const
{
    const scalarField& p = modelPtr_->p().primitiveField();
    const label nCol = nColumns();

    auto tp = tmp<scalarField>::New(nLayers_);
    scalarField& pc = tp.ref();

    forAll(pc, j)
    {
        pc[j] = p[(nLayers_ - 1 - j)*nCol + facei] - pTop_[facei];
    }

    return tp;
}


void Foam::soilRegion::write() const
{
    if (!meshWritten_ && dict_.getOrDefault<bool>("writeMesh", false))
    {
        // The region mesh is built at run time: write its full topology to
        // constant/<region>/polyMesh once (a plain write puts only the
        // points into the time directory, and post-processing tools such as
        // ParaView then cannot find the region)
        fvMesh& mesh = const_cast<fvMesh&>(*meshPtr_);
        mesh.setInstance(mesh.time().constant());
        meshPtr_->write();
        meshWritten_ = true;
    }

    modelPtr_->writeGeneration();
    modelPtr_->writeState();
    writeBedState();
}


// ************************************************************************* //
