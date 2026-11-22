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

#include "heightFieldBed.H"
#include "fac.H"
#include "immersedBoundary.H"
#include "ibWallModel.H"
#include "addToRunTimeSelectionTable.H"
#include "faCFD.H"
#include "globalIndex.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "flipOp.H"
#include "syncTools.H"
#include "processorFaPatch.H"
#include "uniformDimensionedFields.H"
#include "gravityMeshObject.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heightFieldBed, 0);
    addToRunTimeSelectionTable(bedModel, heightFieldBed, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::heightFieldBed::makeColumnMap()
{
    const faMesh& aMesh = sediment_.aMesh();
    const vectorField& Cf = aMesh.areaCentres().primitiveField();

    // Floor-face centres of all processors, projected on z = 0. The floor
    // has O(N^(2/3)) faces; it is gathered once (static map).
    List<pointField> allCentres(Pstream::nProcs());
    allCentres[Pstream::myProcNo()] = Cf;
    allCentres[Pstream::myProcNo()].replace(vector::Z, 0);
    Pstream::allGatherList(allCentres);

    const globalIndex globalFaces(Cf.size());

    pointField centres(globalFaces.totalSize());
    forAll(allCentres, proci)
    {
        SubList<point>(centres, allCentres[proci].size(),
            globalFaces.localStart(proci)) = allCentres[proci];
    }
    allCentres.clear();

    // Nearest floor-face centre (in plan) of every cell centre
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

    const pointField& C = mesh_.C().primitiveField();
    const scalar span2 = magSqr(bb.span());

    labelList globalCol(C.size(), -1);
    forAll(C, celli)
    {
        const point p(C[celli].x(), C[celli].y(), 0);
        const pointIndexHit hit = tree.findNearest(p, span2);
        globalCol[celli] = hit.index();
    }

    // Map from the owners of the floor faces to the cells' processors
    List<Map<label>> compactMap;
    mapDistributeBase mapBase(globalFaces, globalCol, compactMap);
    columnMap_.reset(new mapDistribute(std::move(mapBase)));

    cellColumn_.transfer(globalCol);
}


template<class Type>
Foam::List<Type> Foam::heightFieldBed::columnValues
(
    const UList<Type>& faceValues
) const
{
    List<Type> values(faceValues);
    columnMap_->distribute(values);
    return values;
}


void Foam::heightFieldBed::collectShearStress(const immersedBoundary& ib)
{
    const labelList& ibc = ib.ibCells();
    const vectorField& tauW = ib.wallModel().tauW();
    const scalarField& Aw = ib.ibWallArea();
    const labelList& nearestSolid = ib.nearestSolid();
    const mapDistribute& map = *columnMap_;

    // Area-weighted contributions in compact addressing, from the IB
    // cells whose nearest surface is the bed. Cells beside an immersed
    // solid carry the solid's wall model, not the bed's: when the bed
    // passes through such a cell (a pile foot, a pipe on the bed) the bed
    // shear of its column is evaluated with the wall law from the cell
    // velocity and its height above the bed (without it the
    // columns around a pile get no shear and the bed load piles up)
    List<vector> sumTauA(map.constructSize(), Zero);
    List<scalar> sumA(map.constructSize(), Zero);

    // Two-phase solver: only wet columns (water fraction at or above the
    // threshold) drive the bed; a dry or air-covered column gets no shear
    // and is inactive for the transport. By default (wetTest column) the
    // column is wet when the IB cell or the cell above it is wet: the cut
    // IB cell fills only through its faces, and under a flow along the bed
    // (no vertical velocity, dry IB cells upstream) it keeps the water
    // fraction of the initial bed for seconds while the cell above runs
    // full (dam break 2026-11-17: alpha 0.1-0.2 under 1.5 m/s of water,
    // bed inactive until scour reclassified the cells). wetTest ibCell
    // restores the test on the IB cell alone.
    const volScalarField* alphaPtr =
        mesh_.findObject<volScalarField>("alpha.water");
    const scalar wetThreshold =
        dict_.getOrDefault<scalar>("wetThreshold", 0.5);
    const word wetTest = dict_.getOrDefault<word>("wetTest", "column");
    if (wetTest != "column" && wetTest != "ibCell")
    {
        FatalIOErrorInFunction(dict_)
            << "Unknown wetTest " << wetTest << "; use column or ibCell"
            << exit(FatalIOError);
    }
    const bool columnTest = (wetTest == "column");
    const vectorField& Cc = mesh_.C().primitiveField();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const vectorField& Sf = mesh_.faceAreas();
    auto wet = [&](const label celli) -> bool
    {
        if (!alphaPtr)
        {
            return true;
        }
        const scalarField& alpha = alphaPtr->primitiveField();
        if (alpha[celli] >= wetThreshold)
        {
            return true;
        }
        if (columnTest)
        {
            // The neighbour across the upper (near-horizontal) face; a
            // face on a processor boundary keeps the IB-cell test
            for (const label facei : mesh_.cells()[celli])
            {
                if (facei >= mesh_.nInternalFaces())
                {
                    continue;
                }
                const label other =
                    (own[facei] == celli ? nei[facei] : own[facei]);
                if
                (
                    Cc[other].z() > Cc[celli].z()
                 && mag(Sf[facei].z()) > 0.9*mag(Sf[facei])
                 && alpha[other] >= wetThreshold
                )
                {
                    return true;
                }
            }
        }
        return false;
    };
    label nDry = 0;

    bool anyForeign = false;
    forAll(ibc, i)
    {
        if (nearestSolid[ibc[i]] != -1)
        {
            anyForeign = true;
            continue;
        }
        if (!wet(ibc[i]))
        {
            ++nDry;
            continue;
        }

        const label col = cellColumn_[ibc[i]];
        sumTauA[col] += tauW[i]*Aw[i];
        sumA[col] += Aw[i];
    }

    // Collective (the bed distance gathers the columns): all ranks decide
    reduce(anyForeign, orOp<bool>());

    if (anyForeign && mesh_.foundObject<volVectorField>("U"))
    {
        scalarField dBed;
        vectorField nBed;
        cellDistance(dBed, nBed);
        const vectorField& U =
            mesh_.lookupObject<volVectorField>("U").primitiveField();
        const scalarField nu(ib.wallModel().nuOfIBCells());
        const ibWallLaw& law = ib.wallModel().law();

        forAll(ibc, i)
        {
            const label celli = ibc[i];
            if (nearestSolid[celli] == -1 || !wet(celli))
            {
                continue;
            }
            // First layer above the bed: a cut cell or a cell whose face
            // lies on the bed (centre at half the extent); the next layer
            // is at 1.5 extents
            const scalar y = dBed[celli];
            if (y <= 0 || y > 0.75*ib.cellExtent(celli, nBed[celli]))
            {
                continue;
            }
            const vector& n = nBed[celli];
            const vector Ut(U[celli] - (U[celli] & n)*n);
            const scalar magUt = mag(Ut);
            if (magUt < VSMALL)
            {
                continue;
            }
            const scalar uTau = law.uTau(magUt, y, nu[i]);
            const label col = cellColumn_[celli];
            // Weight: the horizontal extent of the cell (its bed area)
            const scalar A = mesh_.V()[celli]/max(2*y, SMALL);
            sumTauA[col] += sqr(uTau)*(Ut/magUt)*A;
            sumA[col] += A;
        }
    }

    // Sum the contributions on the owners of the floor faces
    const label nFaces = sediment_.aMesh().nFaces();

    mapDistributeBase::distribute
    (
        UPstream::commsTypes::nonBlocking,
        List<labelPair>(),
        nFaces,
        map.constructMap(),
        map.constructHasFlip(),
        map.subMap(),
        map.subHasFlip(),
        sumTauA,
        vector::zero,
        plusEqOp<vector>(),
        flipOp()
    );

    mapDistributeBase::distribute
    (
        UPstream::commsTypes::nonBlocking,
        List<labelPair>(),
        nFaces,
        map.constructMap(),
        map.constructHasFlip(),
        map.subMap(),
        map.subHasFlip(),
        sumA,
        scalar(0),
        plusEqOp<scalar>(),
        flipOp()
    );

    areaScalarField& activeField = sediment_.active();
    scalarField& active = activeField.primitiveFieldRef();
    forAll(active, facei)
    {
        active[facei] = (sumA[facei] > VSMALL ? 1 : 0);
    }
    activeField.correctBoundaryConditions();
    if (alphaPtr)
    {
        reduce(nDry, sumOp<label>());
        Info<< "Bed: " << nDry << " dry IB cells (alpha.water < "
            << wetThreshold << ") carry no shear" << endl;
    }

    areaVectorField& tauField = sediment_.tau();
    vectorField& tau = tauField.primitiveFieldRef();
    forAll(tau, facei)
    {
        tau[facei] =
        (
            sumA[facei] > VSMALL ? sumTauA[facei]/sumA[facei] : vector::zero
        );
        tau[facei].z() = 0;
    }
    tauField.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heightFieldBed::heightFieldBed
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    bedModel(mesh, dict),
    sediment_(mesh, mesh, dict, dict.get<word>("floorPatch")),
    updateDistance_(dict.getOrDefault<scalar>("updateDistance", 5e-4)),
    floorPatch_(mesh.boundaryMesh()[dict.get<word>("floorPatch")]),
    prescribedShear_
    (
        PatchFunction1<vector>::NewIfPresent
        (
            floorPatch_,
            "prescribedShear",
            dict
        )
    ),
    prescribedElevation_
    (
        PatchFunction1<scalar>::NewIfPresent
        (
            floorPatch_,
            "prescribedElevation",
            dict
        )
    ),
    prescribedLiquefiedDepth_
    (
        PatchFunction1<scalar>::NewIfPresent
        (
            floorPatch_,
            "prescribedLiquefiedDepth",
            dict
        )
    ),
    soil_(nullptr),
    rhoW_(1000),
    fsReferenceLevel_(0),
    fsReferenceLevelSet_(false),
    fsAtmosphericPressure_(0),
    reportSolids_(dict.getOrDefault<wordList>("reportSolidForces", wordList())),
    pipesDict_(dict.subOrEmptyDict("liquefiedPipes")),
    pipes_(),
    pipeSolids_(),
    pipeLines_(),
    pipeLengths_(),
    structure_(nullptr),
    structureSolids_(),
    structureOrigins_(),
    structureRs_(),
    structureHydroForce_(true),
    structureRefSet_(false),
    structureDPrev_(Zero),
    structureThetaPrev_(Zero),
    structureV_(Zero),
    structureOmega_(Zero),
    maskDirty_(true),
    spansDirty_(true)
{
    makeColumnMap();

    // Two-phase reference state: bed entries, overridden by the soil's
    rhoW_ = dict.getOrDefault<scalar>("waterDensity", 1000);
    fsReferenceLevelSet_ = dict.found("referenceLevel");
    fsReferenceLevel_ = dict.getOrDefault<scalar>("referenceLevel", 0);
    fsAtmosphericPressure_ = dict.getOrDefault<scalar>("atmosphericPressure", 0);

    if (dict.found("soil"))
    {
        const dictionary& soilDict = dict.subDict("soil");
        rhoW_ = soilDict.getOrDefault<scalar>("waterDensity", rhoW_);
        if (soilDict.found("referenceLevel"))
        {
            fsReferenceLevelSet_ = true;
            fsReferenceLevel_ = soilDict.get<scalar>("referenceLevel");
        }
        fsAtmosphericPressure_ = soilDict.getOrDefault<scalar>
        (
            "atmosphericPressure", fsAtmosphericPressure_
        );
        soil_.reset(new soilRegion(mesh, floorPatch_.name(), soilDict));
    }

    if (pipesDict_.size() && !soil_)
    {
        FatalIOErrorInFunction(dict)
            << "liquefiedPipes needs a soil region" << exit(FatalIOError);
    }

    if (dict.found("structure"))
    {
        if (!soil_)
        {
            FatalIOErrorInFunction(dict)
                << "structure needs a soil region" << exit(FatalIOError);
        }
        const dictionary& sd = dict.subDict("structure");
        structure_.reset(new rigidStructure("structure", sd));
        structureSolids_ = sd.getOrDefault<wordList>("solids", wordList());
        structureHydroForce_ = sd.getOrDefault<bool>("hydrodynamicForce", true);
        soil_->attachStructure(structure_.get());

        // Restart: the structure continues from its written state
        if (structure_->readState(mesh.time()))
        {
            structureDPrev_ = structure_->displacement();
            structureThetaPrev_ = structure_->rotation();
        }
    }

    classified();

    // The flow mesh must extend below the bed: a column whose bed is at or
    // below the floor patch has no solid cell, hence no wall-modelled cell
    // and no bed shear (it stays flat without any message otherwise)
    {
        const labelList& faceLabels = sediment_.aMesh().faceLabels();
        const vectorField Cf(floorPatch_.faceCentres());
        const scalarField& zb = sediment_.zb().primitiveField();
        label nLow = 0;
        scalar minClearance = GREAT;
        forAll(faceLabels, i)
        {
            const scalar c = zb[i] - Cf[faceLabels[i] - floorPatch_.start()].z();
            minClearance = min(minClearance, c);
            nLow += (c <= SMALL);
        }
        reduce(nLow, sumOp<label>());
        reduce(minClearance, minOp<scalar>());
        if (nLow)
        {
            WarningInFunction
                << "the bed is at or below the floor patch "
                << floorPatch_.name() << " in " << nLow
                << " columns (minimum clearance " << minClearance
                << " m): those columns have no cells below the bed, get no"
                << " bed shear and cannot scour. Extend the flow mesh below"
                << " the bed." << endl;
        }
        else
        {
            Info<< "    Height-field bed: minimum clearance of the bed above"
                << " the floor patch " << minClearance << " m" << endl;
        }
    }

    if (prescribedShear_)
    {
        Info<< "    Height-field bed: prescribed bed shear stress ("
            << prescribedShear_->type() << ") replaces the wall shear" << endl;
    }
    if (prescribedElevation_)
    {
        Info<< "    Height-field bed: prescribed elevation ("
            << prescribedElevation_->type() << ") replaces the Exner step"
            << endl;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::heightFieldBed::floorValues
(
    const PatchFunction1<Type>& func
) const
{
    // Patch-face values at the current time, mapped to the finite-area
    // faces (faceLabels are the poly faces of the floor patch)
    const Field<Type> onPatch(func.value(mesh_.time().value()));
    const labelList& faceLabels = sediment_.aMesh().faceLabels();

    auto tvalues = tmp<Field<Type>>::New(faceLabels.size());
    Field<Type>& values = tvalues.ref();

    forAll(faceLabels, i)
    {
        values[i] = onPatch[faceLabels[i] - floorPatch_.start()];
    }

    return tvalues;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heightFieldBed::cellDistance
(
    scalarField& distance,
    vectorField& normal
) const
{
    const areaScalarField& zb = sediment_.zb();
    const areaVectorField gradZb(fac::grad(zb));

    const List<scalar> zCol(columnValues<scalar>(zb.primitiveField()));
    const List<vector> gCol(columnValues<vector>(gradZb.primitiveField()));

    const pointField& C = mesh_.C().primitiveField();

    distance.resize(C.size());
    normal.resize(C.size());

    forAll(C, celli)
    {
        const label col = cellColumn_[celli];
        const vector g(gCol[col].x(), gCol[col].y(), 0);
        const scalar norm = Foam::sqrt(1 + magSqr(g));

        distance[celli] = (C[celli].z() - zCol[col])/norm;
        normal[celli] = vector(-g.x(), -g.y(), 1)/norm;
    }
}


Foam::tmp<Foam::scalarField> Foam::heightFieldBed::pointLevel() const
{
    return pointLevel(sediment_.zb());
}


Foam::tmp<Foam::scalarField> Foam::heightFieldBed::pointLevel
(
    const areaScalarField& zb
) const
{
    const faMesh& aMesh = sediment_.aMesh();
    const areaVectorField gradZb(fac::grad(zb));

    const List<scalar> zCol(columnValues<scalar>(zb.primitiveField()));
    const List<vector> gCol(columnValues<vector>(gradZb.primitiveField()));
    const List<vector> cCol
    (
        columnValues<vector>(aMesh.areaCentres().primitiveField())
    );

    const pointField& P = mesh_.points();
    const labelListList& pointCells = mesh_.pointCells();

    scalarField sum(P.size(), Zero);
    scalarField count(P.size(), Zero);

    forAll(P, pointi)
    {
        for (const label celli : pointCells[pointi])
        {
            const label col = cellColumn_[celli];
            const vector& g = gCol[col];
            const vector& c = cCol[col];

            const scalar zBed =
                zCol[col]
              + g.x()*(P[pointi].x() - c.x())
              + g.y()*(P[pointi].y() - c.y());

            sum[pointi] += P[pointi].z() - zBed;
            count[pointi] += 1;
        }
    }

    syncTools::syncPointList(mesh_, sum, plusEqOp<scalar>(), scalar(0));
    syncTools::syncPointList(mesh_, count, plusEqOp<scalar>(), scalar(0));

    return sum/max(count, scalar(1));
}


bool Foam::heightFieldBed::evolve(const immersedBoundary& ib, const scalar dt)
{
    if (spansDirty_)
    {
        updateSolidSpans(ib);
    }

    collectShearStress(ib);

    if (prescribedShear_)
    {
        // Testing: prescribed shear on the active faces
        const vectorField values(floorValues(*prescribedShear_));
        areaVectorField& tauField = sediment_.tau();
        vectorField& tau = tauField.primitiveFieldRef();
        const scalarField& active = sediment_.active().primitiveField();

        forAll(tau, facei)
        {
            tau[facei] = (active[facei] > 0.5 ? values[facei] : vector::zero);
            tau[facei].z() = 0;
        }
        tauField.correctBoundaryConditions();
    }

    if (prescribedLiquefiedDepth_)
    {
        // Testing: prescribed liquefied depth (overridden below by the
        // soil region when there is one)
        sediment_.setLiquefiedDepth(floorValues(*prescribedLiquefiedDepth_));
    }

    // Hydrodynamic force and moment (about the origin) on the named solids
    if (reportSolids_.size())
    {
        const scalarField pIB(ibPressure(ib));
        const scalarField rhoIB(ibDensity(ib));
        forAll(ib.solids(), solidi)
        {
            const word& name = ib.solids()[solidi].name();
            if (!reportSolids_.found(name))
            {
                continue;
            }
            const vector F(ib.solidForce(solidi, pIB, rhoIB));
            const vector M(ib.solidMoment(solidi, pIB, rhoIB, point::zero));
            Info<< "Solid force: " << name << " F " << F << " N, M " << M
                << " N m" << endl;
            if (Pstream::master())
            {
                const fileName dir
                (
                    mesh_.time().globalPath()/"postProcessing"/"solidForces"
                );
                mkDir(dir);
                const fileName file(dir/name + ".dat");
                const bool exists = isFile(file);
                OFstream os(file, IOstreamOption(), IOstreamOption::APPEND_ATE);
                if (!exists)
                {
                    os  << "# t [s]  F_x F_y F_z [N]  M_x M_y M_z [N m]"
                        << " about the origin" << nl;
                }
                os  << mesh_.time().value() << ' '
                    << F.x() << ' ' << F.y() << ' ' << F.z() << ' '
                    << M.x() << ' ' << M.y() << ' ' << M.z() << nl;
            }
        }
    }

    // Soil region: pore pressure from the bed pressure, seepage gradient
    // back to the sediment core
    if (soil_)
    {
        if (maskDirty_)
        {
            updateStructureMask(ib);
        }

        // Pressure [Pa] and density of the IB cells for the solid forces
        const scalarField pIB(ibPressure(ib));
        const scalarField rhoIB(ibDensity(ib));

        // Structure: the hydrodynamic force and moment of this step on
        // its solids enter its balance in the soil solve
        if (structure_)
        {
            vector F(Zero);
            vector M(Zero);
            if (structureHydroForce_)
            {
                forAll(ib.solids(), solidi)
                {
                    if (structureSolids_.found(ib.solids()[solidi].name()))
                    {
                        F += ib.solidForce(solidi, pIB, rhoIB);
                        M += ib.solidMoment
                        (
                            solidi, pIB, rhoIB, structure_->refPoint()
                        );
                    }
                }
            }
            structure_->setAppliedLoads(F, M);
        }

        soil_->evolve(collectPressure(ib));

        // Structure: velocity and angular velocity of this step from the
        // new state (applied to the solids by immersedBoundary::evolveBed)
        if (structure_)
        {
            const vector d(structure_->displacement());
            const vector theta(structure_->rotation());
            structureV_ = (d - structureDPrev_)/dt;
            structureOmega_ = (theta - structureThetaPrev_)/dt;
            structureDPrev_ = d;
            structureThetaPrev_ = theta;
            structure_->writeHistory
            (
                mesh_.time().globalPath(), mesh_.time().value()
            );
            if (structure_->free())
            {
                structure_->report();
            }
        }
        sediment_.setSeepageGradient(soil_->topGradient());
        if (!prescribedLiquefiedDepth_)
        {
            sediment_.setLiquefiedDepth(soil_->liquefiedDepth());
        }
        sediment_.setCapacityIndex(soil_->capacityIndex());

        if (sediment_.seepage().horizontalForce())
        {
            // Horizontal gradient of the pore pressure just below the bed
            const faMesh& aMesh = sediment_.aMesh();
            areaScalarField pDepth
            (
                IOobject
                (
                    "bedPoreP",
                    mesh_.time().timeName(),
                    aMesh.thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                aMesh,
                dimensionedScalar(dimPressure, Zero),
                "zeroGradient"
            );
            pDepth.primitiveFieldRef() = soil_->pressureAtDepth();
            pDepth.correctBoundaryConditions();

            const areaVectorField gradP(fac::grad(pDepth));
            sediment_.setSeepageForce(gradP.primitiveField(), soil_->gammaS());
        }

        if (mesh_.time().writeTime())
        {
            soil_->write();
            if (structure_)
            {
                structure_->writeState(mesh_.time());
            }
        }

        // Liquefied pipes at their columns (built from the solid geometry
        // at the first step)
        if (pipesDict_.size() && pipes_.empty())
        {
            buildPipes(ib);
        }
        else
        {
            // Burial from the current solid position and bed elevation
            // (scour or deposition at the pipe)
            forAll(pipes_, i)
            {
                scalar zBed, y, zLow, zHigh;
                probePipe
                (
                    ib, pipeSolids_[i], pipeLines_[i].x(),
                    pipeLines_[i].y(), zBed, y, zLow, zHigh
                );
                pipes_[i].setBurial(zBed - zLow);

                // Vertical force of the flow per unit length
                if (pipeLengths_[i] > 0)
                {
                    label solidi = -1;
                    forAll(ib.solids(), j)
                    {
                        if (ib.solids()[j].name() == pipeSolids_[i])
                        {
                            solidi = j;
                        }
                    }
                    const vector F(ib.solidForce(solidi, pIB, rhoIB));
                    pipes_[i].setExternalForce(F.z()/pipeLengths_[i]);
                }
            }
        }
        for (liquefiedPipe& p : pipes_)
        {
            p.update(dt);
        }
    }

    // Flow spin-up: the bed is frozen before morphStartTime
    const bool morph = sediment_.morphActive();
    if (morph)
    {
        if (prescribedElevation_)
        {
            // Testing: the bed follows the prescribed elevation
            sediment_.setElevation(floorValues(*prescribedElevation_));
        }
        else
        {
            sediment_.advance(dt);
        }
    }

    // After the advance: the surface then holds the shear stress and the
    // bed load of this step and the new z_b (as the z_b field files)
    if (mesh_.time().writeTime())
    {
        sediment_.writeSurface();
        sediment_.writeMorphologyState();
    }

    if (!morph)
    {
        return false;
    }

    const scalar maxChange =
        returnReduce
        (
            max(mag(sediment_.zb().primitiveField() - zbClassified_)()),
            maxOp<scalar>()
        );

    return (maxChange > updateDistance_);
}


void Foam::heightFieldBed::classified()
{
    zbClassified_ = sediment_.zb().primitiveField();

    // Moving solids: the spans of the sediment core follow them
    if (pipes_.size())
    {
        spansDirty_ = true;
    }

    if (soil_)
    {
        soil_->setBed(zbAtFloorPoints());
        maskDirty_ = true;
    }
}


Foam::tmp<Foam::scalarField> Foam::heightFieldBed::erodibleIbCells
(
    const immersedBoundary& ib
) const
{
    // 1 - fixed of the floor face of each IB cell's column
    const scalarField erodibleFaces
    (
        scalar(1) - sediment_.fixed().primitiveField()
    );
    const List<scalar> cols(columnValues<scalar>(erodibleFaces));

    const labelList& ibc = ib.ibCells();
    auto tres = tmp<scalarField>::New(ibc.size(), scalar(1));
    scalarField& res = tres.ref();

    forAll(ibc, i)
    {
        const label col = cellColumn_[ibc[i]];
        if (col >= 0 && col < cols.size())
        {
            res[i] = cols[col];
        }
    }

    return tres;
}


void Foam::heightFieldBed::setSuspendedExchange
(
    const immersedBoundary& ib,
    const scalarField& netDeposition
)
{
    // Area-weighted average of the IB-cell exchange on the floor faces,
    // through the same column map as the shear stress
    const labelList& ibc = ib.ibCells();
    const scalarField& Aw = ib.ibWallArea();
    const mapDistribute& map = *columnMap_;

    if (netDeposition.size() != ibc.size())
    {
        FatalErrorInFunction
            << "Expected " << ibc.size() << " IB-cell values, got "
            << netDeposition.size() << exit(FatalError);
    }

    const labelList& nearestSolid = ib.nearestSolid();

    List<scalar> sumXA(map.constructSize(), Zero);
    List<scalar> sumA(map.constructSize(), Zero);

    forAll(ibc, i)
    {
        if (nearestSolid[ibc[i]] != -1)
        {
            continue;
        }

        const label col = cellColumn_[ibc[i]];
        sumXA[col] += netDeposition[i]*Aw[i];
        sumA[col] += Aw[i];
    }

    const label nFaces = sediment_.aMesh().nFaces();

    for (List<scalar>* fld : {&sumXA, &sumA})
    {
        mapDistributeBase::distribute
        (
            UPstream::commsTypes::nonBlocking,
            List<labelPair>(),
            nFaces,
            map.constructMap(),
            map.constructHasFlip(),
            map.subMap(),
            map.subHasFlip(),
            *fld,
            scalar(0),
            plusEqOp<scalar>(),
            flipOp()
        );
    }

    scalarField faceValues(nFaces, Zero);
    forAll(faceValues, facei)
    {
        if (sumA[facei] > VSMALL)
        {
            faceValues[facei] = sumXA[facei]/sumA[facei];
        }
    }

    sediment_.setSuspendedExchange(faceValues);
}


void Foam::heightFieldBed::updateSolidSpans(const immersedBoundary& ib)
{
    spansDirty_ = false;

    const faMesh& aMesh = sediment_.aMesh();
    const label nFaces = aMesh.nFaces();

    if (ib.solids().empty())
    {
        sediment_.setSolidSpans(List<List<Pair<scalar>>>());
        return;
    }

    // Trace the column line of every floor face from the top of the mesh
    // down to the floor, stepping by the distance to the nearest solid
    // (sphere tracing) with a floor on the step; a sign change is
    // bisected to locate the surface
    const vectorField& Cf = aMesh.areaCentres().primitiveField();
    const boundBox& bb = mesh_.bounds();
    const scalar zTop = bb.max().z();
    const scalar height = bb.max().z() - bb.min().z();
    const scalar dzMin = 1e-4*height;
    const scalar tol = 1e-8*height;

    auto minDistance = [&](const pointField& pts, scalarField& d)
    {
        d.resize(pts.size());
        d = GREAT;
        scalarField ds;
        vectorField ns;
        forAll(ib.solids(), solidi)
        {
            ib.solids()[solidi].distance(pts, ds, ns);
            d = min(d, ds);
        }
    };

    List<DynamicList<Pair<scalar>>> spans(nFaces);
    scalarField zTopSpan(nFaces, GREAT);

    pointField pts(nFaces);
    forAll(pts, facei)
    {
        pts[facei] = point(Cf[facei].x(), Cf[facei].y(), zTop);
    }

    scalarField d;
    minDistance(pts, d);
    boolList done(nFaces, false);
    label nActive = nFaces;

    for (label iter = 0; iter < 100000 && nActive; ++iter)
    {
        pointField ptsNew(pts);
        forAll(pts, facei)
        {
            if (!done[facei])
            {
                ptsNew[facei].z() -= max(mag(d[facei]), dzMin);
            }
        }

        scalarField dNew;
        minDistance(ptsNew, dNew);

        // Faces whose column line crossed a surface in this step
        DynamicList<label> crossed(nFaces/8);
        forAll(pts, facei)
        {
            if (!done[facei] && (dNew[facei] < 0) != (d[facei] < 0))
            {
                crossed.append(facei);
            }
        }

        if (crossed.size())
        {
            pointField lo(crossed.size());
            pointField hi(crossed.size());
            scalarField dLo(crossed.size());
            forAll(crossed, i)
            {
                hi[i] = pts[crossed[i]];
                lo[i] = ptsNew[crossed[i]];
                dLo[i] = dNew[crossed[i]];
            }

            pointField mid(crossed.size());
            scalarField dMid;
            for (label k = 0; k < 40; ++k)
            {
                bool converged = true;
                forAll(mid, i)
                {
                    mid[i] = 0.5*(lo[i] + hi[i]);
                    if (hi[i].z() - lo[i].z() > tol)
                    {
                        converged = false;
                    }
                }
                if (converged)
                {
                    break;
                }

                minDistance(mid, dMid);
                forAll(mid, i)
                {
                    if ((dMid[i] < 0) == (dLo[i] < 0))
                    {
                        lo[i] = mid[i];
                        dLo[i] = dMid[i];
                    }
                    else
                    {
                        hi[i] = mid[i];
                    }
                }
            }

            forAll(crossed, i)
            {
                const label facei = crossed[i];
                const scalar zs = 0.5*(lo[i].z() + hi[i].z());

                if (dNew[facei] < 0)
                {
                    // Entered a solid: its top
                    zTopSpan[facei] = zs;
                }
                else
                {
                    // Left a solid: its bottom (the span is clipped at
                    // the top of the mesh if the trace started inside)
                    spans[facei].append
                    (
                        Pair<scalar>(zs, min(zTopSpan[facei], zTop))
                    );
                    zTopSpan[facei] = GREAT;
                }
            }
        }

        pts = ptsNew;
        d = dNew;

        forAll(pts, facei)
        {
            if (!done[facei] && pts[facei].z() <= Cf[facei].z())
            {
                if (d[facei] < 0)
                {
                    // Solid down to the floor
                    spans[facei].append
                    (
                        Pair<scalar>(Cf[facei].z(), zTopSpan[facei])
                    );
                }
                done[facei] = true;
                --nActive;
            }
        }
    }

    if (nActive)
    {
        FatalErrorInFunction
            << "Tracing the solid spans did not finish for " << nActive
            << " floor faces" << exit(FatalError);
    }

    // Bottom-up order and totals
    List<List<Pair<scalar>>> spanList(nFaces);
    const scalarField& S = aMesh.S().field();
    const scalarField& zb = sediment_.zb().primitiveField();
    label nSpanFaces = 0;
    label maxSpans = 0;
    scalar volume = 0;
    scalar volumeBelowBed = 0;

    forAll(spans, facei)
    {
        if (spans[facei].size())
        {
            ++nSpanFaces;
            maxSpans = max(maxSpans, spans[facei].size());
            spanList[facei].setSize(spans[facei].size());
            forAll(spans[facei], i)
            {
                // Traced top-down: reverse
                spanList[facei][i] = spans[facei][spans[facei].size() - 1 - i];
                const scalar zBot = spanList[facei][i].first();
                const scalar zTopi = spanList[facei][i].second();
                volume += S[facei]*(zTopi - zBot);
                volumeBelowBed +=
                    S[facei]*max(min(zb[facei], zTopi) - zBot, scalar(0));
            }
        }
    }

    reduce(nSpanFaces, sumOp<label>());
    reduce(maxSpans, maxOp<label>());
    reduce(volume, sumOp<scalar>());
    reduce(volumeBelowBed, sumOp<scalar>());

    Info<< "    Bed: immersed solids over " << nSpanFaces
        << " floor faces (up to " << maxSpans << " spans per column), "
        << "solid volume in the columns " << volume
        << " m3, of which below z_b " << volumeBelowBed << " m3" << endl;

    sediment_.setSolidSpans(std::move(spanList));
}


void Foam::heightFieldBed::updateStructureMask(const immersedBoundary& ib)
{
    maskDirty_ = false;

    if (ib.solids().empty())
    {
        return;
    }

    const pointField& C = soil_->mesh().cellCentres();
    scalarField mask(C.size(), Zero);

    // Solids moved by liquefied pipes are not part of the soil mechanics
    // (the reduced model carries the pipe-soil interaction)
    wordHashSet pipeSolids;
    for (const entry& e : pipesDict_)
    {
        if (e.isDict())
        {
            pipeSolids.insert(e.dict().get<word>("solid"));
        }
    }

    // The solids of the structure are masked with their initial geometry
    // (small-displacement soil mechanics: the mask stays, the masked cells
    // carry the rigid displacement); recorded at the first call
    if (structure_ && structureOrigins_.empty())
    {
        if (structureSolids_.empty())
        {
            forAll(ib.solids(), solidi)
            {
                if (!pipeSolids.found(ib.solids()[solidi].name()))
                {
                    structureSolids_.append(ib.solids()[solidi].name());
                }
            }
        }
        structureOrigins_.setSize(structureSolids_.size());
        structureRs_.setSize(structureSolids_.size());
        forAll(structureSolids_, i)
        {
            bool found = false;
            forAll(ib.solids(), solidi)
            {
                if (ib.solids()[solidi].name() == structureSolids_[i])
                {
                    structureOrigins_[i] = ib.solids()[solidi].origin();
                    structureRs_[i] = ib.solids()[solidi].R();
                    found = true;
                }
            }
            if (!found)
            {
                FatalErrorInFunction
                    << "structure solid " << structureSolids_[i]
                    << " is not an immersed solid" << exit(FatalError);
            }
        }
        Info<< "    Structure: solids " << structureSolids_ << endl;
    }

    forAll(ib.solids(), solidi)
    {
        const word& solidName = ib.solids()[solidi].name();
        if (pipeSolids.found(solidName))
        {
            continue;
        }

        scalarField d;
        vectorField n;
        const label si = structureSolids_.find(solidName);
        if (structure_ && si >= 0)
        {
            ib.solids()[solidi].distance
            (
                C, d, n, structureOrigins_[si], structureRs_[si]
            );
        }
        else
        {
            ib.solids()[solidi].distance(C, d, n);
        }
        forAll(mask, celli)
        {
            if (d[celli] < 0)
            {
                mask[celli] = 1;
            }
        }
    }

    soil_->setStructureMask(mask);

    // Reference point of the structure: the volume centroid of its masked
    // cells when the dictionary gave none
    if (structure_ && !structureRefSet_)
    {
        const scalarField& Vc = soil_->mesh().V();
        scalar V = 0;
        vector xV(Zero);
        forAll(mask, celli)
        {
            if (mask[celli] > 0.5)
            {
                V += Vc[celli];
                xV += Vc[celli]*C[celli];
            }
        }
        reduce(V, sumOp<scalar>());
        reduce(xV, sumOp<vector>());
        if (V > VSMALL)
        {
            structure_->setRefPoint(xV/V);
        }
        structureRefSet_ = true;
        Info<< "    Structure: embedded volume " << V
            << " m3, reference point " << structure_->refPoint() << endl;
    }
}


Foam::tmp<Foam::areaScalarField> Foam::heightFieldBed::soilSurface() const
{
    const areaScalarField& zb = sediment_.zb();
    auto tzs = tmp<areaScalarField>::New("soilSurface", zb);
    areaScalarField& zs = tzs.ref();
    scalarField& z = zs.primitiveFieldRef();

    // Covered: z_b inside a solid span of the column (the span list is
    // empty on every processor without solids; all processors take part
    // in the reductions below)
    const List<List<Pair<scalar>>>& spans = sediment_.solidSpans();
    scalarField f(z.size(), 1);
    label nCovered = 0;
    forAll(spans, facei)
    {
        for (const Pair<scalar>& sp : spans[facei])
        {
            if (z[facei] >= sp.first() - SMALL && z[facei] <= sp.second() + SMALL)
            {
                f[facei] = 0;
                ++nCovered;
                break;
            }
        }
    }
    if (returnReduce(nCovered, sumOp<label>()) == 0)
    {
        return tzs;
    }

    // Jacobi passes: an unfilled column takes the average of its filled
    // neighbours across the internal edges and the processor patches
    // (values exchanged with non-blocking buffers: a blocking exchange per
    // patch deadlocks when the patch order differs between processors)
    const faMesh& aMesh = sediment_.aMesh();
    const labelUList& own = aMesh.edgeOwner();
    const labelUList& nei = aMesh.edgeNeighbour();
    const label nInt = aMesh.nInternalEdges();
    const faBoundaryMesh& patches = aMesh.boundary();

    for (label pass = 0; pass < 200; ++pass)
    {
        scalarField sum(z.size(), Zero);
        scalarField cnt(z.size(), Zero);

        for (label e = 0; e < nInt; ++e)
        {
            const label o = own[e];
            const label n = nei[e];
            if (f[o] < 0.5 && f[n] > 0.5)
            {
                sum[o] += z[n];
                cnt[o] += 1;
            }
            if (f[n] < 0.5 && f[o] > 0.5)
            {
                sum[n] += z[o];
                cnt[n] += 1;
            }
        }

        if (Pstream::parRun())
        {
            PstreamBuffers pBufs;
            forAll(patches, patchi)
            {
                const auto* ppp = isA<processorFaPatch>(patches[patchi]);
                if (!ppp)
                {
                    continue;
                }
                const labelUList& ef = ppp->edgeFaces();
                scalarField zp(ef.size());
                scalarField fp(ef.size());
                forAll(ef, i)
                {
                    zp[i] = z[ef[i]];
                    fp[i] = f[ef[i]];
                }
                UOPstream toNbr(ppp->neighbProcNo(), pBufs);
                toNbr << zp << fp;
            }
            pBufs.finishedSends();
            forAll(patches, patchi)
            {
                const auto* ppp = isA<processorFaPatch>(patches[patchi]);
                if (!ppp)
                {
                    continue;
                }
                scalarField zNbr;
                scalarField fNbr;
                {
                    UIPstream fromNbr(ppp->neighbProcNo(), pBufs);
                    fromNbr >> zNbr >> fNbr;
                }
                const labelUList& ef = ppp->edgeFaces();
                forAll(ef, i)
                {
                    const label o = ef[i];
                    if (f[o] < 0.5 && fNbr[i] > 0.5)
                    {
                        sum[o] += zNbr[i];
                        cnt[o] += 1;
                    }
                }
            }
        }

        label nNew = 0;
        forAll(z, facei)
        {
            if (f[facei] < 0.5 && cnt[facei] > 0)
            {
                z[facei] = sum[facei]/cnt[facei];
                f[facei] = 1;
                ++nNew;
            }
        }
        if (returnReduce(nNew, sumOp<label>()) == 0)
        {
            break;
        }
    }
    zs.correctBoundaryConditions();

    return tzs;
}


Foam::tmp<Foam::scalarField> Foam::heightFieldBed::zbAtFloorPoints() const
{
    // phi = z - z_b at the mesh points of the soil surface; the floor
    // points lie at z_floor
    const scalarField phi(pointLevel(soilSurface()));
    const labelList& meshPoints = floorPatch_.meshPoints();
    const pointField& P = mesh_.points();

    auto tzb = tmp<scalarField>::New(meshPoints.size());
    scalarField& zb = tzb.ref();

    forAll(meshPoints, i)
    {
        zb[i] = P[meshPoints[i]].z() - phi[meshPoints[i]];
    }

    return tzb;
}


Foam::tmp<Foam::scalarField> Foam::heightFieldBed::collectPressure
(
    const immersedBoundary& ib
) const
{
    const labelList& ibc = ib.ibCells();
    const scalarField& Aw = ib.ibWallArea();
    const mapDistribute& map = *columnMap_;

    // Two-phase solver: p_rgh and the mixture density of the IB cell give
    // the total pressure at the bed of its column, whatever the phase of
    // the cell, p_b = p_rgh + rho |g| (hRef - z_b) (hydrostatic over the
    // half cell); the soil takes the excess over the still-water
    // hydrostatic pressure at the bed, p_b - p_atm - rho_w |g| (z_ref - z_b),
    // which is rho_w g eta under a wave and -rho_w g (z_ref - z_b) on a bed
    // the water has left. Single-phase solver: p is the kinematic excess
    // pressure and the bed pressure is rho_w p.
    const bool twoPhase = mesh_.foundObject<volScalarField>("p_rgh");
    const volScalarField& p =
        mesh_.lookupObject<volScalarField>(twoPhase ? "p_rgh" : "p");

    List<scalar> sumPA(map.constructSize(), Zero);
    List<scalar> sumRA(map.constructSize(), Zero);
    List<scalar> sumA(map.constructSize(), Zero);

    if (twoPhase)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
        forAll(ibc, i)
        {
            const label col = cellColumn_[ibc[i]];
            sumPA[col] += p[ibc[i]]*Aw[i];
            sumRA[col] += rho[ibc[i]]*Aw[i];
            sumA[col] += Aw[i];
        }
    }
    else
    {
        forAll(ibc, i)
        {
            const label col = cellColumn_[ibc[i]];
            sumPA[col] += p[ibc[i]]*Aw[i];
            sumA[col] += Aw[i];
        }
    }

    const label nFaces = sediment_.aMesh().nFaces();

    for (List<scalar>* fld : {&sumPA, &sumRA, &sumA})
    {
        mapDistributeBase::distribute
        (
            UPstream::commsTypes::nonBlocking,
            List<labelPair>(),
            nFaces,
            map.constructMap(),
            map.constructHasFlip(),
            map.subMap(),
            map.subHasFlip(),
            *fld,
            scalar(0),
            plusEqOp<scalar>(),
            flipOp()
        );
    }

    auto tpBed = tmp<scalarField>::New(nFaces, Zero);
    scalarField& pBed = tpBed.ref();

    if (twoPhase)
    {
        const scalar gMag =
            mag(meshObjects::gravity::New(mesh_.time()).value());
        const scalar hRef =
            mesh_.lookupObject<uniformDimensionedScalarField>("hRef").value();
        const scalar zRef = fsReferenceLevelSet_ ? fsReferenceLevel_ : hRef;
        const scalarField& zb = sediment_.zb().primitiveField();

        forAll(pBed, facei)
        {
            if (sumA[facei] > VSMALL)
            {
                const scalar pRgh = sumPA[facei]/sumA[facei];
                const scalar rho = sumRA[facei]/sumA[facei];
                pBed[facei] =
                    pRgh + rho*gMag*(hRef - zb[facei])
                  - fsAtmosphericPressure_ - rhoW_*gMag*(zRef - zb[facei]);
            }
        }
    }
    else
    {
        forAll(pBed, facei)
        {
            pBed[facei] =
                (sumA[facei] > VSMALL ? rhoW_*sumPA[facei]/sumA[facei] : 0);
        }
    }

    return tpBed;
}


Foam::tmp<Foam::scalarField> Foam::heightFieldBed::ibPressure
(
    const immersedBoundary& ib
) const
{
    const labelList& ibc = ib.ibCells();
    auto tpIB = tmp<scalarField>::New(ibc.size(), Zero);
    scalarField& pIB = tpIB.ref();

    if (!mesh_.foundObject<volScalarField>("p"))
    {
        return tpIB;
    }
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");

    if (mesh_.foundObject<volScalarField>("p_rgh"))
    {
        // Two-phase solver: the total pressure at the wall point of the
        // IB cell, extrapolated from the cell centre with the hydrostatic
        // gradient of the cell's phase (the cell-centre value alone
        // overstates a hydrostatic load by rho g d, 10-20 % of the
        // buoyancy of a cylinder at 8 cells per diameter), less the
        // still-water hydrostatic pressure at that point
        const vector g(meshObjects::gravity::New(mesh_.time()).value());
        const scalar gMag = mag(g);
        const scalar zRef =
        (
            fsReferenceLevelSet_
          ? fsReferenceLevel_
          : mesh_.lookupObject<uniformDimensionedScalarField>("hRef").value()
        );
        const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
        const scalarField& d = ib.distance().primitiveField();
        const vectorField& n = ib.normal().primitiveField();
        const pointField& C = mesh_.cellCentres();

        forAll(ibc, i)
        {
            const label celli = ibc[i];
            const vector xw(C[celli] - d[celli]*n[celli]);
            pIB[i] =
                p[celli] + rho[celli]*(g & (xw - C[celli]))
              - fsAtmosphericPressure_
              - rhoW_*gMag*max(zRef - xw.z(), scalar(0));
        }
    }
    else
    {
        forAll(ibc, i)
        {
            pIB[i] = rhoW_*p[ibc[i]];
        }
    }

    return tpIB;
}


Foam::tmp<Foam::scalarField> Foam::heightFieldBed::ibDensity
(
    const immersedBoundary& ib
) const
{
    const labelList& ibc = ib.ibCells();
    auto trhoIB = tmp<scalarField>::New(ibc.size(), rhoW_);

    if (mesh_.foundObject<volScalarField>("p_rgh"))
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
        scalarField& rhoIB = trhoIB.ref();
        forAll(ibc, i)
        {
            rhoIB[i] = rho[ibc[i]];
        }
    }

    return trhoIB;
}


void Foam::heightFieldBed::probePipe
(
    const immersedBoundary& ib,
    const word& solidName,
    const scalar x,
    const scalar L,
    scalar& zBed,
    scalar& y,
    scalar& zLow,
    scalar& zHigh
) const
{
    label solidi = -1;
    forAll(ib.solids(), i)
    {
        if (ib.solids()[i].name() == solidName)
        {
            solidi = i;
        }
    }
    if (solidi < 0)
    {
        FatalErrorInFunction
            << "No immersed solid " << solidName << exit(FatalError);
    }

    // Bed elevation and floor-face centre of the column nearest x
    const vectorField Cf(floorPatch_.faceCentres());
    scalar dLoc = GREAT;
    label iLoc = -1;
    forAll(Cf, i)
    {
        if (mag(Cf[i].x() - x) < dLoc)
        {
            dLoc = mag(Cf[i].x() - x);
            iLoc = i;
        }
    }
    const scalar dMin = returnReduce(dLoc, minOp<scalar>());
    const label owner = returnReduce
    (
        (iLoc >= 0 && dLoc <= dMin) ? Pstream::myProcNo() : Pstream::nProcs(),
        minOp<label>()
    );
    zBed = 0;
    y = 0;
    if (Pstream::myProcNo() == owner)
    {
        zBed = soil_->bedElevation(iLoc);
        y = Cf[iLoc].y();
    }
    reduce(zBed, sumOp<scalar>());
    reduce(y, sumOp<scalar>());

    // The lowest and highest point of the solid on the column line, from
    // the signed distance of probes below and above it
    pointField probes(2);
    probes[0] = point(x, y, zBed - L);
    probes[1] = point(x, y, zBed + L);
    scalarField d;
    vectorField n;
    ib.solids()[solidi].distance(probes, d, n);
    zLow = probes[0].z() + d[0];
    zHigh = probes[1].z() - d[1];
}


void Foam::heightFieldBed::buildPipes(const immersedBoundary& ib)
{
    for (const entry& e : pipesDict_)
    {
        if (!e.isDict())
        {
            continue;
        }
        const dictionary& pd = e.dict();
        const word solidName(pd.get<word>("solid"));
        const scalar x = pd.get<scalar>("x");
        const scalar L = pd.getOrDefault<scalar>("probeDistance", 1);

        scalar zBed, y, zLow, zHigh;
        probePipe(ib, solidName, x, L, zBed, y, zLow, zHigh);
        const scalar D = zHigh - zLow;
        const scalar e0 = max(zBed - zLow, scalar(0));

        Info<< "    Liquefied pipe " << e.keyword() << " moves solid "
            << solidName << ": from the geometry D " << D << " m, bottom "
            << zLow << " m, bed " << zBed << " m, initial burial " << e0
            << " m" << endl;

        pipes_.append(new liquefiedPipe(*soil_, e.keyword(), pd, D, e0));
        pipeSolids_.append(solidName);
        pipeLines_.append(vector2D(x, L));

        const boundBox& bb = mesh_.bounds();
        pipeLengths_.append
        (
            pd.getOrDefault<bool>("hydrodynamicForce", true)
          ? pd.getOrDefault<scalar>("length", bb.max().y() - bb.min().y())
          : scalar(0)
        );
    }
}


void Foam::heightFieldBed::solidVelocities
(
    HashTable<vector>& velocities
) const
{
    // liquefiedPipe velocities are positive downward
    forAll(pipes_, i)
    {
        velocities.set(pipeSolids_[i], vector(0, 0, -pipes_[i].velocity()));
    }

    // Structure: the velocity of each solid's reference point (origin)
    // from the rigid motion about the structure reference point
    if (structure_)
    {
        forAll(structureSolids_, i)
        {
            velocities.set
            (
                structureSolids_[i],
                structureV_
              + (structureOmega_ ^ (structureOrigins_[i] - structure_->refPoint()))
            );
        }
    }
}


void Foam::heightFieldBed::solidAngularVelocities
(
    HashTable<vector>& omegas
) const
{
    if (structure_)
    {
        forAll(structureSolids_, i)
        {
            omegas.set(structureSolids_[i], structureOmega_);
        }
    }
}


// ************************************************************************* //
