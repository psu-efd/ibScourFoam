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

#include "bedRegionMesh.H"
#include "processorPolyPatch.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "Pair.H"
#include <type_traits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::bedRegionMesh::topPatchName("bedRegionTop");


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- Quad above patch edge (a, b) between the point levels k and k + 1,
//- starting at a and oriented away from point c
static face sideFace
(
    const label a,
    const label b,
    const label k,
    const label nPoints,
    const pointField& points,
    const point& c
)
{
    face f(4);
    f[0] = a + k*nPoints;
    f[1] = b + k*nPoints;
    f[2] = b + (k + 1)*nPoints;
    f[3] = a + (k + 1)*nPoints;

    if ((f.areaNormal(points) & (f.centre(points) - c)) < 0)
    {
        f.flip();
    }

    return f;
}


//- Patch face lifted to point level k, oriented along +z (up) or -z
static face levelFace
(
    const face& localFace,
    const label k,
    const label nPoints,
    const pointField& points,
    const bool up
)
{
    face f(localFace);
    forAll(f, fp)
    {
        f[fp] += k*nPoints;
    }

    if ((f.areaNormal(points).z() > 0) != up)
    {
        f.flip();
    }

    return f;
}


//- Add the patches of a column mesh to a polyMesh (or fvMesh)
template<class MeshType>
static void addColumnPatches
(
    MeshType& mesh,
    const bedRegionMesh::columns& col
)
{
    const polyBoundaryMesh& bm = mesh.boundaryMesh();
    List<polyPatch*> patches(col.patchNames.size());

    forAll(patches, i)
    {
        if (col.patchNbrProcs[i] < 0)
        {
            patches[i] = new polyPatch
            (
                col.patchNames[i],
                col.patchSizes[i],
                col.patchStarts[i],
                i,
                bm,
                polyPatch::typeName
            );
        }
        else
        {
            patches[i] = new processorPolyPatch
            (
                col.patchSizes[i],
                col.patchStarts[i],
                i,
                bm,
                Pstream::myProcNo(),
                col.patchNbrProcs[i]
            );
        }
    }

    if constexpr (std::is_same<MeshType, fvMesh>::value)
    {
        mesh.addFvPatches(patches, true);
    }
    else
    {
        mesh.addPatches(patches, true);
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

Foam::bedRegionMesh::columns Foam::bedRegionMesh::makeColumns
(
    const fvMesh& mesh,
    const word& patchName,
    const scalarList& zLevels,
    const word& bottomPatchName,
    const word& topPatchName
)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const label patchi = pbm.findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch " << patchName << " not found"
            << exit(FatalError);
    }

    const label nLayers = zLevels.size() - 1;

    if (nLayers < 1)
    {
        FatalErrorInFunction
            << "At least two point levels are needed, got " << zLevels
            << exit(FatalError);
    }

    const polyPatch& pp = pbm[patchi];
    const faceList& localFaces = pp.localFaces();
    const edgeList& edges = pp.edges();
    const labelListList& edgeFaces = pp.edgeFaces();
    const label nInternalEdges = pp.nInternalEdges();
    const label nPts = pp.nPoints();
    const label nFcs = pp.size();

    columns col;
    col.nColumns = nFcs;
    col.nLayers = nLayers;
    col.nPointsPerLevel = nPts;

    // Points: the plan view at each level
    pointField& points = col.points;
    points.setSize((nLayers + 1)*nPts);

    forAll(pp.localPoints(), pointi)
    {
        const point& p = pp.localPoints()[pointi];
        for (label k = 0; k <= nLayers; ++k)
        {
            points[pointi + k*nPts] = point(p.x(), p.y(), zLevels[k]);
        }
    }

    // Cell centres (plan view) and the vertical position of the layers
    pointField planCentres(nFcs);
    forAll(localFaces, facei)
    {
        planCentres[facei] = localFaces[facei].centre(points);
        planCentres[facei].z() = 0;
    }

    auto cellCentre = [&](const label facei, const label k) -> point
    {
        point c(planCentres[facei]);
        c.z() = 0.5*(zLevels[k] + zLevels[k + 1]);
        return c;
    };

    // Global numbering of the mesh points (shared points are equal on all
    // processors): identifies edges across processors and orders the faces
    // of the processor patches
    labelList meshPointToGlobal;
    {
        labelList uniquePoints;
        mesh.globalData().mergePoints(meshPointToGlobal, uniquePoints);
    }
    const labelList pointToGlobal(meshPointToGlobal, pp.meshPoints());

    // Side faces of the boundary edges: neighbouring patch of the flow mesh.
    // A physical patch (e.g. a structure) takes precedence over a processor
    // patch that touches the same edge: the neighbouring processor then has
    // no bed face on that edge.
    const labelList& meshEdges = pp.meshEdges();
    const labelListList& meshEdgeFaces = mesh.edgeFaces();

    labelList edgePatch(edges.size(), -1);

    for (label edgei = nInternalEdges; edgei < edges.size(); ++edgei)
    {
        for (const label facei : meshEdgeFaces[meshEdges[edgei]])
        {
            if (!mesh.isInternalFace(facei))
            {
                const label pi = pbm.whichPatch(facei);

                if (pi == patchi)
                {
                    continue;
                }

                if (!isA<processorPolyPatch>(pbm[pi]))
                {
                    edgePatch[edgei] = pi;
                    break;
                }

                if (edgePatch[edgei] < 0)
                {
                    edgePatch[edgei] = pi;
                }
            }
        }

        if (edgePatch[edgei] < 0)
        {
            FatalErrorInFunction
                << "No neighbouring patch for edge " << edgei
                << " of patch " << patchName << exit(FatalError);
        }
    }

    // An edge on a processor boundary is a processor face only if the
    // neighbouring processor also has a face of the patch on it. Otherwise
    // (e.g. the bed meets a structure whose face lies on the other
    // processor) it takes the physical patch that the neighbour sees on
    // the edge.
    if (Pstream::parRun())
    {
        // Candidate edges per neighbouring processor
        Map<DynamicList<label>> procEdges;
        Map<label> procPatch;

        forAll(pbm, pi)
        {
            const auto* ppp = isA<processorPolyPatch>(pbm[pi]);
            if (ppp)
            {
                procEdges.insert(ppp->neighbProcNo(), DynamicList<label>());
                procPatch.insert(ppp->neighbProcNo(), pi);
            }
        }

        for (label edgei = nInternalEdges; edgei < edges.size(); ++edgei)
        {
            const auto* ppp = isA<processorPolyPatch>(pbm[edgePatch[edgei]]);
            if (ppp)
            {
                procEdges[ppp->neighbProcNo()].append(edgei);
            }
        }

        // Send the global point pairs of the candidates
        PstreamBuffers pBufs;

        forAllConstIters(procEdges, iter)
        {
            List<Pair<label>> keys(iter.val().size());
            forAll(keys, i)
            {
                const edge& e = edges[iter.val()[i]];
                keys[i] =
                    Pair<label>
                    (
                        pointToGlobal[e.start()], pointToGlobal[e.end()]
                    );
            }

            UOPstream toNbr(iter.key(), pBufs);
            toNbr << keys;
        }

        pBufs.finishedSends();

        // Answer: -2 the neighbour has a patch face on the edge, >= 0 the
        // non-processor patch it sees there, -1 neither
        const labelListList& pointEdges = mesh.pointEdges();
        const edgeList& allEdges = mesh.edges();

        PstreamBuffers replyBufs;

        forAllConstIters(procEdges, iter)
        {
            List<Pair<label>> keys;
            {
                UIPstream fromNbr(iter.key(), pBufs);
                fromNbr >> keys;
            }

            // Global -> local numbers of the points of the processor patch
            const polyPatch& ppp = pbm[procPatch[iter.key()]];
            Map<label> globalToLocal(2*ppp.nPoints());
            for (const label pointi : ppp.meshPoints())
            {
                globalToLocal.insert(meshPointToGlobal[pointi], pointi);
            }

            labelList answers(keys.size(), -1);

            forAll(keys, i)
            {
                const auto p0 = globalToLocal.cfind(keys[i].first());
                const auto p1 = globalToLocal.cfind(keys[i].second());

                if (!p0.good() || !p1.good())
                {
                    continue;
                }

                label meshEdge = -1;
                for (const label ei : pointEdges[p0.val()])
                {
                    if (allEdges[ei].otherVertex(p0.val()) == p1.val())
                    {
                        meshEdge = ei;
                        break;
                    }
                }

                if (meshEdge < 0)
                {
                    continue;
                }

                for (const label facei : meshEdgeFaces[meshEdge])
                {
                    if (mesh.isInternalFace(facei))
                    {
                        continue;
                    }

                    const label pi = pbm.whichPatch(facei);

                    if (pi == patchi)
                    {
                        answers[i] = -2;
                        break;
                    }

                    if (!isA<processorPolyPatch>(pbm[pi]))
                    {
                        answers[i] = pi;
                    }
                }
            }

            UOPstream toNbr(iter.key(), replyBufs);
            toNbr << answers;
        }

        replyBufs.finishedSends();

        forAllConstIters(procEdges, iter)
        {
            labelList answers;
            {
                UIPstream fromNbr(iter.key(), replyBufs);
                fromNbr >> answers;
            }

            forAll(answers, i)
            {
                const label edgei = iter.val()[i];

                if (answers[i] >= 0)
                {
                    edgePatch[edgei] = answers[i];
                }
                else if (answers[i] == -1)
                {
                    FatalErrorInFunction
                        << "Edge " << edges[edgei] << " of patch "
                        << patchName << " on the boundary with processor "
                        << iter.key() << " has no patch face or physical"
                        << " patch there" << exit(FatalError);
                }
            }
        }
    }

    // Lateral neighbours of each column with a higher face number, sorted
    // (upper-triangular order of the internal faces)
    List<DynamicList<Pair<label>>> upperNbrs(nFcs);   // (nbr face, edge)

    for (label edgei = 0; edgei < nInternalEdges; ++edgei)
    {
        const labelList& ef = edgeFaces[edgei];
        const label lo = min(ef[0], ef[1]);
        const label hi = max(ef[0], ef[1]);
        upperNbrs[lo].append(Pair<label>(hi, edgei));
    }

    forAll(upperNbrs, facei)
    {
        Foam::sort(upperNbrs[facei]);
    }

    // Faces, owner, neighbour. Internal faces: for each cell in increasing
    // order, its lateral faces (sorted by neighbour) then the face to the
    // cell above.
    const label nInternalFaces = nLayers*nInternalEdges + (nLayers - 1)*nFcs;

    const label nBoundaryEdges = edges.size() - nInternalEdges;

    DynamicList<face> faces
    (
        nInternalFaces + 2*nFcs + nLayers*nBoundaryEdges
    );
    DynamicList<label> owner(faces.capacity());
    DynamicList<label> neighbour(nInternalFaces);

    for (label k = 0; k < nLayers; ++k)
    {
        for (label facei = 0; facei < nFcs; ++facei)
        {
            const label own = k*nFcs + facei;
            const point c(cellCentre(facei, k));

            for (const Pair<label>& nbr : upperNbrs[facei])
            {
                const edge& e = edges[nbr.second()];
                faces.append
                (
                    sideFace(e.start(), e.end(), k, nPts, points, c)
                );
                owner.append(own);
                neighbour.append(k*nFcs + nbr.first());
            }

            if (k < nLayers - 1)
            {
                faces.append
                (
                    levelFace(localFaces[facei], k + 1, nPts, points, true)
                );
                owner.append(own);
                neighbour.append((k + 1)*nFcs + facei);
            }
        }
    }

    // Patch sizes and starts, in the order: bottom, top, the non-processor
    // patches of the flow mesh (except the extruded patch), processor patches
    wordList& patchNames = col.patchNames;
    labelList& patchSizes = col.patchSizes;
    labelList& patchStarts = col.patchStarts;
    labelList& patchNbrProcs = col.patchNbrProcs;

    DynamicList<word> names;
    DynamicList<label> sizes;
    DynamicList<label> starts;
    DynamicList<label> nbrProcs;

    // Bottom patch: level 0, normals along -z (the cell lies above)
    names.append(bottomPatchName);
    starts.append(faces.size());
    forAll(localFaces, facei)
    {
        faces.append(levelFace(localFaces[facei], 0, nPts, points, false));
        owner.append(facei);
    }
    sizes.append(nFcs);
    nbrProcs.append(-1);

    // Top patch: level nLayers, normals along +z
    names.append(topPatchName);
    starts.append(faces.size());
    forAll(localFaces, facei)
    {
        faces.append(levelFace(localFaces[facei], nLayers, nPts, points, true));
        owner.append((nLayers - 1)*nFcs + facei);
    }
    sizes.append(nFcs);
    nbrProcs.append(-1);

    // Side patches, one per non-processor patch of the flow mesh (also
    // empty ones, so that all processors have the same patches)
    forAll(pbm, pi)
    {
        if (pi == patchi || isA<processorPolyPatch>(pbm[pi]))
        {
            continue;
        }

        names.append(pbm[pi].name());
        starts.append(faces.size());

        label n = 0;
        for (label edgei = nInternalEdges; edgei < edges.size(); ++edgei)
        {
            if (edgePatch[edgei] == pi)
            {
                const edge& e = edges[edgei];
                const label facei = edgeFaces[edgei][0];

                for (label k = 0; k < nLayers; ++k)
                {
                    faces.append
                    (
                        sideFace
                        (
                            e.start(), e.end(), k, nPts, points,
                            cellCentre(facei, k)
                        )
                    );
                    owner.append(k*nFcs + facei);
                    ++n;
                }
            }
        }

        sizes.append(n);
        nbrProcs.append(-1);
    }

    // Processor patches, faces ordered by the global numbers of their edge
    // points (lower number first, then by layer), starting at the lower
    // point so that the faces of both sides are reverses of each other with
    // the same first point
    forAll(pbm, pi)
    {
        const auto* procPatch = isA<processorPolyPatch>(pbm[pi]);

        if (!procPatch)
        {
            continue;
        }

        DynamicList<label> procEdges;
        DynamicList<Pair<label>> keys;

        for (label edgei = nInternalEdges; edgei < edges.size(); ++edgei)
        {
            if (edgePatch[edgei] == pi)
            {
                const edge& e = edges[edgei];
                const label g0 = pointToGlobal[e.start()];
                const label g1 = pointToGlobal[e.end()];
                procEdges.append(edgei);
                keys.append(Pair<label>(min(g0, g1), max(g0, g1)));
            }
        }

        if (procEdges.empty())
        {
            continue;
        }

        names.append(word::null);
        starts.append(faces.size());

        const labelList order(sortedOrder(keys));

        for (const label i : order)
        {
            const label edgei = procEdges[i];
            const edge& e = edges[edgei];
            const bool forward =
                pointToGlobal[e.start()] < pointToGlobal[e.end()];
            const label a = (forward ? e.start() : e.end());
            const label b = (forward ? e.end() : e.start());
            const label facei = edgeFaces[edgei][0];

            for (label k = 0; k < nLayers; ++k)
            {
                faces.append
                (
                    sideFace(a, b, k, nPts, points, cellCentre(facei, k))
                );
                owner.append(k*nFcs + facei);
            }
        }

        sizes.append(nLayers*procEdges.size());
        nbrProcs.append(procPatch->neighbProcNo());
    }

    col.faces = faceList(std::move(faces));
    col.owner = labelList(std::move(owner));
    col.neighbour = labelList(std::move(neighbour));
    patchNames = wordList(std::move(names));
    patchSizes = labelList(std::move(sizes));
    patchStarts = labelList(std::move(starts));
    patchNbrProcs = labelList(std::move(nbrProcs));

    return col;
}


Foam::autoPtr<Foam::polyMesh> Foam::bedRegionMesh::New
(
    const fvMesh& mesh,
    const word& patchName,
    const word& regionName
)
{
    const polyPatch& pp = mesh.boundaryMesh()[patchName];

    // Layer thickness of the order of the face size (geometry only)
    const scalar h =
        Foam::sqrt
        (
            returnReduce(sum(mag(pp.faceAreas())), sumOp<scalar>())
           /max(returnReduce(pp.size(), sumOp<label>()), 1)
        );

    columns col
    (
        makeColumns
        (
            mesh, patchName, scalarList({scalar(0), h}), patchName, topPatchName
        )
    );

    auto meshPtr = autoPtr<polyMesh>::New
    (
        IOobject
        (
            regionName,
            mesh.time().constant(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        std::move(col.points),
        std::move(col.faces),
        std::move(col.owner),
        std::move(col.neighbour),
        true
    );

    addColumnPatches(*meshPtr, col);

    Info<< "    Bed region mesh " << regionName << ": "
        << returnReduce(col.nColumns, sumOp<label>())
        << " bed faces from patch "
        << patchName << endl;

    return meshPtr;
}


Foam::autoPtr<Foam::fvMesh> Foam::bedRegionMesh::NewFvMesh
(
    const fvMesh& mesh,
    const word& patchName,
    const word& regionName,
    const scalarList& zLevels,
    const word& bottomPatchName,
    const word& topPatchName
)
{
    columns col
    (
        makeColumns(mesh, patchName, zLevels, bottomPatchName, topPatchName)
    );

    auto meshPtr = autoPtr<fvMesh>::New
    (
        IOobject
        (
            regionName,
            mesh.time().constant(),
            mesh.time(),
            IOobject::READ_IF_PRESENT,  // schemes and solution of the region
            IOobject::NO_WRITE
        ),
        std::move(col.points),
        std::move(col.faces),
        std::move(col.owner),
        std::move(col.neighbour),
        true
    );

    addColumnPatches(*meshPtr, col);

    Info<< "    Column mesh " << regionName << ": "
        << returnReduce(col.nColumns, sumOp<label>()) << " columns of "
        << col.nLayers << " cells from patch " << patchName << endl;

    return meshPtr;
}


// ************************************************************************* //
