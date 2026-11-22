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

Application
    Test-bedRegionMesh

Description
    Build the static plan-view bed mesh (bedRegionMesh) of a flat patch,
    check it, and compare the geometry of the finite-area meshes built on
    it and directly on the patch of the flow mesh: face areas, edge
    lengths and delta coefficients must agree to round-off, edge for edge.

Usage
    \verbatim
    Test-bedRegionMesh -patch bottom [-parallel]
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"
#include "bedRegionMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- faMesh definition with one faPatch per non-coupled neighbouring patch
dictionary faDefinition(const polyMesh& mesh, const word& patchName)
{
    dictionary definition;
    definition.add("polyMeshPatches", wordList(1, patchName));

    dictionary boundary;
    for (const polyPatch& pp : mesh.boundaryMesh())
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

    return definition;
}


//- Largest relative difference of two lists of equal length
scalar maxRelDiff(const scalarField& a, const scalarField& b)
{
    scalar d = 0;
    forAll(a, i)
    {
        d = max(d, mag(a[i] - b[i])/max(mag(a[i]), VSMALL));
    }
    return returnReduce(d, maxOp<scalar>());
}

} // End namespace Foam


int main(int argc, char *argv[])
{
    argList::addOption("patch", "name", "flat bed patch (default bottom)");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word patchName(args.getOrDefault<word>("patch", "bottom"));

    autoPtr<polyMesh> bedMesh
    (
        bedRegionMesh::New(mesh, patchName, "bedRegion")
    );

    Info<< nl << "checkMesh of the bed region mesh:" << endl;
    const bool failed = bedMesh->checkMesh(true);
    Info<< (failed ? "    FAILED" : "    OK") << nl << endl;

    Info<< "Patches of the bed region mesh:" << nl;
    label nProcFaces = 0;
    for (const polyPatch& pp : bedMesh->boundaryMesh())
    {
        // Processor patches differ between processors: no reductions
        if (pp.coupled())
        {
            nProcFaces += pp.size();
            continue;
        }
        Info<< "    " << pp.name() << " " << pp.type() << " "
            << returnReduce(pp.size(), sumOp<label>()) << endl;
    }
    Info<< "    processor faces "
        << returnReduce(nProcFaces, sumOp<label>()) << endl;

    Info<< "    shared points "
        << bedMesh->globalData().nGlobalPoints() << " (bed region), "
        << mesh.globalData().nGlobalPoints() << " (flow mesh)" << endl;
    Info<< "faMesh on the flow patch ..." << endl;
    faMesh aFlow(mesh, faDefinition(mesh, patchName));
    Info<< "faMesh on the bed region ..." << endl;
    faMesh aBed(*bedMesh, faDefinition(*bedMesh, patchName));
    Info<< "faMeshes done" << endl;

    Info<< nl << "faMesh faces " << returnReduce(aFlow.nFaces(), sumOp<label>())
        << " / " << returnReduce(aBed.nFaces(), sumOp<label>())
        << ", internal edges "
        << returnReduce(aFlow.nInternalEdges(), sumOp<label>()) << " / "
        << returnReduce(aBed.nInternalEdges(), sumOp<label>()) << endl;

    bool same =
        aFlow.nFaces() == aBed.nFaces()
     && aFlow.nInternalEdges() == aBed.nInternalEdges();
    reduce(same, andOp<bool>());

    if (!same)
    {
        FatalErrorInFunction << "Topology differs" << exit(FatalError);
    }

    Info<< "max relative difference (flow patch vs bed region):" << nl
        << "    face area      " << maxRelDiff(aFlow.S().field(), aBed.S().field()) << nl
        << "    |Le| internal  "
        << maxRelDiff(aFlow.magLe().primitiveField(), aBed.magLe().primitiveField()) << nl
        << "    deltaCoeffs    "
        << maxRelDiff
           (
               aFlow.deltaCoeffs().primitiveField(),
               aBed.deltaCoeffs().primitiveField()
           ) << nl
        << "    face centre x  "
        << maxRelDiff
           (
               aFlow.areaCentres().primitiveField().component(vector::X),
               aBed.areaCentres().primitiveField().component(vector::X)
           ) << endl;

    // Boundary edges per faPatch, matched by name
    Info<< "boundary edges (flow / bed region), deltaCoeffs difference:" << nl;
    forAll(aFlow.boundary(), pi)
    {
        const faPatch& fp = aFlow.boundary()[pi];

        if (fp.coupled())
        {
            // Processor faPatches differ between processors: no reductions
            continue;
        }

        const label bi = aBed.boundary().findPatchID(fp.name());
        const label nFlow = returnReduce(fp.size(), sumOp<label>());

        if (bi < 0)
        {
            Info<< "    " << fp.name() << " " << nFlow << " / missing" << endl;
            continue;
        }
        const faPatch& bp = aBed.boundary()[bi];
        scalar d = -1;
        if (fp.size() == bp.size())
        {
            d = maxRelDiff
            (
                aFlow.deltaCoeffs().boundaryField()[pi],
                aBed.deltaCoeffs().boundaryField()[bi]
            );
        }
        else
        {
            d = returnReduce(scalar(GREAT), maxOp<scalar>());
        }
        Info<< "    " << fp.name() << " " << nFlow << " / "
            << returnReduce(bp.size(), sumOp<label>()) << "  " << d << endl;
    }

    label nProcFlow = 0, nProcBed = 0;
    for (const faPatch& fp : aFlow.boundary()) if (fp.coupled()) nProcFlow += fp.size();
    for (const faPatch& fp : aBed.boundary()) if (fp.coupled()) nProcBed += fp.size();
    Info<< "    processor edges " << returnReduce(nProcFlow, sumOp<label>())
        << " / " << returnReduce(nProcBed, sumOp<label>()) << endl;

    Info<< nl << "End" << nl << endl;
    return 0;
}


// ************************************************************************* //
