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

Class
    Foam::bedRegionMesh

Description
    Column meshes built in memory on the faces of a patch of the flow mesh.

    The faces of the patch, projected on z = 0, are extruded along z into
    columns of nLayers cells (one column per patch face). Cell
    layer*nFaces + facei is layer 'layer' (0 at the bottom) of the column
    of patch face facei. The points of level k (0 .. nLayers) are
    (x, y, zLevels[k]); the levels are uniform at construction and may be
    moved afterwards (polyMesh::movePoints). Each boundary edge of the patch
    becomes a stack of side faces in the patch of the flow mesh that meets
    the patch along that edge (inlet, outlet, a structure, ...), and edges
    on processor boundaries become processor faces, ordered by the global
    point numbering of the flow mesh (edge first, then layer) so that both
    sides match. The bottom and top faces form two patches with the given
    names (normals along -z and +z).

    Face i of the bottom (and of the top) patch is face i of the patch of
    the flow mesh, so data are exchanged without interpolation.

    Two uses:
    - the static plan-view mesh of a moving bed patch (one layer, bottom
      patch named as the bed patch, top patch 'bedRegionTop'), for the
      sediment core (bedSediment) of a body-fitted, deforming flow mesh;
    - the soil region under a height-field bed (soilRegionMesh): nLayers
      cells per column, bottom patch 'soilBottom', top patch 'soilTop' at
      the bed elevation.

SourceFiles
    bedRegionMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_bedRegionMesh_H
#define Foam_bedRegionMesh_H

#include "fvMesh.H"
#include "polyMesh.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class bedRegionMesh Declaration
\*---------------------------------------------------------------------------*/

class bedRegionMesh
{
public:

    //- Name of the extruded (top) patch of the plan-view mesh
    static const word topPatchName;

    //- Topology and geometry of a column mesh, before the patches are added
    struct columns
    {
        pointField points;
        faceList faces;
        labelList owner;
        labelList neighbour;
        wordList patchNames;
        labelList patchSizes;
        labelList patchStarts;
        labelList patchNbrProcs;    // -1: not a processor patch
        label nColumns;             // patch faces on this processor
        label nLayers;
        label nPointsPerLevel;
    };


    // Static Member Functions

        //- Build the columns of patch patchName of mesh with the point
        //- levels zLevels (size nLayers + 1, increasing)
        static columns makeColumns
        (
            const fvMesh& mesh,
            const word& patchName,
            const scalarList& zLevels,
            const word& bottomPatchName,
            const word& topPatchName
        );

        //- Build the plan-view mesh (one layer between z = 0 and a layer
        //- thickness of the order of the face size) of patch patchName
        //- of mesh, registered with the time of mesh under regionName
        static autoPtr<polyMesh> New
        (
            const fvMesh& mesh,
            const word& patchName,
            const word& regionName
        );

        //- Build a finite-volume column mesh of patch patchName of mesh,
        //- registered with the time of mesh under regionName
        static autoPtr<fvMesh> NewFvMesh
        (
            const fvMesh& mesh,
            const word& patchName,
            const word& regionName,
            const scalarList& zLevels,
            const word& bottomPatchName,
            const word& topPatchName
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
