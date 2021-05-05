 /*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
     \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
      \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
 -------------------------------------------------------------------------------
 License
     This file is part of OpenFOAM.

     OpenFOAM is free software: you can redistribute it and/or modify it
     under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
     for more details.

     You should have received a copy of the GNU General Public License
     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

 Class
     Foam::cutTriSurfaceMesh

 Description
     IOoject and searching on triSurface

     Note: when constructing from dictionary has optional parameters:
         - scale     : scaling factor.
         - tolerance : relative tolerance for doing intersections
                       (see triangle::intersection)
         - minQuality: discard triangles with low quality when getting normal

 SourceFiles
     cutTriSurfaceMesh.C

 \*---------------------------------------------------------------------------*/

#ifndef cutTriSurfaceMesh_H
#define cutTriSurfaceMesh_H

#include "treeBoundBox.H"
#include "triSurfaceSearch.H"
#include "surfaceIntersection.H"
#include "indexedOctree.H"
#include "treeDataTriSurface.H"
#include "edgeIntersections.H"
#include "treeDataEdge.H"
#include "treeDataPrimitivePatch.H"
#include "DynamicList.H"
#include "triSurface.H"
#include "meshSearch.H"
#include "polyMesh.H"
#include "SubField.H"
#include "dualPatch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cutTriSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class cutTriSurfaceMesh
:
 public triSurface
{
public:
 
    typedef treeDataPrimitivePatch<polyPatch> patchTreeData;
    typedef indexedOctree<patchTreeData> patchSearch;
    typedef PrimitivePatch<face, ::Foam::List, pointField, point> PPatch;

private:
 // Private member data

     //- include faces
     mutable labelList* includeFacesPtr_;

     //- include edges
     mutable labelList* includeEdgesPtr_;

     //- include points
     mutable labelList* includePointsPtr_;

     //- boundary faces label
     mutable labelList* boundaryFacesPtr_;
     mutable labelListList* boundaryFacesPatchIDPtr_;
     mutable List<Pair<point>>* boundaryFacesNearestPatchFaceInfoPtr_;

     //- boundary edges label
     mutable labelList* boundaryEdgesPtr_;
     mutable labelListList* boundaryEdgesPatchIDPtr_;
     mutable List<Pair<point>>* boundaryEdgesNearestPatchFaceInfoPtr_;

     //- outside points and its nearest boundarypoints or includepoints
     mutable labelList* outsidePointsPtr_;
     mutable labelListList* outsidePointsStencilPtr_;
     mutable scalarListList* outsidePointsStencilWeightPtr_;

     //- include faces
     mutable labelList* globalIncludeFacesPtr_;

     //- include edges
     mutable labelList* globalIncludeEdgesPtr_;

     //- include points
     mutable labelList* globalIncludePointsPtr_;

     //- boundary faces label
     mutable labelList* globalBoundaryFacesPtr_;
     mutable labelListList* globalBoundaryFacesPatchIDPtr_;
     mutable labelListList* globalBoundaryFacesPatchIDProcPtr_;
     mutable List<Pair<point>>* globalBoundaryFacesNearestPatchFaceInfoPtr_;

     //- boundary edges label
     mutable labelList* globalBoundaryEdgesPtr_;
     mutable labelListList* globalBoundaryEdgesPatchIDPtr_;
     mutable labelListList* globalBoundaryEdgesPatchIDProcPtr_;
     mutable List<Pair<point>>* globalBoundaryEdgesNearestPatchFaceInfoPtr_;

    //- cutTriSurfaceMesh related information
    //- dualPatch
    mutable dualPatch* dualPatchPtr_;
    
    //- Include dual faces 
    mutable labelList* includeDualFacesPtr_;
    
    //- Include dual edges 
    mutable labelList* includeDualEdgesPtr_;
    
    //- Include dual points 
    mutable labelList* includeDualPointsPtr_;    
    
    //- boundary dual edges 
    mutable labelList* boundaryDualEdgesPtr_;    
    
    //- boundary dual edge patch ID
    mutable labelListList* boundaryDualEdgesPatchIDPtr_;      

     //- Disallow default bitwise copy construct
     cutTriSurfaceMesh(const cutTriSurfaceMesh&);

     //- Disallow default bitwise assignment
     void operator=(const cutTriSurfaceMesh&);


    void makeBoundaryFacesAndEdges
    (
        const polyMesh& mesh
    ) const;
    
    void makeBoundaryFacesAndEdges
    (
        const polyMesh& mesh,
        const PtrList<triSurface>& addSurfs
    ) const;
    
    
    
    void makeGlobal
    (
        const polyMesh& mesh
    ) const;

    void makeGlobal
    (
        const polyMesh& mesh,
        const PtrList<triSurface>& addSurfs
    ) const;


    void makeOutsidePointsStencil() const;
    
    void makeDualPatch() const;    

    labelList findNearestBoundaryPoints
    (
        const label& pointID,
        const labelHashSet& includePointsSet
    ) const;

    template<class T>
    void sendToMaster
    (
        List<T>& Values
    ) const;

public:

 //- Runtime type information
 TypeName("cutTriSurfaceMesh");

 //- Return the mesh sub-directory name (usually "triSurface")
 static word meshSubDir;


 // Constructors

     //- Construct from triSurface
     cutTriSurfaceMesh
     (
         const triSurface& surf,
         const polyMesh& mesh
     );
     
      cutTriSurfaceMesh
     (
         const triSurface& surf,
         const polyMesh& mesh,
         const PtrList<triSurface>& addSurfs
     );



 //- Destructor
 virtual ~cutTriSurfaceMesh();

     //- Clear storage
     void clearOut();


 // Member Functions
     //- include faces
     const labelList& includeFaces() const
     {
         return *includeFacesPtr_;
     }

     //- include edges
     const labelList& includeEdges() const
     {
         return *includeEdgesPtr_;
     }

     //- include points
     const labelList& includePoints() const
     {
         return *includePointsPtr_;
     }

     //- boundary faces label
     const labelList& boundaryFaces() const
     {
         return *boundaryFacesPtr_;
     }
     const labelListList& boundaryFacesPatchID() const
     {
         return *boundaryFacesPatchIDPtr_;
     }
     const List<Pair<point>>& boundaryFacesNearestPatchFaceInfo() const
     {
         return *boundaryFacesNearestPatchFaceInfoPtr_;
     }

     //- boundary edges label
     const labelList& boundaryEdges() const
     {
         return *boundaryEdgesPtr_;
     }
     const labelListList& boundaryEdgesPatchID() const
     {
         return *boundaryEdgesPatchIDPtr_;
     }
     const List<Pair<point>>& boundaryEdgesNearestPatchFaceInfo() const
     {
         return *boundaryEdgesNearestPatchFaceInfoPtr_;
     }

     const labelList& outsidePoints() const
     {
         return *outsidePointsPtr_;
     }
     const labelListList& outsidePointsStencil() const
     {
         return *outsidePointsStencilPtr_;
     }
     const scalarListList& outsidePointsStencilWeight() const
     {
         return *outsidePointsStencilWeightPtr_;
     }

     const labelList& globalIncludeFaces() const
     {
         return *globalIncludeFacesPtr_;
     }

     //- include edges
     const labelList& globalIncludeEdges() const
     {
         return *globalIncludeEdgesPtr_;
     }

     //- include points
     const labelList& globalIncludePoints() const
     {
         return *globalIncludePointsPtr_;
     }

     //- boundary faces label
     const labelList& globalBoundaryFaces() const
     {
         return *globalBoundaryFacesPtr_;
     }
     const labelListList& globalBoundaryFacesPatchID() const
     {
         return *globalBoundaryFacesPatchIDPtr_;
     }
     const labelListList& globalBoundaryFacesPatchIDProc() const
     {
         return *globalBoundaryFacesPatchIDProcPtr_;
     }
     const List<Pair<point>>& globalBoundaryFacesNearestPatchFaceInfo() const
     {
         return *globalBoundaryFacesNearestPatchFaceInfoPtr_;
     }

     //- boundary edges label
     const labelList& globalBoundaryEdges() const
     {
         return *globalBoundaryEdgesPtr_;
     }
     const labelListList& globalBoundaryEdgesPatchID() const
     {
         return *globalBoundaryEdgesPatchIDPtr_;
     }
     const labelListList& globalBoundaryEdgesPatchIDProc() const
     {
         return *globalBoundaryEdgesPatchIDProcPtr_;
     }
     const List<Pair<point>>& globalBoundaryEdgesNearestPatchFaceInfo() const
     {
         return *globalBoundaryEdgesNearestPatchFaceInfoPtr_;
     }

    labelList markFaces
    (
        const triSurfaceSearch& querySurf,
        const triSurfaceSearch& queryPatchSurf
    )const;
 
 
 
    bool ifInStl
    (
        const vector& C,
        const triSurfaceSearch& querySurf
    ) const;
    
    
    
    
    
//- cutTriSurfaceMesh related information

     const dualPatch& newDualPatch() const
     {
         if(!dualPatchPtr_)
         {
            makeDualPatch();
         }
         return *dualPatchPtr_;
     }
         
     //- include faces
     const labelList& includeDualFaces() const
     {
         if(!includeDualFacesPtr_)
         {
            makeDualPatch();
         }         
         return *includeDualFacesPtr_;
     }

     //- include edges
     const labelList& includeDualEdges() const
     {
         if(!includeDualEdgesPtr_)
         {
            makeDualPatch();
         }         
         return *includeDualEdgesPtr_;
     }

     //- include points
     const labelList& includeDualPoints() const
     {
         if(!includeDualPointsPtr_)
         {
            makeDualPatch();
         }         
         return *includeDualPointsPtr_;
     }
     
     //- boundary edges label
     const labelList& boundaryDualEdges() const
     {
         if(!boundaryDualEdgesPtr_)
         {
            makeDualPatch();
         }         
         return *boundaryDualEdgesPtr_;
     }
     
     const labelListList& boundaryDualEdgesPatchID() const
     {
         if(!boundaryDualEdgesPatchIDPtr_)
         {
            makeDualPatch();
         }         
         return *boundaryDualEdgesPatchIDPtr_;
     }    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
