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
     Foam::dualPatch

 Description
     dual mesh for primitivePacth

 SourceFiles
     dualPatch.C

 \*---------------------------------------------------------------------------*/

#ifndef dualPatch_H
#define dualPatch_H

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
#include "PrimitivePatchInterpolation.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class dualPatch Declaration
\*---------------------------------------------------------------------------*/


class dualPatch
{
public:
    typedef treeDataPrimitivePatch<polyPatch> patchTreeData;
    typedef indexedOctree<patchTreeData> patchSearch;
    typedef PrimitivePatch<face, ::Foam::List, pointField, point> PPatch;
    typedef PrimitivePatchInterpolation< PrimitivePatch<face, ::Foam::List, pointField, point> >  PPatchInterpolation;

private:
 // Private member data
     
     //- old ppatch
     mutable PPatch* patchPtr_;  
     
     //- dual ppatch
     mutable PPatch* dualPatchPtr_;
     
     //- dual ppatch interpolation
     mutable PPatchInterpolation* dualPatchInterpPtr_;     
     
     //- addressing between dual patch and original patch
     //- dual face to original points (global)
     mutable labelList* mapFromDualFacesToPointsPtr_;
    
     //- original face centres to dual points (global)
     mutable labelList* mapFromFacesToDualPointsPtr_;   

     //- Disallow default bitwise copy construct
     dualPatch(const dualPatch&);

     //- Disallow default bitwise assignment
     void operator=(const dualPatch&);

    void makeDualPatch();


    template<class T>
    void sendToMaster
    (
        List<T>& Values
    ) const;

public:

     //- Runtime type information
     TypeName("dualPatch");

     //- Return the mesh sub-directory name (usually "triSurface")
     static word meshSubDir;


 // Constructors

     //- Construct from triSurface
     dualPatch
     (
         const triSurface& surf
     );
     
     //- Construct from faces and points
     dualPatch
     (
         const faceList& faces,
         const pointField& points
     ); 

 //- Destructor
 virtual ~dualPatch();

     //- Clear storage
     void clearOut();


 // Member Functions
     const PPatch& newPatch() const
     {
         return *dualPatchPtr_;
     }
     const PPatch& oldPatch() const
     {
         return *patchPtr_;
     }
     const PPatchInterpolation& dualPatchInterp() const
     {
        if(!dualPatchInterpPtr_)
        {
            dualPatchInterpPtr_ = new PPatchInterpolation(newPatch());
        }
        return *dualPatchInterpPtr_;
     }
     const faceList& dualFaces() const
     {
         return *dualPatchPtr_;
     }
     const pointField& dualPoints() const
     {
         return dualPatchPtr_->points();
     }
     const labelList& mapFromFacesToDualPoints() const
     {
         return *mapFromFacesToDualPointsPtr_;
     }
     const labelList& mapFromDualFacesToPoints() const
     {
         return *mapFromDualFacesToPointsPtr_;
     }
     
    // not used
    void makeCTSM
    (
        const labelList& includeFaces,
        const labelList& includeEdges,
        const labelList& includePoints,
        const labelList& boundaryEdges,
        const labelList& boundaryEdgesPatchID,
        
        labelList& includeDualFaces,
        labelList& includeDualEdges,
        labelList& includeDualPoints,
        labelList& boundaryDualEdges,
        labelList& boundaryDualEdgesPatchID        
        
    ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
