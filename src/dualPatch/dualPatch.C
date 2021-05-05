/*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
     \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
      \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "dualPatch.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "EdgeMap.H"
#include "triSurfaceFields.H"
#include "Time.H"
#include "PatchTools.H"
#include "processorCyclicPolyPatch.H"
#include "vtkSurfaceWriter.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
 defineTypeNameAndDebug(dualPatch, 0);
 //addToRunTimeSelectionTable(searchableSurface, dualPatch, dict);
 //word dualPatch::meshSubDir = "triSurface";
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dualPatch::dualPatch
(
    const triSurface& surf
)
:
patchPtr_(nullptr),
dualPatchPtr_(nullptr),
dualPatchInterpPtr_(nullptr),
mapFromDualFacesToPointsPtr_(nullptr),
mapFromFacesToDualPointsPtr_(nullptr)
{
    // this is for version 5.0
    faceList newFaces(surf.size());
    forAll(newFaces,I)
    {
        newFaces[I] = surf[I].triFaceFace();
    }
    patchPtr_ = new PPatch(newFaces,surf.points());
    
    //otherwise in dev or higher
    // patchPtr_ = new PPatch(surf.faces(,surf.points());
    makeDualPatch();
}

Foam::dualPatch::dualPatch
(
    const faceList& faces,
    const pointField& points
)
:
patchPtr_(new PPatch(faces,points)),
dualPatchPtr_(nullptr),
dualPatchInterpPtr_(nullptr),
mapFromDualFacesToPointsPtr_(nullptr),
mapFromFacesToDualPointsPtr_(nullptr)
{
    makeDualPatch();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

void Foam::dualPatch::makeDualPatch()
{
    PPatch surf(oldPatch());
    // find all boundary points
    const labelListList& edgeFaces = surf.edgeFaces();
    const edgeList& edges = surf.edges();
    labelHashSet boundaryPointsSet;
    forAll(edgeFaces,I)
    {
        if(edgeFaces[I].size()==1)
        {
            forAll(edges[I],II)
            {
                if(!boundaryPointsSet.found(edges[I][II]))
                {
                    boundaryPointsSet.insert(edges[I][II]);
                }
            }
        }
    }
    
    labelList boundaryPoints = boundaryPointsSet.toc();
    
    // new polygon Faces for dual mesh
    DynamicList<face> dualFaces;
    const labelListList& pointFaces = surf.pointFaces();
    const labelListList& faceFaces = surf.faceFaces();
    const pointField& dualPoints = surf.faceCentres();
    const pointField& pointNormals = surf.pointNormals();
    
    DynamicList<label> mFDFTP;
    
    forAll(pointFaces,I)
    {
        if(boundaryPointsSet.found(I)) continue;
        
        DynamicList<label> dualFace;
        const labelList& curFaces = pointFaces[I];
      
        labelHashSet unVisitedSet(curFaces);
        
        label faceID = curFaces[0];
        
        do
        {
            unVisitedSet.erase(faceID);
            dualFace.append(faceID);
            // find next faceID
            const labelList& neighFaces = faceFaces[faceID];
            forAll(neighFaces,II)
            {
                if(unVisitedSet.found(neighFaces[II]))
                {
                    faceID = neighFaces[II];
                    break;
                }
            }
        }
        while(unVisitedSet.size()>SMALL);
        dualFace.shrink();
        face newDualFace(dualFace);
        // check reverse
        scalar re = pointNormals[I]&newDualFace.normal(dualPoints);
        if(re<0)
        {
            newDualFace = newDualFace.reverseFace();
        }
        
        dualFaces.append(newDualFace);
        mFDFTP.append(I);
    }
    
    dualFaces.shrink();
    dualPatchPtr_ = new PPatch (dualFaces,dualPoints);
    mFDFTP.shrink();
    mapFromDualFacesToPointsPtr_ = new labelList(mFDFTP);
    
    mapFromFacesToDualPointsPtr_ = new labelList(dualPoints.size(),-1);
    labelList& mFFTDP = *mapFromFacesToDualPointsPtr_; // local points
    forAll(mFFTDP,I)
    {
        mFFTDP[I] = dualPatchPtr_->whichPoint(I); // from global to local
        if(mFFTDP[I]<0 and debug)
        {
            WarningIn
            (
                "dualPatch::makeDualPatch() const"
            )   << " mapFromFacesToDualPoints: face "<<I
                <<" does not have corresponding dual point"
                <<endl;
        }
    }  
}

Foam::dualPatch::~dualPatch()
{
    clearOut();
}

void Foam::dualPatch::clearOut()
{
    deleteDemandDrivenData(patchPtr_);
    deleteDemandDrivenData(dualPatchPtr_);
    deleteDemandDrivenData(dualPatchInterpPtr_);
    deleteDemandDrivenData(mapFromDualFacesToPointsPtr_);
    deleteDemandDrivenData(mapFromFacesToDualPointsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// not used
void Foam::dualPatch::makeCTSM
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
    
) const
{
    
}


// ************************************************************************* //
