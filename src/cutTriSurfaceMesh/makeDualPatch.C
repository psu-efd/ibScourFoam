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
     
Description
    make information necessary to make dual mesh 

SourceFiles
    makeDualPatch.C

\*---------------------------------------------------------------------------*/
#include "vtkSurfaceWriter.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cutTriSurfaceMesh::makeDualPatch() const 
{
    dualPatchPtr_ = new dualPatch(*this);
    dualPatch& patch = *dualPatchPtr_;
    PPatch surf = patch.newPatch();
    

    const labelList& includePoints = globalIncludePoints();// local pointID
    const labelList& boundaryEdges = globalBoundaryEdges();
    const labelListList& boundaryEdgesPatchID = globalBoundaryEdgesPatchID();
    
    const labelList& mFDFTP = patch.dualPatch::mapFromDualFacesToPoints();

    // reverse 
    labelList mFDFTP_reversed(this->localPoints().size(),-1);
    forAll(mFDFTP,I)
    {
        mFDFTP_reversed[mFDFTP[I]] = I;
    }
    
    // boundary edge points
    labelHashSet boundaryEdgesPointsSet;
    DynamicList<label> boundaryEdgesPoints;
    DynamicList<labelList> boundaryEdgesPointsPatchID;
    const edgeList& edges = this->edges();
    forAll(boundaryEdges,I)
    {
        label edgeID = boundaryEdges[I];
        forAll(edges[edgeID],II)
        {
            label pointID = edges[edgeID][II];
            if(!boundaryEdgesPointsSet.found(pointID))
            {
                boundaryEdgesPoints.append(pointID);
                boundaryEdgesPointsSet.insert(pointID);
                boundaryEdgesPointsPatchID.append(boundaryEdgesPatchID[I]);
            }
        }
    }  
    
    boundaryEdgesPoints.shrink();
    boundaryEdgesPointsPatchID.shrink();
    // find includeDualFaces, all face centers should be include points but not on edge
    labelHashSet includePointsSet(includePoints);
    labelHashSet includeDualFacesSet;
    
    forAll(includePoints,I)
    {
        label pointID = includePoints[I];
        if(!boundaryEdgesPointsSet.found(pointID))
        {
            label faceID = mFDFTP_reversed[pointID];
            if(faceID<0) Info<<I<<" mFDFTP_reversed"<<endl;
            includeDualFacesSet.insert(faceID);
        }
    }
    includeDualFacesPtr_ = new labelList(includeDualFacesSet.toc());
    
    // find boundaryDualEdges
    labelHashSet boundaryDualEdgesSet;
    DynamicList<label> boundaryDualEdges;
    DynamicList<labelList> boundaryDualEdgesPatchID;
        
    const labelListList& faceEdges = surf.faceEdges();
    const labelListList& edgeFaces = surf.edgeFaces();
    forAll(boundaryEdgesPoints,I)
    {    
        label pointID = boundaryEdgesPoints[I];
        label faceID = mFDFTP_reversed[pointID];
        const labelList& curEdges = faceEdges[faceID];
        // find edges neighbouring includeDualFaces
        forAll(curEdges,II)
        {
            label edgeID = curEdges[II];
            labelList curFaces = edgeFaces[edgeID];
            if(curFaces.size()<2) continue;
            label neighFaceID = curFaces[0];
            if(faceID == neighFaceID) neighFaceID = curFaces[1];
            if(includeDualFacesSet.found(neighFaceID))
            {
                boundaryDualEdgesSet.insert(edgeID);
                boundaryDualEdges.append(edgeID);
                boundaryDualEdgesPatchID.append(boundaryEdgesPointsPatchID[I]);
            }
        }
    }
    
    boundaryDualEdges.shrink();
    boundaryDualEdgesPatchID.shrink();

    boundaryDualEdgesPtr_ = new labelList(boundaryDualEdges);
    boundaryDualEdgesPatchIDPtr_ = new labelListList(boundaryDualEdgesPatchID);
    
    // find includeDualEdges and includeDualPoints
    const labelList& includeDualFaces = *includeDualFacesPtr_;
    labelHashSet includeDualEdgesSet;
    labelHashSet includeDualPointsSet;
    forAll(includeDualFaces,I)
    {
        label faceID =  includeDualFaces[I];
        const labelList& faceEdges = surf.faceEdges()[faceID];

        forAll(faceEdges,II)
        {
            label edgeID = faceEdges[II];
            if(!includeDualEdgesSet.found(edgeID))
            {
                includeDualEdgesSet.insert(edgeID);
            }
        }

        const face& myFaces = surf.localFaces()[faceID];

        forAll(myFaces,II)
        {
            label pointID = myFaces[II];
            // use local point
            if(!includeDualPointsSet.found(pointID))
            {
                includeDualPointsSet.insert(pointID);
            }
        }        
    }
    
    includeDualPointsPtr_ = new labelList(includeDualPointsSet.toc());    
    includeDualEdgesPtr_ = new labelList(includeDualEdgesSet.toc());        
}



// ************************************************************************* //

