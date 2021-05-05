/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

SourceFiles
    immersedBoundaryFvMeshPostEvaluation.C
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvMesh.H"
#include "MeshWave.H"
#include "SortableList.H"
#include "fvCFD.H"
#include "turbulenceModel.H"
#include "vectorTools.H"
#include "fileOperation.H"
#include "turbulentTransportModel.H"
#include "vtkSurfaceWriter.H"
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::immersedBoundaryFvMesh::makeDualMesh()const
{
    dualMeshListPtr_ = new PtrList<fvMesh> (objectsList().size());
    oldToNewDualEdgeMapListPrt_ = new PtrList<labelList>(objectDictList().size());
    newToOldDualEdgeMapListPrt_ = new PtrList<labelList>(objectDictList().size());

    forAll(objectsList(),objectID)
    {

        if(dual(objectID))
        {
            makeDualMesh(objectID);
        }
    }
}


void Foam::immersedBoundaryFvMesh::makeDualMesh(const label& objectID)const
{
    const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];
    const PtrList<triSurface>& addSurfs= addObjectsList();
    const PtrList<word>& addObjectsName = addObjectsNameList();
    const dualPatch& dPatch = cTSM.newDualPatch();
    const PPatch& surf = dPatch.newPatch();
    const labelList& includeFaces = cTSM.includeDualFaces();
    const labelList& includeEdges = cTSM.includeDualEdges();
    const labelList& includePoints = cTSM.includeDualPoints();// local pointID
    const labelList& boundaryEdges = cTSM.boundaryDualEdges();
    const labelListList& boundaryEdgesPatchID = cTSM.boundaryDualEdgesPatchID();
    
    
    
    const faceList& surf_faces = surf;
    pointField surf_points = surf.points();
    
    // If IB is set as sediment bed, 
    if(objectDictList()[objectID].found("sediment"))
    {
        forAll(surf_points,I)
        {
            surf_points[I].z() = 0.0;
        }
    }
    
    PPatch surf_org(surf_faces,surf_points);
    
    const pointField& oldPoints = surf_org.localPoints();

    const pointField& pointNormals = surf_org.pointNormals();

    //make new points and new bottom faces
    pointField newPoints(includePoints.size());
    faceList newBottomFaces(includeFaces.size());
    labelList newBottomOwners(includeFaces.size(),-1);
    labelList oldToNewPtMap(oldPoints.size(),-1);
    forAll(includePoints,I)// locally points
    {
        oldToNewPtMap[includePoints[I]] = I;
        newPoints[I] = oldPoints[includePoints[I]];
    }
    labelList oldToNewFcMap(surf_org.size(),-1);
    forAll(includeFaces,I)
    {
        face myFace = surf_org.localFaces()[includeFaces[I]];
        forAll(myFace,II)
        {
            myFace[II] = oldToNewPtMap[myFace[II]];
        }
        newBottomFaces[I] = myFace;
        oldToNewFcMap[includeFaces[I]] = I;
        newBottomOwners[I] = I;
    }
    

    //make new added points and new top faces
    pointField newAddPoints(newPoints);
    faceList newTopFaces(newBottomFaces);

    //0: extrude according to face/point normal
    //1: extrude according to a defined value
    //2: extrude according to gravity, applied in sediment
    label extrudeType(objectDictList()[objectID].lookupOrDefault<label>("extrudeType",2));

    scalar extrudeDirectionMag(objectDictList()[objectID].lookupOrDefault<scalar>("extrudeDirectionMag",0.02));

    vectorField extrudeDirection(newPoints.size(),pTraits<vector>::zero);
    
    vectorField cellCentresTranslation(surf_org.size(),pTraits<vector>::zero);
    
    if(extrudeType == 2 and !objectDictList()[objectID].found("sediment"))
    {
        extrudeType = 0;
    } 
    
    if(extrudeType == 0)
    {
        forAll(pointNormals,I)
        {
            if(I>-1 and oldToNewPtMap[I]<extrudeDirection.size() and oldToNewPtMap[I]>-1)
            {
                extrudeDirection[oldToNewPtMap[I]] = pointNormals[I];
            }
        }
        
        cellCentresTranslation = surf_org.faceNormals()/mag(surf_org.faceNormals());
    }
    else if(extrudeType == 1)
    {
        extrudeDirection = objectDictList()[objectID].lookupOrDefault<vector>
            ("extrudeDirection",pTraits<vector>::zero);
            
        cellCentresTranslation = objectDictList()[objectID].lookupOrDefault<vector>
            ("extrudeDirection",pTraits<vector>::zero);
    }
    else if(extrudeType == 2)
    {
        // assuming plane is horizontal    
        const dictionary& dict = objectDictList()[objectID].subDict("sediment");
        vector gravity(dict.lookup("gravity"));
        scalar ed = 1;
        if(!ibTriNormalsFlipList()[objectID])
        {
            ed = -1;
        }
        extrudeDirection = ed*gravity;
        cellCentresTranslation = ed*gravity;
    }
    
    vectorField translation(extrudeDirection/mag(extrudeDirection)*extrudeDirectionMag); // point translation 

    cellCentresTranslation = 
        cellCentresTranslation/mag(cellCentresTranslation)*extrudeDirectionMag; // cell centres translation

    newAddPoints = newAddPoints + translation;
    labelList newTopOwners(newBottomOwners);

    // points
    newPoints.append(newAddPoints);

    const label diff_pointID = newAddPoints.size();
    
    forAll(newTopFaces,I)
    {
        face& myFace = newTopFaces[I];

        forAll(myFace,II)
        {
            myFace[II] +=diff_pointID;
        }
        // reverse top face normal
        myFace = myFace.reverseFace();
    }

    // build faces from includeEdges, build owner and neighbour for these faces
    // includeEgdges also include boundaryEdges
    faceList newInternalFaces(includeEdges.size()-boundaryEdges.size());
    labelList newInternalOwners(newInternalFaces.size(),-1);
    labelList newInternalNeighbours(newInternalFaces.size(),-1);
    labelHashSet boundaryEdgesSet(boundaryEdges);
    label newNInternalFaces = 0;
    labelList pointCounter(newPoints.size(),0);

    labelList oldToNewEdgeMap(surf_org.nEdges(),-1);
    // internal faces
    forAll(includeEdges,I)
    {
        label edgeID = includeEdges[I];
        if(!boundaryEdgesSet.found(edgeID))
        {
            const edge& e = surf_org.edges()[edgeID];
            labelList myFace(4); // alway has four points
            myFace = -1;
            // point ID needs to be converted to new pointID
            label pointID0 = oldToNewPtMap[e.start()];
            label pointID1 = oldToNewPtMap[e.end()];
            if(pointID0<0 or pointID1<0)
            {
                WarningIn
                (
                    "immersedBoundaryFvMesh::extrudedMesh() const"
                )   << "include edge "<<edgeID
                    <<" does not have oldToNewPtMap for one of its ends "
                    <<pointID0<<" "<<pointID1<<endl;
            }

            myFace[0] = pointID0;
            myFace[1] = pointID1;
            myFace[2] = pointID1+diff_pointID;
            myFace[3] = pointID0+diff_pointID;
            pointCounter[pointID0] += 1;
            pointCounter[pointID1] += 1;
            pointCounter[pointID1+diff_pointID] += 1;
            pointCounter[pointID0+diff_pointID] += 1;
            face newMyFace(myFace);

            newInternalFaces[newNInternalFaces] = newMyFace;
            
            // face ID needs to be converted to new faceID
            label faceID0 = surf_org.edgeFaces()[edgeID][0];
            label faceID1 = surf_org.edgeFaces()[edgeID][1];
            faceID0 = oldToNewFcMap[faceID0];
            faceID1 = oldToNewFcMap[faceID1];
            if(faceID0<0 or faceID1<0)
            {
                WarningIn
                (
                    "immersedBoundaryFvMesh::extrudedMesh() const"
                )   << "include edge or internal face "<<edgeID
                    <<" does not have oldToNewFcMap for one of neighbours "
                    <<faceID0<<" "<<faceID1<<endl;
            }
            newInternalOwners[newNInternalFaces] = min(faceID0,faceID1);
            newInternalNeighbours[newNInternalFaces] = max(faceID0,faceID1);
            oldToNewEdgeMap[edgeID] = newNInternalFaces;
            newNInternalFaces++;
        }
    }

    // boundary faces
     label nAddsurf = addSurfs.size();
    List<DynamicList<face>> newPatchFaces(this->boundary().size()+nAddsurf);
    List<DynamicList<label>> newPatchOwners(this->boundary().size()+nAddsurf);
    labelList newPatchNFaces(this->boundary().size()+nAddsurf,0);
    labelHashSet includeFacesSet(includeFaces);
    List<DynamicList<label>> oldToNewEdgeMapPatches(this->boundary().size()+nAddsurf);

    

    forAll(boundaryEdges,I)
    {
        label edgeID = boundaryEdges[I];
        labelList patchIDs = boundaryEdgesPatchID[I];
        label patchID=0;
        
        if(patchIDs.size()>1 and max(patchIDs) > this->boundary().size()-1)
        {
           patchID = this->boundary().size()+nAddsurf-1;
        }
        else
        {
           patchID = patchIDs[0];
        }

        const edge& e = surf_org.edges()[edgeID];
        labelList myFace(4); // alway has four points
        myFace = -1;
        // point ID needs to be converted to new pointID
        label pointID0 = oldToNewPtMap[e.start()];
        label pointID1 = oldToNewPtMap[e.end()];
        if(pointID0<0 or pointID1<0)
        {
            Info<<pointID0<<" "<<pointID1<<endl;
        }
        myFace[0] = pointID0;
        myFace[1] = pointID1;
        myFace[2] = pointID1+diff_pointID;
        myFace[3] = pointID0+diff_pointID;
        pointCounter[pointID0] += 1;
        pointCounter[pointID1] += 1;
        pointCounter[pointID1+diff_pointID] += 1;
        pointCounter[pointID0+diff_pointID] += 1;
        face newMyFace(myFace);
        newPatchFaces[patchID].append(newMyFace);
        oldToNewEdgeMapPatches[patchID].append(edgeID);
        // face ID needs to be converted to new faceID
        label faceID0 = surf_org.edgeFaces()[edgeID][0];
        label faceID1 = surf_org.edgeFaces()[edgeID][1];
        if(faceID0<0 or faceID1<0)
        {
            WarningIn
            (
                "immersedBoundaryFvMesh::extrudedMesh() const"
            )   << "boundary edge or boundary face "<<edgeID
                <<" does not have oldToNewFcMap for one of neighbours "
                <<faceID0<<" "<<faceID1<<endl;
        }
        if(includeFacesSet.found(faceID0))
        {
            newPatchOwners[patchID].append(oldToNewFcMap[faceID0]);
        }
        else if(includeFacesSet.found(faceID1))
        {
            newPatchOwners[patchID].append(oldToNewFcMap[faceID1]);
        }
        else
        {
            Info<<"something goes wrong in making extrudedMesh"<<endl;
        }
        newPatchNFaces[patchID]++;
    }

    forAll(pointCounter,I)
    {
        if(pointCounter[I]>6 and debug == 2)
        {
              Info<<"pointCounter is too large "<<I<<endl;
          }
    }

    // faces
    DynamicList<face> newFaces;
    newFaces.append(newInternalFaces);

    // owners and neighbours
    DynamicList<label> newOwners;
    newOwners.append(newInternalOwners);

    DynamicList<label> newNeighbours;
    newNeighbours.append(newInternalNeighbours);

    // pre-add boundary faces and owners
    DynamicList<label> startFaceIDList;
    DynamicList<label> faceSizeList;
    DynamicList<word> patchNameList;
    DynamicList<label> oldPatchIDList;
    startFaceIDList.append(newFaces.size());
    faceSizeList.append(newBottomFaces.size()+newTopFaces.size());
    patchNameList.append("bottomAndTop");
    oldPatchIDList.append(-1);

    newFaces.append(newBottomFaces);
    newFaces.append(newTopFaces);
    newOwners.append(newBottomOwners);
    newOwners.append(newTopOwners);

    forAll(newPatchFaces, oldPatchID)
    {
        newPatchFaces[oldPatchID].shrink();
        newPatchOwners[oldPatchID].shrink();
        if(newPatchFaces[oldPatchID].size()>0)
        {
            if (oldPatchID < this->boundary().size())
            {
                patchNameList.append(this->boundaryMesh()[oldPatchID].name());
            }
            else
            {
                patchNameList.append(addObjectsName[oldPatchID-this->boundary().size()]);
            }
            
            startFaceIDList.append(newFaces.size());
            faceSizeList.append(newPatchFaces[oldPatchID].size());
            newFaces.append(newPatchFaces[oldPatchID]);
            newOwners.append(newPatchOwners[oldPatchID]);
            oldPatchIDList.append(oldPatchID);
        }
    }
 
    // make sure face normal is correct
    newFaces.shrink();
    newOwners.shrink();
    newNeighbours.shrink();
    startFaceIDList.shrink();
    faceSizeList.shrink();
    patchNameList.shrink();
    oldPatchIDList.shrink();
   
    const pointField& faceCentres = surf_org.faceCentres();
    forAll(newFaces, I)
    {
        // face normal direction
        face& myFace = newFaces[I];
        
        vector normal = myFace.normal(newPoints);
        vector centre = myFace.centre(newPoints);

        
        vector own = faceCentres[includeFaces[newOwners[I]]]
            +cellCentresTranslation[includeFaces[newOwners[I]]]*0.5;

        scalar re = normal&(centre-own);
        if(re<0)
        {
            myFace = myFace.reverseFace();
        }
    }

    Xfer<pointField> XferNewPoints(newPoints);
    
    Foam::fvMesh mesh
    (
        IOobject
        (
            objectNames(objectID)+"_dual",
            time().constant(),
            time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        XferNewPoints,
        newFaces.xfer(),
        newOwners.xfer(),
        newNeighbours.xfer()
    );
        

    // make new patches
    List<polyPatch*> newPatchPtrList(patchNameList.size());
    label newPatchID = 0;

    // make empty boundary bottomAndTop, now wall
    emptyPolyPatch emptyPatch
    (
        patchNameList[0], // name
        startFaceIDList[0], // size
        faceSizeList[0],
        newPatchID,
        mesh.boundaryMesh(),
        "empty"
    );
    
  
    
    // clone empty patch as the 0th in the new patch list
    newPatchPtrList[newPatchID] = emptyPatch.clone
        (
            mesh.boundaryMesh(),
            newPatchID,
            faceSizeList[0],
            startFaceIDList[0]
        ).ptr();
        
    

    forAll(newPatchPtrList,patchI)
    {
        
        if(patchI>0)
        {
            const label& oldPatchID = oldPatchIDList[patchI];
            const label& startFaceID = startFaceIDList[patchI];
            const label& faceSize = faceSizeList[patchI];
            if (oldPatchID <  this->boundary().size())
            {
                const polyPatch& oldPatch = this->boundaryMesh()[oldPatchID];
                newPatchPtrList[patchI] = oldPatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchI,
                        faceSize,
                        startFaceID
                    ).ptr();
               
             }
             else
             {
                wallPolyPatch wallPatch
                (
                    patchNameList[patchI], // name
                    startFaceIDList[patchI], // size
                    faceSizeList[patchI],
                    patchI,
                    mesh.boundaryMesh(),
                    "wall"
                );
                    newPatchPtrList[patchI] = wallPatch.clone
                    (
                        mesh.boundaryMesh(),
                        patchI,
                        faceSizeList[patchI],
                        startFaceIDList[patchI]
                    ).ptr();
                     
             }  
            oldToNewEdgeMapPatches[oldPatchID].shrink();    
            forAll(oldToNewEdgeMapPatches[oldPatchID],I)
            {
                oldToNewEdgeMap[oldToNewEdgeMapPatches[oldPatchID][I]] = startFaceID+I;
            }
        }
    }
    
    mesh.addPatches(newPatchPtrList);
    
    

    mesh.write();

    oldToNewDualEdgeMapListPrt_->set
    (
        objectID,
        new labelList(oldToNewEdgeMap)
    );

    labelList newToOldEdgeMap(mesh.nFaces(),-1);

    forAll(oldToNewEdgeMap,I)
    {
        if(oldToNewEdgeMap[I]>-1)
        {
            newToOldEdgeMap[oldToNewEdgeMap[I]]=I;
        }
    }

    newToOldDualEdgeMapListPrt_->set
    (
        objectID,
        new labelList(newToOldEdgeMap)
    );
    

    if(debug)
    {
        Info<<mesh.checkMesh(true)<<endl;
    }

    // The following export is very necessary for completeness
    mesh.lookupObject<pointZoneMesh>("pointZones").write();
    mesh.lookupObject<faceZoneMesh>("faceZones").write();
    mesh.lookupObject<cellZoneMesh>("cellZones").write();

    if(Pstream::parRun() and Pstream::myProcNo()==Pstream::masterNo())
    {
        parallelWriteMesh(mesh.lookupObject<pointIOField>("points"));
        parallelWriteMesh(mesh.lookupObject<faceCompactIOList>("faces"));
        parallelWriteMesh(mesh.lookupObject<labelIOList>("owner"));
        parallelWriteMesh(mesh.lookupObject<labelIOList>("neighbour"));
        parallelWriteMesh(mesh.lookupObject<polyBoundaryMesh>("boundary"));
        parallelWriteMesh(mesh.lookupObject<pointZoneMesh>("pointZones"));
        parallelWriteMesh(mesh.lookupObject<faceZoneMesh>("faceZones"));
        parallelWriteMesh(mesh.lookupObject<cellZoneMesh>("cellZones"));
    }
    // Read again and save it to pointer
    // Cannot directly use the mesh made above
    dualMeshListPtr_->set
    ( 
        objectID,
        new fvMesh
        (
            IOobject
            (
                objectNames(objectID)+"_dual",
                time().constant(),
                time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            XferNewPoints,
            newFaces.xfer(),
            newOwners.xfer(),
            newNeighbours.xfer()
        )
    );
    
}


void Foam::immersedBoundaryFvMesh::makeDualMeshAddressing(label objectID) const
{
    const double Oldtime1=time().elapsedCpuTime();

    const triSurface& surf = objectsList()[objectID];
    const pointField& pts = surf.localPoints(); 
    const labelListList& pointEdges = surf.pointEdges(); 
    const labelListList& faceEdges = surf.faceEdges(); 
    const labelListList& pointFaces = surf.pointFaces();   
    const labelListList& edgeFaces = surf.edgeFaces(); 
    const labelList&  boundaryPoints = surf.boundaryPoints(); 
    const pointField& faceCentres = surf.faceCentres();	
    pointField newPoints(faceCentres);
    labelHashSet boundaryPointsSet(boundaryPoints);
    labelListList dualMeshAddressing(pts.size());
  
    faceList newFaces(pts.size());
    forAll(pts,pointID)
    {
        if(!boundaryPointsSet.found(pointID))
        {
             const labelList edges = pointEdges[pointID];
             const labelList faces = pointFaces[pointID];
             DynamicList<label> newFaces(faces.size());            
             labelHashSet edgesSet(edges);
             label nFaces = 0;
             label lastEdge = edges[0]; 
             label lastFace = edgeFaces[lastEdge][0];                    
             newFaces.append(lastFace);
             nFaces++;  
   

             bool isAllFaceSweeped = false;
             while (!isAllFaceSweeped)
             {  
                  const labelList& myEdges = faceEdges[lastFace];

                  forAll(myEdges,edgeID)
                  {
                        label eI = myEdges[edgeID];
                        
    
                        if((eI != lastEdge) and (edgesSet.found(eI)))
                        {
                            lastEdge = eI;
         
                            const labelList& myEdgeFaces = edgeFaces[eI];
       
                            forAll(myEdgeFaces,fI)
                            {
                                if(lastFace != myEdgeFaces[fI])
                                {
                                                                                                  
                                     labelHashSet   newFacesSet(newFaces);
                                     if(newFacesSet.found(myEdgeFaces[fI]))
                                     {
                                        isAllFaceSweeped = true;
        
                                     }
                                     else
                                     {
                                        newFaces.append(myEdgeFaces[fI]);
      
                                        lastFace = myEdgeFaces[fI];
                                        nFaces++;
                                     } 
                                     break; 
                                }
                                
                            }
                            break;
                        }
                  }
             }
       
             newFaces.shrink();
             dualMeshAddressing[pointID] = newFaces;

             dualMeshAddressingListPtr_->set
             (
                objectID,
                new labelListList(dualMeshAddressing)
             );
        }
        
    } 
    const double Oldtime2=time().elapsedCpuTime();
    Info<<"dualMeshAddr Executation Time = "<<Oldtime2-Oldtime1<< " s"<<endl;
    Info<< " dualMeshAddr is finished "  <<endl; 
}
    
    
 
// ************************************************************************* //

