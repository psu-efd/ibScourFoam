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

#include "cutTriSurfaceMesh.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "EdgeMap.H"
#include "triSurfaceFields.H"
#include "triSurfaceTools.H"
#include "Time.H"
#include "PatchTools.H"
#include "processorCyclicPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
 defineTypeNameAndDebug(cutTriSurfaceMesh, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutTriSurfaceMesh::cutTriSurfaceMesh
(
    const triSurface& surf,
    const polyMesh& mesh
)
:
triSurface(surf),
includeFacesPtr_(nullptr),
includeEdgesPtr_(nullptr),
includePointsPtr_(nullptr),
boundaryFacesPtr_(nullptr),
boundaryFacesPatchIDPtr_(nullptr),
boundaryFacesNearestPatchFaceInfoPtr_(nullptr),
boundaryEdgesPtr_(nullptr),
boundaryEdgesPatchIDPtr_(nullptr),
boundaryEdgesNearestPatchFaceInfoPtr_(nullptr),
outsidePointsPtr_(nullptr),
outsidePointsStencilPtr_(nullptr),
outsidePointsStencilWeightPtr_(nullptr),
globalIncludeFacesPtr_(nullptr),
globalIncludeEdgesPtr_(nullptr),
globalIncludePointsPtr_(nullptr),
globalBoundaryFacesPtr_(nullptr),
globalBoundaryFacesPatchIDPtr_(nullptr),
globalBoundaryFacesPatchIDProcPtr_(nullptr),
globalBoundaryFacesNearestPatchFaceInfoPtr_(nullptr),
globalBoundaryEdgesPtr_(nullptr),
globalBoundaryEdgesPatchIDPtr_(nullptr),
globalBoundaryEdgesPatchIDProcPtr_(nullptr),
globalBoundaryEdgesNearestPatchFaceInfoPtr_(nullptr),
dualPatchPtr_(nullptr),
includeDualFacesPtr_(nullptr),
includeDualEdgesPtr_(nullptr),
includeDualPointsPtr_(nullptr),
boundaryDualEdgesPtr_(nullptr),
boundaryDualEdgesPatchIDPtr_(nullptr)
{
    makeBoundaryFacesAndEdges(mesh);
    const double Oldtime0=mesh.time().elapsedCpuTime();
    makeGlobal(mesh); //for mpi

    const double Oldtime1=mesh.time().elapsedCpuTime();
    makeOutsidePointsStencil();
    const double Oldtime2=mesh.time().elapsedCpuTime();

    if(debug)
    {
        Info<<"makeGlobal Time = "<<Oldtime1-Oldtime0<< " s"<<endl;
        Info<<"makeOutsidePointsStencil Time = "<<Oldtime2-Oldtime1<< " s"<<endl;
    }

}







Foam::cutTriSurfaceMesh::cutTriSurfaceMesh
(
    const triSurface& surf,
    const polyMesh& mesh,
    const PtrList<triSurface>& addSurfs
)
:
triSurface(surf),
includeFacesPtr_(nullptr),
includeEdgesPtr_(nullptr),
includePointsPtr_(nullptr),
boundaryFacesPtr_(nullptr),
boundaryFacesPatchIDPtr_(nullptr),
boundaryFacesNearestPatchFaceInfoPtr_(nullptr),
boundaryEdgesPtr_(nullptr),
boundaryEdgesPatchIDPtr_(nullptr),
boundaryEdgesNearestPatchFaceInfoPtr_(nullptr),
outsidePointsPtr_(nullptr),
outsidePointsStencilPtr_(nullptr),
outsidePointsStencilWeightPtr_(nullptr),
globalIncludeFacesPtr_(nullptr),
globalIncludeEdgesPtr_(nullptr),
globalIncludePointsPtr_(nullptr),
globalBoundaryFacesPtr_(nullptr),
globalBoundaryFacesPatchIDPtr_(nullptr),
globalBoundaryFacesPatchIDProcPtr_(nullptr),
globalBoundaryFacesNearestPatchFaceInfoPtr_(nullptr),
globalBoundaryEdgesPtr_(nullptr),
globalBoundaryEdgesPatchIDPtr_(nullptr),
globalBoundaryEdgesPatchIDProcPtr_(nullptr),
globalBoundaryEdgesNearestPatchFaceInfoPtr_(nullptr),
dualPatchPtr_(nullptr),
includeDualFacesPtr_(nullptr),
includeDualEdgesPtr_(nullptr),
includeDualPointsPtr_(nullptr),
boundaryDualEdgesPtr_(nullptr),
boundaryDualEdgesPatchIDPtr_(nullptr)
{
    makeBoundaryFacesAndEdges(mesh,addSurfs);
    const double Oldtime0=mesh.time().elapsedCpuTime();
    makeGlobal(mesh,addSurfs); //for mpi

    const double Oldtime1=mesh.time().elapsedCpuTime();
    
    makeOutsidePointsStencil();
    
    const double Oldtime2=mesh.time().elapsedCpuTime();

    if(debug)
    {
        Info<<"makeGlobal Time = "<<Oldtime1-Oldtime0<< " s"<<endl;
        Info<<"makeOutsidePointsStencil Time = "<<Oldtime2-Oldtime1<< " s"<<endl;
    }


}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

void Foam::cutTriSurfaceMesh::makeGlobal
(
    const polyMesh& procMesh
) const
{
 

    if (!Pstream::parRun())
    {
        globalIncludeFacesPtr_ = new labelList(includeFaces());
        globalIncludeEdgesPtr_ = new labelList(includeEdges());
        globalIncludePointsPtr_ = new labelList(includePoints());
        globalBoundaryFacesPtr_ = new labelList(boundaryFaces());
        globalBoundaryFacesPatchIDPtr_ = new labelListList(boundaryFacesPatchID());
        labelListList globalBoundaryFacesPatchIDProc(boundaryFacesPatchID().size());
        forAll(globalBoundaryFacesPatchIDProc,I)
        {
            globalBoundaryFacesPatchIDProc[I].setSize(boundaryFacesPatchID()[I].size(),0);
        }
        globalBoundaryFacesPatchIDProcPtr_ = new labelListList(globalBoundaryFacesPatchIDProc);

        globalBoundaryEdgesPtr_ = new labelList(boundaryEdges());
        globalBoundaryEdgesPatchIDPtr_ = new labelListList(boundaryEdgesPatchID());
        labelListList globalBoundaryEdgesPatchIDProc(boundaryEdgesPatchID().size());
        forAll(globalBoundaryEdgesPatchIDProc,I)
        {
            globalBoundaryEdgesPatchIDProc[I].setSize(boundaryEdgesPatchID()[I].size(),0);
        }
        globalBoundaryEdgesPatchIDProcPtr_ = new labelListList(globalBoundaryEdgesPatchIDProc);

        return;
    }
    const labelList& boundaryFaces = this->boundaryFaces();
    labelListList& boundaryFacesPatchID = *boundaryFacesPatchIDPtr_;

    labelHashSet procPatchListSet;

    // replace processorCyclicPolyPatchID to its referPatchID
    labelList procPatchIDReplace(procMesh.boundaryMesh().size());

    forAll(procMesh.boundaryMesh(),patchID)
    {
        procPatchIDReplace[patchID]=patchID;
        if
        (
            isA<processorCyclicPolyPatch>
            (
                procMesh.boundaryMesh()[patchID]
            )
        )
        {
            const processorCyclicPolyPatch& pcPP =
                        refCast<const processorCyclicPolyPatch>(procMesh.boundaryMesh()[patchID]);
            procPatchIDReplace[patchID]=pcPP.referPatchID();
        }
        else if
        (
            isA<processorPolyPatch>
            (
                procMesh.boundaryMesh()[patchID]
            )
        )
        {
            procPatchListSet.insert(patchID);
        }
    }

    DynamicList<label> procBoundaryFaces;
    DynamicList<labelList> procBoundaryFacesPatchID;


    forAll(boundaryFaces,I)
    {
        label faceID = boundaryFaces[I];
        labelList patchIDs = boundaryFacesPatchID[I];
        bool unChanged = false;
        labelHashSet newPatchIDSet;
        forAll(patchIDs, ID)
        {
            patchIDs[ID]=procPatchIDReplace[patchIDs[ID]];
            if(!procPatchListSet.found(patchIDs[ID]))
            {
                newPatchIDSet.insert(patchIDs[ID]);
                unChanged = true;
            }
        }
        if(unChanged)
        {
            procBoundaryFaces.append(faceID);
            procBoundaryFacesPatchID.append(newPatchIDSet.toc());
        
        }

    }
   
    procBoundaryFaces.shrink();
    procBoundaryFacesPatchID.shrink();

    
    // transfer to each processor, not only master
    List<labelList> procBoundaryFacesList(Pstream::nProcs());
    List<labelListList> procBoundaryFacesPatchIDList(Pstream::nProcs());


    procBoundaryFacesList[Pstream::myProcNo()]=procBoundaryFaces;
    procBoundaryFacesPatchIDList[Pstream::myProcNo()]=procBoundaryFacesPatchID;



    Pstream::gatherList(procBoundaryFacesList);
    Pstream::scatterList(procBoundaryFacesList);

    Pstream::gatherList(procBoundaryFacesPatchIDList);
    Pstream::scatterList(procBoundaryFacesPatchIDList);



    DynamicList<label> globalBoundaryFaces;
    DynamicList<labelList> globalBoundaryFacesPatchID;
    DynamicList<labelList> globalBoundaryFacesPatchIDProc;


    HashTable<label,label> globalBoundaryFacesSet;
    label globalBoundaryFacesCounts=0;

    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        forAll(procBoundaryFacesList[procI],I)
        {
            label faceID = procBoundaryFacesList[procI][I];
            labelList PatchIDs = procBoundaryFacesPatchIDList[procI][I];

            // make sure remove duplicate faceID
            if(!globalBoundaryFacesSet.found(faceID))
            {
                globalBoundaryFacesSet.insert(faceID,globalBoundaryFacesCounts);
                globalBoundaryFacesCounts++;
                globalBoundaryFaces.append(faceID);
                globalBoundaryFacesPatchID.append(PatchIDs);
                labelList PatchIDsProc(PatchIDs.size(),procI);
                globalBoundaryFacesPatchIDProc.append(PatchIDsProc);

            }
            else // add duplicate face patch IDs and procs
            {
                label newID = globalBoundaryFacesSet.find(faceID)();

                DynamicList<label> tempDList;
                tempDList.append(globalBoundaryFacesPatchID[newID]);
                tempDList.append(PatchIDs);
                tempDList.shrink();
                globalBoundaryFacesPatchID[newID]=tempDList;
                labelList PatchIDsProc(PatchIDs.size(),procI);
                DynamicList<label> tempDList1;
                tempDList1.append(globalBoundaryFacesPatchIDProc[newID]);
                tempDList1.append(PatchIDsProc);
                tempDList1.shrink();
                globalBoundaryFacesPatchIDProc[newID]=tempDList1;
            }
        }
    }
    
    globalBoundaryFaces.shrink();
    globalBoundaryFacesPatchID.shrink();
    globalBoundaryFacesPatchIDProc.shrink();
    globalBoundaryFacesPtr_ = new labelList(globalBoundaryFaces);
    globalBoundaryFacesPatchIDPtr_ = new labelListList(globalBoundaryFacesPatchID);
    globalBoundaryFacesPatchIDProcPtr_ = new labelListList(globalBoundaryFacesPatchIDProc);

    labelHashSet globalIncludeFacesSet;

    label startID=0;
    List<labelHashSet> globalIncludeFacesSetList;
    label nFace=0;
    const triSurface& surf(*this);

    while(nFace<(surf.size()-globalBoundaryFaces.size()))
    {   
        labelHashSet nbFaceIDSet;
        labelList old_nbFaceID;
        old_nbFaceID.append(startID);
        labelHashSet includeFacesLocalSet;
        const labelListList& faceFaces = surf.faceFaces();

        while(old_nbFaceID.size()>0)
        {
            nbFaceIDSet.clear();
            forAll(old_nbFaceID, faceI)
            {
                startID=old_nbFaceID[faceI];
                forAll(faceFaces[startID], faceII)
                {
                    label faceID = faceFaces[startID][faceII];
                    if(!globalIncludeFacesSet.found(faceID) and !globalBoundaryFacesSet.found(faceID) )
                    {
                        nFace++;
                        includeFacesLocalSet.insert(faceID);
                        globalIncludeFacesSet.insert(faceID);
                        nbFaceIDSet.insert(faceID);
                    }
                }
            }
            old_nbFaceID=nbFaceIDSet.toc();
        }

        globalIncludeFacesSetList.append(includeFacesLocalSet);
        forAll(surf,faceID)
        {
            if(!globalIncludeFacesSet.found(faceID) and !globalBoundaryFacesSet.found(faceID) )
            {
                startID=faceID;
                if(includeFacesLocalSet.size()<SMALL)
                {
                        globalIncludeFacesSetList[globalIncludeFacesSetList.size()-1].insert(faceID);
                        globalIncludeFacesSet.insert(faceID);
                        nFace++;
                }
                break;
            }
        }
    }
    globalIncludeFacesSet.clear();
    labelList globalIncludeFaces;

    forAll(globalIncludeFacesSetList,I)
    {
        labelHashSet tempSet=globalIncludeFacesSetList[I];
        label faceID = tempSet.toc()[0];
        scalar indicator = -1;

        const face& myFaces = surf[faceID];

        forAll(myFaces,pointI)
        {

            if(procMesh.findCell(surf.points()[myFaces[pointI]])>-1)
            {
                indicator=1;
            }
        }
     
        if(indicator>0)
        {
            globalIncludeFaces.append(tempSet.toc());
        }
    }

    labelListList procInlcudeFaces(Pstream::nProcs());
    procInlcudeFaces[Pstream::myProcNo()]=globalIncludeFaces;
    Pstream::gatherList(procInlcudeFaces);
    Pstream::scatterList(procInlcudeFaces);
    forAll(procInlcudeFaces,I)
    {
        labelHashSet gFacesSet(globalIncludeFaces);
        if(Pstream::myProcNo()!=I and procInlcudeFaces[I].size()>0)
        {
            if(!gFacesSet.found(procInlcudeFaces[I][0]))
            {
                globalIncludeFaces.append(procInlcudeFaces[I]);
            }
        }
    }

    globalIncludeFacesPtr_ =  new labelList(globalIncludeFaces);

    // find globalIncludeEdges, and globalIncludePoints
    labelHashSet globalIncludeEdgesSet;
    labelHashSet globalIncludePointsSet;

    forAll(globalIncludeFaces,I)
    {
        label faceID =  globalIncludeFaces[I];
        const labelList& faceEdges = surf.faceEdges()[faceID];

        forAll(faceEdges,II)
        {
            label edgeID = faceEdges[II];
            if(!globalIncludeEdgesSet.found(edgeID))
            {
                globalIncludeEdgesSet.insert(edgeID);
            }
        }

        const face& myFaces = surf.localFaces()[faceID];

        forAll(myFaces,II)
        {
            label pointID = myFaces[II];
            // use local point

            if(!globalIncludePointsSet.found(pointID))
            {
                globalIncludePointsSet.insert(pointID);
            }
        }
    }

    globalIncludeEdgesPtr_ =  new labelList(globalIncludeEdgesSet.toc());
    globalIncludePointsPtr_ =  new labelList(globalIncludePointsSet.toc());

    // find globalBoundaryEdges and globalBoundaryEdgesPatchID
    labelHashSet globalBoundaryEdgesSet;
    labelList globalBoundaryEdges(surf.nEdges());
    labelListList globalBoundaryEdgesPatchID(surf.nEdges());
    labelListList globalBoundaryEdgesPatchIDProc(surf.nEdges());
    label nBoundaryEdges = 0;

    forAll(globalBoundaryFaces,I)
    {
        label faceID =  globalBoundaryFaces[I];

        const labelList& faceEdges = surf.faceEdges()[faceID];
        forAll(faceEdges,II)
        {
            label edgeID = faceEdges[II];
            if(globalIncludeEdgesSet.found(edgeID) and !globalBoundaryEdgesSet.found(edgeID))
            {
                globalBoundaryEdgesSet.insert(edgeID);
                globalBoundaryEdges[nBoundaryEdges]=edgeID;
                globalBoundaryEdgesPatchID[nBoundaryEdges]=globalBoundaryFacesPatchID[I];
                globalBoundaryEdgesPatchIDProc[nBoundaryEdges]=globalBoundaryFacesPatchIDProc[I];
                nBoundaryEdges++;
            }
        }
    }

    globalBoundaryEdges.setSize(nBoundaryEdges);
    globalBoundaryEdgesPatchID.setSize(nBoundaryEdges);
    globalBoundaryEdgesPatchIDProc.setSize(nBoundaryEdges);
    globalBoundaryEdgesPtr_ = new labelList(globalBoundaryEdges);
    globalBoundaryEdgesPatchIDPtr_ = new labelListList(globalBoundaryEdgesPatchID);
    globalBoundaryEdgesPatchIDProcPtr_ = new labelListList(globalBoundaryEdgesPatchIDProc);

}












void Foam::cutTriSurfaceMesh::makeGlobal
(
    const polyMesh& procMesh,
    const PtrList<triSurface>& addSurfs

) const
{
    // this is used for mpi
    // basic idea is to change boundary faces with processorPolyPatch to include faces
    // global information are made from globalBoundaryFaces

    if (!Pstream::parRun())
    {
        globalIncludeFacesPtr_ = new labelList(includeFaces());
        globalIncludeEdgesPtr_ = new labelList(includeEdges());
        globalIncludePointsPtr_ = new labelList(includePoints());
        globalBoundaryFacesPtr_ = new labelList(boundaryFaces());
        globalBoundaryFacesPatchIDPtr_ = new labelListList(boundaryFacesPatchID());
        labelListList globalBoundaryFacesPatchIDProc(boundaryFacesPatchID().size());
               
        forAll(globalBoundaryFacesPatchIDProc,I)
        {
            globalBoundaryFacesPatchIDProc[I].setSize(boundaryFacesPatchID()[I].size(),0);
        }
        globalBoundaryFacesPatchIDProcPtr_ = new labelListList(globalBoundaryFacesPatchIDProc);

        globalBoundaryEdgesPtr_ = new labelList(boundaryEdges());
        globalBoundaryEdgesPatchIDPtr_ = new labelListList(boundaryEdgesPatchID());
        labelListList globalBoundaryEdgesPatchIDProc(boundaryEdgesPatchID().size());
        forAll(globalBoundaryEdgesPatchIDProc,I)
        {
            globalBoundaryEdgesPatchIDProc[I].setSize(boundaryEdgesPatchID()[I].size(),0);
        }
        globalBoundaryEdgesPatchIDProcPtr_ = new labelListList(globalBoundaryEdgesPatchIDProc);

        return;
    }
    const labelList& boundaryFaces = this->boundaryFaces();
    labelListList& boundaryFacesPatchID = *boundaryFacesPatchIDPtr_;

    // check processor interface for each procMesh
    labelHashSet procPatchListSet;

    // replace processorCyclicPolyPatchID to its referPatchID
    labelList procPatchIDReplace(procMesh.boundaryMesh().size());

    forAll(procMesh.boundaryMesh(),patchID)
    {
        procPatchIDReplace[patchID]=patchID;
        if
        (
            isA<processorCyclicPolyPatch>
            (
                procMesh.boundaryMesh()[patchID]
            )
        )
        {
            const processorCyclicPolyPatch& pcPP =
                        refCast<const processorCyclicPolyPatch>(procMesh.boundaryMesh()[patchID]);
            procPatchIDReplace[patchID]=pcPP.referPatchID();
        }
        else if
        (
            isA<processorPolyPatch>
            (
                procMesh.boundaryMesh()[patchID]
            )
        )
        {
            procPatchListSet.insert(patchID);
        }
    }

    DynamicList<label> procBoundaryFaces;
    DynamicList<labelList> procBoundaryFacesPatchID;
    

    forAll(boundaryFaces,I)
    {
        label faceID = boundaryFaces[I];
        labelList patchIDs = boundaryFacesPatchID[I];
        bool unChanged = false;
        labelHashSet newPatchIDSet;
        forAll(patchIDs, I)
        { 
            if (patchIDs[I] < procMesh.boundaryMesh().size())
            {
                patchIDs[I]=procPatchIDReplace[patchIDs[I]];
                if(!procPatchListSet.found(patchIDs[I]))
                {
                    newPatchIDSet.insert(patchIDs[I]);
                    unChanged = true;
                }
            }
            else
            {
                newPatchIDSet.insert(patchIDs[I]);
                unChanged = true;

            }
        }
        
        if(unChanged)
        {
            procBoundaryFaces.append(faceID);
            procBoundaryFacesPatchID.append(newPatchIDSet.toc());
           
        }
    }
    
    procBoundaryFaces.shrink();
    procBoundaryFacesPatchID.shrink();

    


    // transfer to each processor, not only master
    List<labelList> procBoundaryFacesList(Pstream::nProcs());
    List<labelListList> procBoundaryFacesPatchIDList(Pstream::nProcs());
 

    procBoundaryFacesList[Pstream::myProcNo()]=procBoundaryFaces;
    procBoundaryFacesPatchIDList[Pstream::myProcNo()]=procBoundaryFacesPatchID;



    Pstream::gatherList(procBoundaryFacesList);
    Pstream::scatterList(procBoundaryFacesList);

    Pstream::gatherList(procBoundaryFacesPatchIDList);
    Pstream::scatterList(procBoundaryFacesPatchIDList);


     DynamicList<label> globalBoundaryFaces;
    DynamicList<labelList> globalBoundaryFacesPatchID;
    DynamicList<labelList> globalBoundaryFacesPatchIDProc;
   

    HashTable<label,label> globalBoundaryFacesSet;
    label globalBoundaryFacesCounts=0;

    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        forAll(procBoundaryFacesList[procI],I)
        {
            label faceID = procBoundaryFacesList[procI][I];
            labelList PatchIDs = procBoundaryFacesPatchIDList[procI][I];
           
            // make sure remove duplicate faceID
            if(!globalBoundaryFacesSet.found(faceID))
            {
                globalBoundaryFacesSet.insert(faceID,globalBoundaryFacesCounts);
                globalBoundaryFacesCounts++;
                globalBoundaryFaces.append(faceID);
                globalBoundaryFacesPatchID.append(PatchIDs);
                labelList PatchIDsProc(PatchIDs.size(),procI);
                globalBoundaryFacesPatchIDProc.append(PatchIDsProc);
               
            }
            else // add duplicate face patch IDs and procs
            {
                label newID = globalBoundaryFacesSet.find(faceID)();

                DynamicList<label> tempDList;
                tempDList.append(globalBoundaryFacesPatchID[newID]);
                tempDList.append(PatchIDs);
                tempDList.shrink();
                globalBoundaryFacesPatchID[newID]=tempDList;
                labelList PatchIDsProc(PatchIDs.size(),procI);
                DynamicList<label> tempDList1;
                
                tempDList1.append(globalBoundaryFacesPatchIDProc[newID]);
                tempDList1.append(PatchIDsProc);
                tempDList1.shrink();
                globalBoundaryFacesPatchIDProc[newID]=tempDList1;
            }
        }
    }
    globalBoundaryFaces.shrink();
    globalBoundaryFacesPatchID.shrink();
    globalBoundaryFacesPatchIDProc.shrink();
    globalBoundaryFacesPtr_ = new labelList(globalBoundaryFaces);
    globalBoundaryFacesPatchIDPtr_ = new labelListList(globalBoundaryFacesPatchID);
    globalBoundaryFacesPatchIDProcPtr_ = new labelListList(globalBoundaryFacesPatchIDProc);

   
  
    // the following are the same as in makeBoundaryFacesAndEdges

    // find globalIncludeFaces
    labelHashSet globalIncludeFacesSet;

    label startID=0;
    List<labelHashSet> globalIncludeFacesSetList;
    label nFace=0;
    const triSurface& surf(*this);

    while(nFace<(surf.size()-globalBoundaryFaces.size()))
    {    
        labelHashSet nbFaceIDSet;
        labelList old_nbFaceID;
        old_nbFaceID.append(startID);
        labelHashSet includeFacesLocalSet;
        const labelListList& faceFaces = surf.faceFaces();

        while(old_nbFaceID.size()>0)
        {
            nbFaceIDSet.clear();
            forAll(old_nbFaceID, faceI)
            {
                startID=old_nbFaceID[faceI];
                forAll(faceFaces[startID], faceII)
                {
                    label faceID = faceFaces[startID][faceII];
                    if(!globalIncludeFacesSet.found(faceID) and !globalBoundaryFacesSet.found(faceID) )
                    {
                        nFace++;
                        includeFacesLocalSet.insert(faceID);
                        globalIncludeFacesSet.insert(faceID);
                        nbFaceIDSet.insert(faceID);
                    }
                }
            }
            old_nbFaceID=nbFaceIDSet.toc();
        }

        globalIncludeFacesSetList.append(includeFacesLocalSet);
        forAll(surf,faceID)
        {
            if(!globalIncludeFacesSet.found(faceID) and !globalBoundaryFacesSet.found(faceID) )
            {
                startID=faceID;
                if(includeFacesLocalSet.size()<SMALL)
                {
                        globalIncludeFacesSetList[globalIncludeFacesSetList.size()-1].insert(faceID);
                        globalIncludeFacesSet.insert(faceID);
                        nFace++;
                }
                break;
            }
        }
    }
    globalIncludeFacesSet.clear();
    labelList globalIncludeFaces;

    forAll(globalIncludeFacesSetList,I)
    {
        labelHashSet tempSet=globalIncludeFacesSetList[I];
        label faceID = tempSet.toc()[0];
        scalar indicator1 = -1;
        scalar indicator2 = -1;
        const face& myFaces = surf[faceID];
         
        forAll(addSurfs, surfID)
        {
        
            triSurface addSurf(addSurfs[surfID]);
            
            triSurfaceSearch queryAddSurf(addSurf);
         
            forAll(myFaces,pointI)
            {
            
                vector pt = surf.points()[myFaces[pointI]];
                bool faceInStl = ifInStl(pt,queryAddSurf);
                if(!faceInStl) 
                {
                    indicator1=1;
                    
                }
               
            }
        }
        forAll(myFaces,pointI)
        {

            if(procMesh.findCell(surf.points()[myFaces[pointI]])>-1 )
            {
                indicator2=1;
            }
           
        }
        
        
        

        if(indicator1>0 and indicator2>0)
        {
            globalIncludeFaces.append(tempSet.toc());
        }
    }

    labelListList procInlcudeFaces(Pstream::nProcs());
    procInlcudeFaces[Pstream::myProcNo()]=globalIncludeFaces;
    Pstream::gatherList(procInlcudeFaces);
    Pstream::scatterList(procInlcudeFaces);
    forAll(procInlcudeFaces,I)
    {
        labelHashSet gFacesSet(globalIncludeFaces);
        if(Pstream::myProcNo()!=I and procInlcudeFaces[I].size()>0)
        {
            if(!gFacesSet.found(procInlcudeFaces[I][0]))
            {
                globalIncludeFaces.append(procInlcudeFaces[I]);
            }
        }
    }


    
    
    globalIncludeFacesPtr_ =  new labelList(globalIncludeFaces);

    // find globalIncludeEdges, and globalIncludePoints
    labelHashSet globalIncludeEdgesSet;
    labelHashSet globalIncludePointsSet;

    forAll(globalIncludeFaces,I)
    {
        label faceID =  globalIncludeFaces[I];
        const labelList& faceEdges = surf.faceEdges()[faceID];

        forAll(faceEdges,II)
        {
            label edgeID = faceEdges[II];
            if(!globalIncludeEdgesSet.found(edgeID))
            {
                globalIncludeEdgesSet.insert(edgeID);
            }
        }

        const face& myFaces = surf.localFaces()[faceID];

        forAll(myFaces,II)
        {
            label pointID = myFaces[II];
            // use local point
            //label localPointID = surf.whichPoint(pointID);
            if(!globalIncludePointsSet.found(pointID))
            {
                globalIncludePointsSet.insert(pointID);
            }
        }
    }

    globalIncludeEdgesPtr_ =  new labelList(globalIncludeEdgesSet.toc());
    globalIncludePointsPtr_ =  new labelList(globalIncludePointsSet.toc());

    // find globalBoundaryEdges and globalBoundaryEdgesPatchID
    labelHashSet globalBoundaryEdgesSet;
    labelList globalBoundaryEdges(surf.nEdges());
    labelListList globalBoundaryEdgesPatchID(surf.nEdges());
    labelListList globalBoundaryEdgesPatchIDProc(surf.nEdges());

    label nBoundaryEdges = 0;

    forAll(globalBoundaryFaces,I)
    {
        label faceID =  globalBoundaryFaces[I];

        const labelList& faceEdges = surf.faceEdges()[faceID];
        forAll(faceEdges,II)
        {
            label edgeID = faceEdges[II];
            if(globalIncludeEdgesSet.found(edgeID) and !globalBoundaryEdgesSet.found(edgeID))
            {
                globalBoundaryEdgesSet.insert(edgeID);
                globalBoundaryEdges[nBoundaryEdges]=edgeID;
                globalBoundaryEdgesPatchID[nBoundaryEdges]=globalBoundaryFacesPatchID[I];
                globalBoundaryEdgesPatchIDProc[nBoundaryEdges]=globalBoundaryFacesPatchIDProc[I];
           
                nBoundaryEdges++;
            }
        }
    }

    globalBoundaryEdges.setSize(nBoundaryEdges);
    globalBoundaryEdgesPatchID.setSize(nBoundaryEdges);
    globalBoundaryEdgesPatchIDProc.setSize(nBoundaryEdges);
    globalBoundaryEdgesPtr_ = new labelList(globalBoundaryEdges);
    globalBoundaryEdgesPatchIDPtr_ = new labelListList(globalBoundaryEdgesPatchID);
    globalBoundaryEdgesPatchIDProcPtr_ = new labelListList(globalBoundaryEdgesPatchIDProc);


}




















void Foam::cutTriSurfaceMesh::makeOutsidePointsStencil() const
{
    if(debug==2)
    {
        Pout<<"makeOutsidePointsStencil"<<endl;
    }
    DynamicList<label> outsidePoints; // local
    DynamicList<labelList> outsidePointsStencil; // local
    DynamicList<scalarList> outsidePointsStencilWeight; // local
    labelHashSet includePointsSet(globalIncludePoints()); // local

    const labelList& gBE=globalBoundaryEdges();
    labelHashSet boundaryPointsSet;
    forAll(gBE,I)
    {
        label edgeID=gBE[I];
        edge e = this->edges()[edgeID];
        if(!boundaryPointsSet.found(e[0]))boundaryPointsSet.insert(e[0]);
        if(!boundaryPointsSet.found(e[1]))boundaryPointsSet.insert(e[1]);
    }

    labelHashSet outsidePointsSet;


    // new algorithm
    labelHashSet visitedSet(includePointsSet);
    //labelHashSet oldVisitSet(visitedSet);
    // find stencil for boundaryPoints at first, the same as before
    forAllIter(labelHashSet, boundaryPointsSet, iter)
    {
        label pointID = iter.key();
        outsidePointsSet.insert(pointID);
        outsidePoints.append(pointID);
        labelList stencil =
            findNearestBoundaryPoints
            (
                pointID, // local
                includePointsSet //local
            );
        outsidePointsStencil.append(stencil);
        scalarField stencilWeight(stencil.size(),0);
        // make stencil weights
        scalar sumW=0;
        forAll(stencilWeight,I)
        {
            stencilWeight[I]=mag(this->localPoints()[pointID]-this->localPoints()[stencil[I]]);
            stencilWeight[I] = stencilWeight[I]*stencilWeight[I];
            sumW= sumW+stencilWeight[I];
        }
        stencilWeight=stencilWeight/(sumW+SMALL);
        outsidePointsStencilWeight.append(stencilWeight);
        if(!visitedSet.found(pointID)) visitedSet.insert(pointID);
    }
    // find neighbours of boundary points until all visited
    const labelListList& pointEdges = this->pointEdges();
    labelHashSet nextNeighbourListSet;
    forAllIter(labelHashSet, boundaryPointsSet, iter)
    {
        label pointID = iter.key();
        forAll(pointEdges[pointID],II)
        {
            label edgeID = pointEdges[pointID][II];
            edge e = this->edges()[edgeID];
            label otherPointID=e.otherVertex(pointID);
            if(!visitedSet.found(otherPointID) and !nextNeighbourListSet.found(otherPointID))
            {
                nextNeighbourListSet.insert(otherPointID);
            }
        }
    }
    while(nextNeighbourListSet.size()>0)
    {

        forAllIter(labelHashSet, nextNeighbourListSet, iter)
        {
            label pointID = iter.key();
            outsidePointsSet.insert(pointID);
            outsidePoints.append(pointID);
            labelList stencil =
                findNearestBoundaryPoints
                (
                    pointID, // local
                    visitedSet //local
                );
            outsidePointsStencil.append(stencil);
            scalarField stencilWeight(stencil.size(),0);
            // make stencil weights
            scalar sumW=0;
            forAll(stencilWeight,I)
            {
                stencilWeight[I]=mag(this->localPoints()[pointID]-this->localPoints()[stencil[I]]);
                stencilWeight[I] = stencilWeight[I]*stencilWeight[I];
                sumW= sumW+stencilWeight[I];
            }
            stencilWeight=stencilWeight/(sumW+SMALL);
            outsidePointsStencilWeight.append(stencilWeight);
            if(!visitedSet.found(pointID)) visitedSet.insert(pointID);
        }

        labelHashSet oldNLST(nextNeighbourListSet);
        nextNeighbourListSet.clear();
        forAllIter(labelHashSet, oldNLST, iter)
        {
            label pointID = iter.key();
            forAll(pointEdges[pointID],II)
            {
                label edgeID = pointEdges[pointID][II];
                edge e = this->edges()[edgeID];
                label otherPointID=e.otherVertex(pointID);
                if(!visitedSet.found(otherPointID) and !nextNeighbourListSet.found(otherPointID))
                {
                    nextNeighbourListSet.insert(otherPointID);
                }
            }
        }
    }
    outsidePoints.shrink();
    outsidePointsStencil.shrink();
    outsidePointsStencilWeight.shrink();
    outsidePointsPtr_ = new labelList (outsidePoints);
    outsidePointsStencilPtr_ = new labelListList (outsidePointsStencil);
    outsidePointsStencilWeightPtr_ = new scalarListList (outsidePointsStencilWeight);
}

Foam::labelList Foam::cutTriSurfaceMesh::findNearestBoundaryPoints
(
    const label& pointID, //local
    const labelHashSet& includePointsSet //local
) const
{
    const labelListList& pointEdges = this->pointEdges();
    labelHashSet visitedSet;
    label nextPointID = pointID;

    labelHashSet stencilSet;

    labelHashSet nextNeighbourListSet;
    visitedSet.insert(pointID);

    forAll(pointEdges[nextPointID],II)
    {
        label edgeID = pointEdges[nextPointID][II];
        edge e = this->edges()[edgeID];

        label otherPointID=e.otherVertex(nextPointID);

        if(!visitedSet.found(otherPointID) and
            !nextNeighbourListSet.found(otherPointID))
        {
            nextNeighbourListSet.insert(otherPointID);
        }
    }

    labelList neighbourList = nextNeighbourListSet.toc();
    while(visitedSet.size()<pointEdges.size() and stencilSet.size()<5)
    {
        nextNeighbourListSet.clear();
        forAll(neighbourList,I)
        {
            nextPointID=neighbourList[I];

            if(!visitedSet.found(nextPointID))
            {
                visitedSet.insert(nextPointID);

                if(includePointsSet.found(nextPointID) and !stencilSet.found(nextPointID))
                {
                    stencilSet.insert(nextPointID);
                }

                forAll(pointEdges[nextPointID],II)
                {
                    label edgeID = pointEdges[nextPointID][II];
                    edge e = this->edges()[edgeID];
                    label otherPointID=e.otherVertex(nextPointID);

                    if(!visitedSet.found(otherPointID) and
                        !nextNeighbourListSet.found(otherPointID))
                    {
                        nextNeighbourListSet.insert(otherPointID);
                    }
                }
            }
        }
        neighbourList=nextNeighbourListSet.toc();
    }

    return stencilSet.toc();
}

template<class T>
void Foam::cutTriSurfaceMesh::sendToMaster
(
    List<T>& Values
) const
{
    if (Pstream::parRun())
    {
        if (Pstream::myProcNo() != Pstream::masterNo())
        {
            // Slave sends data to master
            // Do not send empty list
           
            // create the stream to send to master
            OPstream oStream
            (
                Pstream::commsTypes::blocking,
                  Pstream::masterNo()
            );
               oStream << Values[Pstream::myProcNo()];
            
            
        }
        else
        {
            // Master receives data from slaves
            for(label procI = 1; procI < Pstream::nProcs(); procI++)
            {
                // Do not receive empty list
               
                 // create the input stream from processor procI
                 IPstream iStream(Pstream::commsTypes::blocking, procI);
                 iStream >> Values[procI];
                

            }

        }
    }
}
void Foam::cutTriSurfaceMesh::makeBoundaryFacesAndEdges
(
    const polyMesh& mesh
) const
{
    if(debug==2)
    {
        Info<<"makeBoundaryFacesAndEdges"<<endl;
    }
    // make triSurface for each patch
    const double Oldtime1=mesh.time().elapsedCpuTime();

    const polyBoundaryMesh& patches= mesh.boundaryMesh();

    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();

    PtrList<triSurface> patchSurfPtrList(patches.size());

    triFaceList totalTriFaces;

    forAll(patches,patchI)
    {
        const polyPatch& patch = patches[patchI];

        DynamicList<face> dynNewTriFaces(patch.size()*3);

        for (label localFaceID = 0; localFaceID < patch.size(); localFaceID++)
        {
            label globalFaceID = patch.start()+localFaceID;
            face oldFace = faces[globalFaceID];
            oldFace.triangles(points,dynNewTriFaces);
        }
        triFaceList newTriFaces(dynNewTriFaces.shrink());
        totalTriFaces.append(newTriFaces);

        triSurface patchSurf(newTriFaces,points);

        boolList includeMap(patchSurf.size(), true);

        labelList pointMap;
        labelList faceMap;

        patchSurfPtrList.set
        (
            patchI,
            new triSurface
            (
                patchSurf.subsetMesh
                (
                    includeMap,
                    pointMap,
                    faceMap
                )
            )
        );
    }

    const triSurface& surf(*this);

    triSurfaceSearch querySurf(surf);
    
    const double Oldtime2=mesh.time().elapsedCpuTime();
   
     
    // find faces cut by each patch
    List<DynamicList<label>> cutFacesList(surf.size());

    forAll(mesh.boundaryMesh(),patchI)
    {
        triSurface patchSurf = patchSurfPtrList[patchI];

        triSurfaceSearch queryPatchSurf(patchSurf);

        labelList cutFaces = markFaces(querySurf,queryPatchSurf);

        forAll(cutFaces,faceI)
        {
            cutFacesList[cutFaces[faceI]].append(patchI);
        }
    }


    // find boundaryFaces
    labelList boundaryFaces(cutFacesList.size());
    labelListList boundaryFacesPatchID(cutFacesList.size());
    label nBoundaryFaces =0;
    forAll(cutFacesList,faceI)
    {
        cutFacesList[faceI].shrink();
        if(cutFacesList[faceI].size()>0)
        {
            boundaryFaces[nBoundaryFaces]=faceI;
            boundaryFacesPatchID[nBoundaryFaces]=cutFacesList[faceI];
            nBoundaryFaces++;
        }
    }
    boundaryFaces.setSize(nBoundaryFaces);
    boundaryFacesPatchID.setSize(nBoundaryFaces);
    boundaryFacesPtr_ = new labelList(boundaryFaces);
    boundaryFacesPatchIDPtr_ = new labelListList(boundaryFacesPatchID);
    const double Oldtime3=mesh.time().elapsedCpuTime();
   
    // find includeFaces
    labelHashSet boundaryFacesSet(boundaryFaces);
    labelHashSet includeFacesSet;

    label startID=0;
    List<labelHashSet> includeFacesSetList;
    label nFace=0;

    while(nFace<(surf.size()-boundaryFaces.size()))
    {
        labelHashSet nbFaceIDSet;
        labelList old_nbFaceID;
        old_nbFaceID.append(startID);
        labelHashSet includeFacesLocalSet;
        const labelListList& faceFaces = surf.faceFaces();

        while(old_nbFaceID.size()>0)
        {
            nbFaceIDSet.clear();
            forAll(old_nbFaceID, faceI)
            {
                startID=old_nbFaceID[faceI];
                forAll(faceFaces[startID], faceII)
                {
                    label faceID = faceFaces[startID][faceII];
                    if(!includeFacesSet.found(faceID) and !boundaryFacesSet.found(faceID) )
                    {
                        nFace++;
                        includeFacesLocalSet.insert(faceID);
                        includeFacesSet.insert(faceID);
                        nbFaceIDSet.insert(faceID);
                    }
                }
            }
            old_nbFaceID=nbFaceIDSet.toc();
        }

        includeFacesSetList.append(includeFacesLocalSet);
        forAll(surf,faceID)
        {
            if(!includeFacesSet.found(faceID) and !boundaryFacesSet.found(faceID) )
            {
                startID=faceID;
                if(includeFacesLocalSet.size()<SMALL )
                {
                        includeFacesSetList[includeFacesSetList.size()-1].insert(faceID);
                        includeFacesSet.insert(faceID);
                        nFace++;
                }
                break;
            }
        }
    }
    includeFacesSet.clear();
    labelList includeFaces;

    forAll(includeFacesSetList,I)
    {
        labelHashSet tempSet=includeFacesSetList[I];
        label faceID = tempSet.toc()[0];
        scalar indicator = 1;

        const face& myFaces = surf[faceID];

        forAll(myFaces,pointI)
        {
            if(mesh.findCell(surf.points()[myFaces[pointI]])<0)
            {
                indicator=-1;
            }
        }
        if(indicator>0)
        {
            includeFaces.append(tempSet.toc());
        }
    }

    includeFacesPtr_ =  new labelList(includeFaces);
         const double Oldtime31=mesh.time().elapsedCpuTime();
    



    // find includeEdges, and includePoints
    labelHashSet includeEdgesSet;
    labelHashSet includePointsSet;

    forAll(includeFaces,I)
    {
        label faceID =  includeFaces[I];
        const labelList& faceEdges = surf.faceEdges()[faceID];

        forAll(faceEdges,II)
        {
            label edgeID = faceEdges[II];
            if(!includeEdgesSet.found(edgeID))
            {
                includeEdgesSet.insert(edgeID);
            }
        }

        const face& myFaces = surf.localFaces()[faceID];

        forAll(myFaces,II)
        {
            label pointID = myFaces[II];
            // use local point
            //label localPointID = surf.whichPoint(pointID);
            if(!includePointsSet.found(pointID))
            {
                includePointsSet.insert(pointID);
            }
        }
    }

    includeEdgesPtr_ =  new labelList(includeEdgesSet.toc());
    includePointsPtr_ =  new labelList(includePointsSet.toc());
     const double Oldtime4=mesh.time().elapsedCpuTime();
    
    // find boundaryEdges and boundaryEdgesPatchID
    labelHashSet boundaryEdgesSet;
    labelList boundaryEdges(surf.nEdges());
    labelListList boundaryEdgesPatchID(surf.nEdges());
    label nBoundaryEdges = 0;

    labelList faceEdgesMap(boundaryFaces.size(),-1);
    forAll(boundaryFaces,I)
    {
        label faceID =  boundaryFaces[I];

        const labelList& faceEdges = surf.faceEdges()[faceID];
        forAll(faceEdges,II)
        {
            label edgeID = faceEdges[II];
            if(includeEdgesSet.found(edgeID) and !boundaryEdgesSet.found(edgeID))
            {
                boundaryEdgesSet.insert(edgeID);
                boundaryEdges[nBoundaryEdges]=edgeID;
                boundaryEdgesPatchID[nBoundaryEdges]=boundaryFacesPatchID[I];
                faceEdgesMap[I]=nBoundaryEdges;
                nBoundaryEdges++;
            }
        }
    }

    boundaryEdges.setSize(nBoundaryEdges);
    boundaryEdgesPatchID.setSize(nBoundaryEdges);
    boundaryEdgesPtr_ = new labelList(boundaryEdges);
    boundaryEdgesPatchIDPtr_ = new labelListList(boundaryEdgesPatchID);
    const double Oldtime5=mesh.time().elapsedCpuTime();
    if(debug)
    {
        Info<<" make triSurface for each patch "<<Oldtime2-Oldtime1<< " s"<<endl;
        Info<<"find boundaryFaces "<<Oldtime3-Oldtime2<< " s"<<endl;

        Info<<" inernalFace "<<Oldtime31-Oldtime3<< " s"<<endl;
        Info<<" inernalEdge "<<Oldtime4-Oldtime31<< " s"<<endl;

        Info<<" boundaryEdge "<<Oldtime5-Oldtime4<< " s"<<endl;
    }
 
}







void Foam::cutTriSurfaceMesh::makeBoundaryFacesAndEdges
(
    const polyMesh& mesh,
    const PtrList<triSurface>& addSurfs
) const
{
    if(debug==2)
    {
        Info<<"makeBoundaryFacesAndEdges"<<endl;
    }
    // make triSurface for each patch
    const polyBoundaryMesh& patches= mesh.boundaryMesh();

    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    label nAddsurf = addSurfs.size();
    
    PtrList<triSurface> patchSurfPtrList(patches.size()+nAddsurf);

    triFaceList totalTriFaces;

    forAll(patches,patchI)
    {
        const polyPatch& patch = patches[patchI];

        DynamicList<face> dynNewTriFaces(patch.size()*3);

        for (label localFaceID = 0; localFaceID < patch.size(); localFaceID++)
        {
            label globalFaceID = patch.start()+localFaceID;
            face oldFace = faces[globalFaceID];
            oldFace.triangles(points,dynNewTriFaces);
        }
        triFaceList newTriFaces(dynNewTriFaces.shrink());
        totalTriFaces.append(newTriFaces);

        triSurface patchSurf(newTriFaces,points);

        boolList includeMap(patchSurf.size(), true);

        labelList pointMap;
        labelList faceMap;

        patchSurfPtrList.set
        (
            patchI,
            new triSurface
            (
                patchSurf.subsetMesh
                (
                    includeMap,
                    pointMap,
                    faceMap
                )
            )
        );
    }
    
    
    forAll(addSurfs,triSurfI)
    {
        triSurface patchSurf(addSurfs[triSurfI]); 
        
        boolList includeMap(patchSurf.size(), true);

        labelList pointMap;
        labelList faceMap;
        patchSurfPtrList.set
        (
            patches.size()+triSurfI,
            new triSurface
            (
                patchSurf.subsetMesh
                (
                    includeMap,
                    pointMap,
                    faceMap
                )
            )
        );
        Info<<" add object to cut triSurface "<<endl;
    }
    


    const triSurface& surf(*this);

    triSurfaceSearch querySurf(surf);

    // find faces cut by each patch
    List<DynamicList<label>> cutFacesList(surf.size());

    for (label patchI = 0; patchI < patchSurfPtrList.size(); patchI++)
    {
        triSurface patchSurf = patchSurfPtrList[patchI];
    
        triSurfaceSearch queryPatchSurf(patchSurf);

        labelList cutFaces = markFaces(querySurf,queryPatchSurf);
        
        forAll(cutFaces,faceI)
        {
            cutFacesList[cutFaces[faceI]].append(patchI);
        }
        
    }
    
    // find boundaryFaces
    labelList boundaryFaces(cutFacesList.size());
    labelListList boundaryFacesPatchID(cutFacesList.size());
    label nBoundaryFaces =0;
    forAll(cutFacesList,faceI)
    {
        cutFacesList[faceI].shrink();
        if(cutFacesList[faceI].size()>0)
        {
            boundaryFaces[nBoundaryFaces]=faceI;
            boundaryFacesPatchID[nBoundaryFaces]=cutFacesList[faceI];
            nBoundaryFaces++;
        }
    }
    boundaryFaces.setSize(nBoundaryFaces);
    boundaryFacesPatchID.setSize(nBoundaryFaces);
    boundaryFacesPtr_ = new labelList(boundaryFaces);
    boundaryFacesPatchIDPtr_ = new labelListList(boundaryFacesPatchID);

    // find includeFaces
    labelHashSet boundaryFacesSet(boundaryFaces);
    labelHashSet includeFacesSet;

    label startID=0;
    List<labelHashSet> includeFacesSetList;
    label nFace=0;
    
    const labelListList& faceFaces = surf.faceFaces();
    while(nFace<(surf.size()-boundaryFaces.size()))
    {
        labelHashSet nbFaceIDSet;
        labelList old_nbFaceID;
        old_nbFaceID.append(startID);
        labelHashSet includeFacesLocalSet;
        
        while(old_nbFaceID.size()>0)
        {
            nbFaceIDSet.clear();
            forAll(old_nbFaceID, faceI)
            {
                startID=old_nbFaceID[faceI];
                forAll(faceFaces[startID], faceII)
                {
                    label faceID = faceFaces[startID][faceII];
                    if(!includeFacesSet.found(faceID) and !boundaryFacesSet.found(faceID) )
                    {
                        nFace++;
                        includeFacesLocalSet.insert(faceID);
                        includeFacesSet.insert(faceID);
                        nbFaceIDSet.insert(faceID);
                    }
                }
            }
            old_nbFaceID=nbFaceIDSet.toc();
        }

        includeFacesSetList.append(includeFacesLocalSet);
        forAll(surf,faceID)
        {
            if(!includeFacesSet.found(faceID) and !boundaryFacesSet.found(faceID) )
            {
                startID=faceID;
                if(includeFacesLocalSet.size()<SMALL )
                {
                        includeFacesSetList[includeFacesSetList.size()-1].insert(faceID);
                        includeFacesSet.insert(faceID);
                        nFace++;
                }
                break;
            }
        }
    }
    includeFacesSet.clear();
    labelList includeFaces;

    forAll(includeFacesSetList,I)
    {
        labelHashSet tempSet=includeFacesSetList[I];
        label faceID = tempSet.toc()[0];
        scalar indicator = 1;

        const face& myFaces = surf[faceID];
      
        forAll(addSurfs, surfID)
        {
        
            triSurface addSurf(addSurfs[surfID]);
            
            triSurfaceSearch queryAddSurf(addSurf);
            
            forAll(myFaces,pointI)
            {
            
                vector pt = surf.points()[myFaces[pointI]];
                bool faceInStl = ifInStl(pt,queryAddSurf);
                if(faceInStl) 
                {
                    indicator=-1;
                  
                }
            }
          
        }
        
      
        forAll(myFaces,pointI)
        {

            if(mesh.findCell(surf.points()[myFaces[pointI]])<0 )
            {
                indicator=-1;
               
            }
        }
        
        
        
        
        if(indicator>0)
        {
            includeFaces.append(tempSet.toc());
        }
       
    }
    
    
  
    includeFacesPtr_ =  new labelList(includeFaces);
    // find includeEdges, and includePoints
    labelHashSet includeEdgesSet;
    labelHashSet includePointsSet;

    forAll(includeFaces,I)
    {
        label faceID =  includeFaces[I];
        const labelList& faceEdges = surf.faceEdges()[faceID];

        forAll(faceEdges,II)
        {
            label edgeID = faceEdges[II];
            if(!includeEdgesSet.found(edgeID))
            {
                includeEdgesSet.insert(edgeID);
            }
        }

        const face& myFaces = surf.localFaces()[faceID];

        forAll(myFaces,II)
        {
            label pointID = myFaces[II];
            // use local point
  
            if(!includePointsSet.found(pointID))
            {
                includePointsSet.insert(pointID);
            }
        }
    }

    includeEdgesPtr_ =  new labelList(includeEdgesSet.toc());
    includePointsPtr_ =  new labelList(includePointsSet.toc());

    // find boundaryEdges and boundaryEdgesPatchID
    labelHashSet boundaryEdgesSet;
    labelList boundaryEdges(surf.nEdges());
    labelListList boundaryEdgesPatchID(surf.nEdges());
    label nBoundaryEdges = 0;

    labelList faceEdgesMap(boundaryFaces.size(),-1);
    forAll(boundaryFaces,I)
    {
        label faceID =  boundaryFaces[I];

        const labelList& faceEdges = surf.faceEdges()[faceID];
        forAll(faceEdges,II)
        {
            label edgeID = faceEdges[II];
            if(includeEdgesSet.found(edgeID) and !boundaryEdgesSet.found(edgeID))
            {
                boundaryEdgesSet.insert(edgeID);
                boundaryEdges[nBoundaryEdges]=edgeID;
                boundaryEdgesPatchID[nBoundaryEdges]=boundaryFacesPatchID[I];
                faceEdgesMap[I]=nBoundaryEdges;
                nBoundaryEdges++;
            }
        }
    }

    boundaryEdges.setSize(nBoundaryEdges);
    boundaryEdgesPatchID.setSize(nBoundaryEdges);
    boundaryEdgesPtr_ = new labelList(boundaryEdges);
    boundaryEdgesPatchIDPtr_ = new labelListList(boundaryEdgesPatchID);

}









Foam::labelList Foam::cutTriSurfaceMesh::markFaces
(
    const triSurfaceSearch& querySurf,
    const triSurfaceSearch& queryPatchSurf
)const
{
    const triSurface& surf(querySurf.surface());

    labelHashSet cutFacesSet;
    label nCutFaces = 0;

    forAll(surf.edges(), edgeI)
    {
        const edge& e = surf.edges()[edgeI];

        const point& p0 = surf.localPoints()[e.start()];
        const point& p1 = surf.localPoints()[e.end()];

        pointIndexHit pHit(queryPatchSurf.tree().findLineAny(p0, p1));

        if (pHit.hit())
        {
            const labelList& myFaces = surf.edgeFaces()[edgeI];

            forAll(myFaces, myFaceI)
            {
                label faceI = myFaces[myFaceI];

                if (!cutFacesSet.found(faceI))
                {
                    cutFacesSet.insert(faceI);
                 
                    nCutFaces++;
                }
            }
        }
     }

    return cutFacesSet.toc();
}




// determine if a point is in stl or not
bool Foam::cutTriSurfaceMesh::ifInStl
(
    const vector& C,
    const triSurfaceSearch& querySurf
) const
{
    bool ifInStl = true;

   
    scalar Delta = querySurf.tree().bb().typDim();
   

    vector span
        (
            50*Delta,
            50*Delta,
            50*Delta
        );

    pointIndexHit pih = querySurf.nearest(C, span);

    
    if (pih.hit())
    {
        point nearestPoint = pih.hitPoint();

        vector nearestPointNormal =
            triSurfaceTools::surfaceNormal
            (
                querySurf.surface(),
                pih.index(),
                pih.hitPoint()
            );

        scalar Indicator = nearestPointNormal & ( nearestPoint - C );

        if (Indicator>0)
        {
            ifInStl = true; // cellID inside
        }
        else
        {
            ifInStl = false; // cellID outside
        }
    }
    else
    { 
        FatalErrorIn
                (
                    "Face in the STL has problem"
                )   << "Can't find nearest triSurface point for point "
                    << C << ", "
                    << "span = " << span
                    << "\nYou could try to increase the search span. "
                    << abort(FatalError);
    }

    return ifInStl;
}











Foam::cutTriSurfaceMesh::~cutTriSurfaceMesh()
{
    clearOut();
}


void Foam::cutTriSurfaceMesh::clearOut()
{
    deleteDemandDrivenData(includeFacesPtr_);
    deleteDemandDrivenData(includeEdgesPtr_);
    deleteDemandDrivenData(includePointsPtr_);
    deleteDemandDrivenData(boundaryFacesPtr_);
    deleteDemandDrivenData(boundaryFacesPatchIDPtr_);
    deleteDemandDrivenData(boundaryFacesNearestPatchFaceInfoPtr_);
    deleteDemandDrivenData(boundaryEdgesPtr_);
    deleteDemandDrivenData(boundaryEdgesPatchIDPtr_);
    deleteDemandDrivenData(boundaryEdgesNearestPatchFaceInfoPtr_);
    deleteDemandDrivenData(outsidePointsPtr_);
     deleteDemandDrivenData(outsidePointsStencilPtr_);
     deleteDemandDrivenData(outsidePointsStencilWeightPtr_);
    deleteDemandDrivenData(globalIncludeFacesPtr_);
    deleteDemandDrivenData(globalIncludeEdgesPtr_);
    deleteDemandDrivenData(globalIncludePointsPtr_);
    deleteDemandDrivenData(globalBoundaryFacesPtr_);
    deleteDemandDrivenData(globalBoundaryFacesPatchIDPtr_);
    deleteDemandDrivenData(globalBoundaryFacesPatchIDProcPtr_);
    deleteDemandDrivenData(globalBoundaryFacesNearestPatchFaceInfoPtr_);
    deleteDemandDrivenData(globalBoundaryEdgesPtr_);
    deleteDemandDrivenData(globalBoundaryEdgesPatchIDPtr_);
    deleteDemandDrivenData(globalBoundaryEdgesPatchIDProcPtr_);
    deleteDemandDrivenData(globalBoundaryEdgesNearestPatchFaceInfoPtr_);
    
    deleteDemandDrivenData(dualPatchPtr_);
    deleteDemandDrivenData(includeDualFacesPtr_);
    deleteDemandDrivenData(includeDualEdgesPtr_);
    deleteDemandDrivenData(includeDualPointsPtr_);
    deleteDemandDrivenData(boundaryDualEdgesPtr_);        
    deleteDemandDrivenData(boundaryDualEdgesPatchIDPtr_);          
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
#include "makeDualPatch.C"

// ************************************************************************* //
