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

not used any more

SourceFiles
    immersedBoundaryFvMeshMapping.C
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::immersedBoundaryFvMesh::makeTriAddressing()const
{
    cutTriSurfaceMeshListPtr_ = new PtrList<cutTriSurfaceMesh>(objectsList().size());
    PtrList<cutTriSurfaceMesh>& cutTriSurfaceMeshList = *cutTriSurfaceMeshListPtr_;

    ibCellsToTriAddrListPtr_ = new PtrList<labelListList>(objectsList().size());
    ibCellsToTriWeightsListPtr_ = new PtrList<scalarListList>(objectsList().size());
    ibCellsToTriMatrixListPtr_ = new PtrList<PtrList<scalarSquareMatrix>>(objectsList().size());
    ghostCellsToTriAddrListPtr_ = new PtrList<labelListList>(objectsList().size());
    ghostCellsToTriWeightsListPtr_ = new PtrList<scalarListList>(objectsList().size());
    ghostCellsToTriMatrixListPtr_ = new PtrList<PtrList<scalarSquareMatrix>>(objectsList().size());

    ibCellsToTriEdgeAddrListPtr_ = new PtrList<labelListList>(objectsList().size());
    ibCellsToTriEdgeWeightsListPtr_ = new PtrList<scalarListList>(objectsList().size());
    ghostCellsToTriEdgeAddrListPtr_ = new PtrList<labelListList>(objectsList().size());
    ghostCellsToTriEdgeWeightsListPtr_ = new PtrList<scalarListList>(objectsList().size());

    proTriFacesInMeshListPtr_ = new PtrList<labelListList>(objectsList().size());

    
    label sedimentObject = 0;
     forAll(objectsList(),objectID)
    {
        if( objectDictList()[objectID].found("sediment")) sedimentObject = objectID;
        
    }
    bool multipleObject_ =  ibProperties().lookupOrDefault<bool>("multipleObject",false);
    forAll(objectsList(),objectID)
    {
        const double Oldtime1=time().elapsedCpuTime();
        if( objectDictList()[objectID].found("sediment") and multipleObject_)
        {
            cutTriSurfaceMeshList.set
            (
                objectID,
                new cutTriSurfaceMesh
                (
                    objectsList()[sedimentObject],
                    *this,
                    addObjectsList()
                )
            );    
        }
       else
        {
        
            cutTriSurfaceMeshList.set
            (
                objectID,
                new cutTriSurfaceMesh
                (
                    objectsList()[objectID],
                    *this
                )
            );
        }    
        
        
        const double Oldtime2=time().elapsedCpuTime();

        // makeTriEdgeAddressing has to be used together with(after) makeTriAddressing
        if(IBtypeList()[objectID] == "classic")
        {
            makeTriAddressing(objectID,false);// for ib cells
            const double Oldtime21=time().elapsedCpuTime();

            const double Oldtime22=time().elapsedCpuTime();
            Info<<"makeTriEdgeAddressing Executation Time = "<<Oldtime22-Oldtime21<< " s"<<endl;
        }
        if(IBtypeList()[objectID] == "ghost-cell" or IBtypeList()[objectID] == "mix")
        {
            makeTriAddressing(objectID,true);// for ghost cells
            makeTriEdgeAddressing(objectID,true);
        }
        const double Oldtime3=time().elapsedCpuTime();
        if(debug)
        {
             Info<<"cutTriSurfaceMeshList Executation Time = "<<Oldtime2-Oldtime1<< " s"<<endl;
             Info<<"MakeTriAddressing Executation Time = "<<Oldtime3-Oldtime2<< " s"<<endl;
        }
    }

}

// For each edge, it includes its two(one) neighbour faces' addressing.
// That's why it requires makeTriAddressing first.
void Foam::immersedBoundaryFvMesh::makeTriEdgeAddressing
(
    const label& objectID,
    bool ifGhost
)const
{
    DynamicList<label> hf; // hit faces
    DynamicList<vector> ibp; // ib/hit points

    if(Pstream::parRun())
    {
        // for mpi, sending hf and ibp to each processor
        List<DynamicList<label>> prochf(Pstream::nProcs());
        List<DynamicList<vector>> proibp(Pstream::nProcs());

        if(ifGhost)
        {
            prochf[Pstream::myProcNo()] = ghostHitFacesList()[objectID];
            proibp[Pstream::myProcNo()] = ghostHitPointsList()[objectID];
        }
        else
        {
            prochf[Pstream::myProcNo()] = ibHitFacesList()[objectID];
            prochf[Pstream::myProcNo()].shrink();
            proibp[Pstream::myProcNo()] = ibHitPointsList()[objectID];
            proibp[Pstream::myProcNo()].shrink();
        }

        Pstream::gatherList(prochf);
        Pstream::scatterList(prochf);
        Pstream::gatherList(proibp);
        Pstream::scatterList(proibp);


        forAll(proibp,procI)
        {
            
            ibp.append(proibp[procI]);
            hf.append(prochf[procI]);
        
        }
        
    }
    else
    {
        if(ifGhost)
        {
            hf.append(ghostHitFacesList()[objectID]);
            ibp.append(ghostHitPointsList()[objectID]);
        }
        else
        {
            hf.append(ibHitFacesList()[objectID]);
            ibp.append(ibHitPointsList()[objectID]);
        }
    }
    hf.shrink();
    ibp.shrink();
    
    const triSurface& triPatch = objectsList()[objectID];
    const edgeList& edges = triPatch.edges();
    const labelListList& edgeFaces = triPatch.edgeFaces();
    const labelListList& faceFaces = triPatch.faceFaces();
    labelListList triAddr;

    if(ifGhost)
    {
        triAddr = ghostCellsToTriAddrList()[objectID];
    }
    else
    {
        triAddr = ibCellsToTriAddrList()[objectID];
    }

    labelListList triEdgeAddr(edges.size());
    scalarListList triEdgeWeights(edges.size());

    const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];

    DynamicList<label> triFacesInMesh;
    triFacesInMesh.append(cTSM.globalIncludeFaces());
    triFacesInMesh.append(cTSM.globalBoundaryFaces());
    triFacesInMesh.shrink();


    // only include edges that belong to includeFaces
    DynamicList<label> triEdgesInMesh;
    triEdgesInMesh.append(cTSM.globalIncludeEdges());

    triEdgesInMesh.shrink();
    labelHashSet triFacesInMeshSet(triFacesInMesh);
    
    // for each include edge, find its neighbour faces
    // and their addressing
    forAll(triEdgesInMesh,I)
    {
        label edgeID = triEdgesInMesh[I];

        DynamicList<label> curDynamicAddr;
        labelHashSet faceNeiSet;
        labelList edgeFace = edgeFaces[edgeID];
        forAll(edgeFace,II)
        {
            const labelList& curFaces = faceFaces[edgeFace[II]];
            forAll(curFaces,III)
            {
                if(!faceNeiSet.found(curFaces[III]))
                {
                    faceNeiSet.insert(curFaces[III]);
                }
            }
        }

        forAllIter(labelHashSet, faceNeiSet, iter)
        {
            label faceID = iter.key();
            if(triFacesInMeshSet.found(faceID))
            {
                curDynamicAddr.append(triAddr[faceID]);
            }
        }
         
        curDynamicAddr.shrink(); 
        labelHashSet curAddrSet;
        forAll(curDynamicAddr,II)
        {
            if(!curAddrSet.found(curDynamicAddr[II])) curAddrSet.insert(curDynamicAddr[II]);
        }

        labelList& curAddr = triEdgeAddr[edgeID];
        curAddr = curAddrSet.toc();

        scalarList& curW = triEdgeWeights[edgeID];
        curW.setSize(curAddr.size());
        scalar sumW = 0.0;

        point edgeCentre = edges[edgeID].centre(triPatch.localPoints());
        scalar edgeMag = edges[edgeID].mag(triPatch.localPoints());
        scalar sigma = max(0.001,edgeMag);
        forAll (curAddr, II)
        {
            vector Length = edgeCentre - ibp[curAddr[II]];
            //Gaussian Smoothing Filter
            curW[II]=   1.0/2.0/constant::mathematical::pi/sqr(sigma)
                        *Foam::exp(-magSqr(Length)/2.0/sqr(sigma))+SMALL;
            sumW += curW[II];
        }
        
        forAll (curW, II)
		{
		    curW[II] /= sumW;
		}

	    }


	    if(ifGhost)
	    {
		ghostCellsToTriEdgeAddrListPtr_->set
		(
		    objectID,
		    new labelListList
		    (
			triEdgeAddr
		    )
		);
		ghostCellsToTriEdgeWeightsListPtr_->set
		(
		    objectID,
		    new scalarListList
		    (
			triEdgeWeights
		    )
		);
	    }
	    else
	    {
		ibCellsToTriEdgeAddrListPtr_->set
		(
		    objectID,
		    new labelListList
		    (
			triEdgeAddr
		    )
		);
		ibCellsToTriEdgeWeightsListPtr_->set
		(
		    objectID,
		    new scalarListList
		    (
			triEdgeWeights
		    )
		);
	}
}

// It calculates addressing between each triFace and hit points.
// If ifGhost is true, it refers to the addressing between
// triFace and hit points based on ghost cell center.
void Foam::immersedBoundaryFvMesh::makeTriAddressing
(
    const label& objectID,
    bool ifGhost
)const
{

    // hit Point to triFace center
    // Get reference to tri patch and hit faces
    const triSurface& triPatch = objectsList()[objectID];
    const vectorField& triCentres = triPatch.faceCentres();

    DynamicList<label> hf; // hit faces
    DynamicList<vector> ibp; // hit/ib points
    
     const double Oldtime1=time().elapsedCpuTime(); 
    
    if(Pstream::parRun())
    {
	    // for mpi, sending hf and ibp to each processor
	    List<DynamicList<label>> prochf(Pstream::nProcs());
	    List<DynamicList<vector>> proibp(Pstream::nProcs());

	    if(ifGhost)
	    {
	        prochf[Pstream::myProcNo()] = ghostHitFacesList()[objectID];
	        proibp[Pstream::myProcNo()] = ghostHitPointsList()[objectID];
	    }
	    else
	    {
	        prochf[Pstream::myProcNo()] = ibHitFacesList()[objectID];
                    prochf[Pstream::myProcNo()].shrink();
	        proibp[Pstream::myProcNo()] = ibHitPointsList()[objectID];
                    proibp[Pstream::myProcNo()].shrink();
	    }

	    Pstream::gatherList(prochf);
	    Pstream::scatterList(prochf);
	    Pstream::gatherList(proibp);
	    Pstream::scatterList(proibp);

	    
        forAll(proibp,procI)
        {
	           
            ibp.append(proibp[procI]);
            hf.append(prochf[procI]);
	        
        }
	    
        }
        else
        {
	    if(ifGhost)
	    {
	        hf.append(ghostHitFacesList()[objectID]);
	        ibp.append(ghostHitPointsList()[objectID]);
	    }
	    else
	    {
	        hf.append(ibHitFacesList()[objectID]);
	        ibp.append(ibHitPointsList()[objectID]);
	    }
    }
    hf.shrink();
    ibp.shrink();  
    const double Oldtime2=time().elapsedCpuTime();
    // mark hit point ID inside each triFace
    List<DynamicList<label>> hitTris(triPatch.size());

    forAll (hf, hfI)
    {
	    hitTris[hf[hfI]].append(hfI);
    }
    labelListList addr(triPatch.size());
    scalarListList w(triPatch.size());
    
    labelListList globalAddr(triPatch.size());
    scalarListList globalW(triPatch.size());

    DynamicList<label> triFacesInMesh;

    const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];


    triFacesInMesh.append(cTSM.includeFaces());
    triFacesInMesh.append(cTSM.boundaryFaces());
    
    triFacesInMesh.shrink();
        
    if(Pstream::parRun())
    {

        List<DynamicList<label>> proTri(Pstream::nProcs());
        proTri[Pstream::myProcNo()]=triFacesInMesh;
        proTri[Pstream::myProcNo()].shrink();
        Pstream::gatherList(proTri);   
        Pstream::scatterList(proTri);

           
        proTriFacesInMeshListPtr_->set
        (
         objectID,
             new labelListList
             (
              	proTri
             )
        );
    }

       
    const labelListList& faceFaces = triPatch.faceFaces();

    if(debug)
    {
	    Info<< "triangles: " << triPatch.size() << " triFacesInMesh: " << returnReduce(triFacesInMesh.size(), maxOp<label>()) << " hit: " << returnReduce(hf.size(), maxOp<label>()) << endl;
    }
   
     const double Oldtime3=time().elapsedCpuTime();  
   
    scalar triInterpSize_ = ibProperties().lookupOrDefault<scalar>("triInterpSize",5);
    triInterpSize_ = objectDictList()[objectID].lookupOrDefault<scalar>("triInterpSize",triInterpSize_);

    forAll (triFacesInMesh, tfimI)
    {

	    const label triI = triFacesInMesh[tfimI];

	    scalar sumW = 0;

	    // No direct hit.  Start a neighbourhood search

	    // Record already visited faces
	    labelHashSet visited;

	    // Collect new faces to visit
	    labelHashSet nextToVisit;

	    // Collect new faces to visit
	    //SLList<label> nextToVisit;

	    // Collect IB points for interpolation
	    labelHashSet ibPointsToUse(hitTris[triI].shrink());

	    // Initialise with the original tri
	    nextToVisit.insert(triI);

	    // Loop its neighbours and neighbours' neighbours
	    // until hit points number is larger than triInterpSize_
	    do
	    {
	        const labelList NTV = nextToVisit.toc();
	        nextToVisit.clear();
	        forAll(NTV,I)
	        {
		    const label& next_triI=NTV[I];
		    const labelList& hitPointsInNextTriI=hitTris[next_triI];
		    if (!visited.found(next_triI))
		    {
		        if (hitPointsInNextTriI.size() > 0)
		        {
			    forAll (hitPointsInNextTriI, ptI)
			    {
			        const label& hitPtI=hitPointsInNextTriI[ptI];
			        if(!ibPointsToUse.found(hitPtI))
			        {
				    ibPointsToUse.insert(hitPtI);
			        }
			    }
		        }
		        visited.insert(hitPointsInNextTriI);
		    }
		    const labelList& ff=faceFaces[next_triI];
		    forAll(ff,ffI)
		    {
		        const label ff_I=ff[ffI];
		        if (!visited.found(ff_I) and !nextToVisit.found(ff_I))
		        {
			    nextToVisit.insert(ff_I);
		        }
		    }
	        }
	    } while
	    (
	        ibPointsToUse.size() < triInterpSize_
	     && !nextToVisit.empty()
	    );
	
	    // Found neighbourhood: collect addressing and weights
	    addr[triI] = ibPointsToUse.toc();
	    w[triI].setSize(addr[triI].size());

	    labelList& curAddr = addr[triI];
	    scalarList& curW = w[triI];

	    vector curTriCentre = triCentres[triI];
	    face iFace = triPatch[triI];

	    forAll (curAddr, ibI)
	    {
	        vector Length = curTriCentre - ibp[curAddr[ibI]];
	        
	        curW[ibI] = 1.0/(mag(Length)+SMALL);

	        curW[ibI] = sqrt(curW[ibI]);
	        sumW += curW[ibI];
	    }     
	    // Divide weights by sum distance
	    forAll (curW, ibI)
	    {
	        curW[ibI] /= sumW;
	    }
	    
    }
     const double Oldtime4=time().elapsedCpuTime();
    if(debug)
    {
        Info<<"makeTri Executation Time 1 = "<<Oldtime2-Oldtime1<< " s"<<endl;
        Info<<"makeTri Executation Time 2 = "<<Oldtime3-Oldtime2<< " s"<<endl;
        Info<<"makeTri Executation Time 3 = "<<Oldtime4-Oldtime3<< " s"<<endl;
    }
     


    globalAddr = addr;
    globalW = w; 
   

    
    if(ifGhost)
    {
        ghostCellsToTriAddrListPtr_->set
        (
            objectID,
            new labelListList
            (
                addr
            )
        );
        ghostCellsToTriWeightsListPtr_->set
        (
            objectID,
            new scalarListList
            (
                w
            )
        );

    }
    else
    {
        ibCellsToTriAddrListPtr_->set
        (
            objectID,
            new labelListList
            (
                globalAddr
            )
        );
        ibCellsToTriWeightsListPtr_->set
        (
            objectID,
            new scalarListList
            (
                globalW
            )
        );

    }

}



// global will be done totally
template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::mapFromIBHitToTriFace
(
    const Field<Type>& ibValues,
    label objectID
)const
{
    if (ibValues.size() != ibCellsList()[objectID].size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::Field<Type>\n"
            "immersedBoundaryFvMesh::mapFromIBHitToTriFace\n"
            "(\n"
            "    const Field<Type>& ibValues,\n"
            "    label objectID\n"
            ") const"
        )   << "Field size does not correspond to size of IB points "
            << "triangulated surface for object " << objectNames(objectID) << nl
            << "Field size = " << ibValues.size()
            << " IB points size = " << ibCellsList()[objectID].size()
            << abort(FatalError);
    }

    // for mpi, sending hf and ibp to master

    DynamicList<Type> masterIbValues;

    if(Pstream::parRun())
    {
        List< DynamicList<Type>> proIbValues(Pstream::nProcs());
        proIbValues[Pstream::myProcNo()]=ibValues;
        proIbValues[Pstream::myProcNo()].shrink();
        Pstream::gatherList(proIbValues);
        Pstream::scatterList(proIbValues);
       
        forAll(proIbValues,procI)
        {
            
            masterIbValues.append(proIbValues[procI]);
            
        }

        
    }
    else
    {
        masterIbValues.append(ibValues);
    }
       masterIbValues.shrink();  
    const labelListList& ctfAddr = ibCellsToTriAddrList()[objectID];
    const scalarListList& ctfWeights =ibCellsToTriWeightsList()[objectID];



    Field<Type> ibPsi(ctfAddr.size(), pTraits<Type>::zero);

    // Do interpolation
    forAll (ctfAddr, triI)
    {
        const labelList& curAddr = ctfAddr[triI];
        const scalarList& curWeights = ctfWeights[triI];

        forAll (curAddr, cellI)
        {
            ibPsi[triI] += curWeights[cellI]*masterIbValues[curAddr[cellI]];
           
        }
         
       
        
    }
    
  
    const double Oldtime1=time().elapsedCpuTime();
    
    
   if(Pstream::parRun())
   {
	// for mpi, sending hf and ibp to each processor

        List<DynamicList<Type>> proIbPsi(Pstream::nProcs());
       
        const labelListList& proTri=proTriFacesInMeshList()[objectID];
        proIbPsi[Pstream::myProcNo()]=ibPsi;
        
        proIbPsi[Pstream::myProcNo()].shrink();
        
        Pstream::gatherList(proIbPsi);
        Pstream::scatterList(proIbPsi);
            forAll(proTri,procI)
            {
                   
                    forAll (proTri[procI], triI)
                    {
                        label triID = proTri[procI][triI];
                    	ibPsi[triID] = proIbPsi[procI][triID];
                    }
            }

    } 
      label tmp = Pstream::myProcNo();
        reduce(tmp,sumOp<label>());

    
     const double Oldtime2=time().elapsedCpuTime();
     if(debug)
     {
        Info<<"mapFromIb Executation Time 4 = "<<Oldtime2-Oldtime1<< " s"<<endl;
     }

    return ibPsi;
}

template<class Type>
void Foam::immersedBoundaryFvMesh::solveMatrix
(
    Type& psi,
    const scalarSquareMatrix& matrix,
    const Field<Type>& source,
    const point& point,
    const scalarList& weight

)const
{
    // Loop over field components

    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        scalarField sourceCmpt(source.component(cmpt));

        scalarRectangularMatrix tmpMatrix=matrix;
        scalarRectangularMatrix invMatrix(SVDinv(tmpMatrix,0.0));

        scalarField coeffCmpt(sourceCmpt.size(),0.0);
        for (label i = 0; i < invMatrix.m(); i++)
        {
            for (label j = 0; j < invMatrix.n(); j++)
            {
                coeffCmpt[j]+=matrix[i][j]*sourceCmpt[i];
            }
        }

        scalar X = point.x();
        scalar Y = point.y();
        scalar Z = point.z();

        scalar psiCmpt = 0.0;
        label coeff = 0;
        psiCmpt += coeffCmpt[coeff++]*1.0;
        psiCmpt += coeffCmpt[coeff++]*X;
        psiCmpt += coeffCmpt[coeff++]*Y;
        psiCmpt += coeffCmpt[coeff++]*Z;
        psiCmpt += coeffCmpt[coeff++]*X*Y;
        psiCmpt += coeffCmpt[coeff++]*X*Z;
        psiCmpt += coeffCmpt[coeff++]*Y*Z;
        psiCmpt += coeffCmpt[coeff++]*X*Y*Z;
           setComponent(psi, cmpt)=psiCmpt;
    }

    
}

// global will be done totally
template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::mapFromGHOSTHitToTriFace
(
    const Field<Type>& ghostValues,
    label objectID
)const
{
    if (ghostValues.size() != ghostCellsList()[objectID].size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::Field<Type>\n"
            "immersedBoundaryFvMesh::mapFromGHOSTHitToTriFace\n"
            "(\n"
            "    const Field<Type>& ibValues,\n"
            "    label objectID\n"
            ") const"
        )   << "Field size does not correspond to size of IB points "
            << "triangulated surface for object " << objectNames(objectID) << nl
            << "Field size = " << ghostValues.size()
            << " IB points size = " << ghostCellsList()[objectID].size()
            << abort(FatalError);
    }

    // for mpi, sending hf and ibp to master
    List< DynamicList<Type>> proGhostValues(Pstream::nProcs());
    DynamicList<Type> masterGhostValues;
    masterGhostValues.append(ghostValues);

    if(Pstream::parRun())
    {
        proGhostValues[Pstream::myProcNo()]=ghostValues;

        Pstream::gatherList(proGhostValues);
        Pstream::scatterList(proGhostValues);
        
        forAll(proGhostValues,procI)
        {
            if(procI!=Pstream::myProcNo())
            {
                masterGhostValues.append(proGhostValues[procI].shrink());
            }
        }
        
    }
    masterGhostValues.shrink();
    const labelListList& ctfAddr = ghostCellsToTriAddrList()[objectID];
    const scalarListList& ctfWeights = ghostCellsToTriWeightsList()[objectID];

    Field<Type> ghostPsi(ctfAddr.size(), pTraits<Type>::zero);

    // Do interpolation
    forAll (ctfAddr, triI)
    {
        const labelList& curAddr = ctfAddr[triI];
        const scalarList& curWeights = ctfWeights[triI];

        forAll (curAddr, cellI)
        {
            ghostPsi[triI] += curWeights[cellI]*masterGhostValues[curAddr[cellI]];
        }
    }

    return ghostPsi;
}

// global will be done totally
template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::mapFromIBHitToTriEdge
(
    const Field<Type>& ibValues,
    label objectID
)const
{
    if (ibValues.size() != ibCellsList()[objectID].size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::Field<Type>\n"
            "immersedBoundaryFvMesh::mapFromIBHitToTriEdge\n"
            "(\n"
            "    const Field<Type>& ibValues,\n"
            "    label objectID\n"
            ") const"
        )   << "Field size does not correspond to size of IB points "
            << "triangulated surface for object " << objectNames(objectID) << nl
            << "Field size = " << ibValues.size()
            << " IB points size = " << ibCellsList()[objectID].size()
            << abort(FatalError);
    }

    // for mpi, sending hf and ibp to master

    DynamicList<Type> masterIbValues;


    if(Pstream::parRun())
    {
        List< DynamicList<Type>> proIbValues(Pstream::nProcs());
        proIbValues[Pstream::myProcNo()]=ibValues;
        proIbValues[Pstream::myProcNo()].shrink();
        Pstream::gatherList(proIbValues);
        Pstream::scatterList(proIbValues);
       
        forAll(proIbValues,procI)
        {
            
            masterIbValues.append(proIbValues[procI]);
            
        }
        
    }
    else
    {
        masterIbValues.append(ibValues);
    }
    masterIbValues.shrink();
    const labelListList& ctfAddr = ibCellsToTriEdgeAddrList()[objectID];
    const scalarListList& ctfWeights =ibCellsToTriEdgeWeightsList()[objectID];


    Field<Type> ibPsi(ctfAddr.size(), pTraits<Type>::zero);

    // Do interpolation
    forAll (ctfAddr, triI)
    {
        const labelList& curAddr = ctfAddr[triI];
        const scalarList& curWeights = ctfWeights[triI];

        forAll (curAddr, cellI)
        {
            ibPsi[triI] += curWeights[cellI]*masterIbValues[curAddr[cellI]];
        }
    }

    return ibPsi;
}

// global will be done totally
template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::mapFromGHOSTHitToTriEdge
(
    const Field<Type>& ghostValues,
    label objectID
)const
{
    if (ghostValues.size() != ghostCellsList()[objectID].size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::Field<Type>\n"
            "immersedBoundaryFvMesh::mapFromGHOSTHitToTriEdge\n"
            "(\n"
            "    const Field<Type>& ibValues,\n"
            "    label objectID\n"
            ") const"
        )   << "Field size does not correspond to size of IB points "
            << "triangulated surface for object " << objectNames(objectID) << nl
            << "Field size = " << ghostValues.size()
            << " IB points size = " << ghostCellsList()[objectID].size()
            << abort(FatalError);
    }

    // for mpi, sending hf and ibp to master
    List< DynamicList<Type>> proGhostValues(Pstream::nProcs());
    DynamicList<Type> masterGhostValues;
    masterGhostValues.append(ghostValues);

    if(Pstream::parRun())
    {
        proGhostValues[Pstream::myProcNo()]=ghostValues;
        Pstream::gatherList(proGhostValues);
        Pstream::scatterList(proGhostValues);
        
        forAll(proGhostValues,procI)
        {
            if(procI!=Pstream::myProcNo())
            {
                masterGhostValues.append(proGhostValues[procI].shrink());
            }
        }
        
    }
    masterGhostValues.shrink();
    const labelListList& ctfAddr = ghostCellsToTriEdgeAddrList()[objectID];
    const scalarListList& ctfWeights = ghostCellsToTriEdgeWeightsList()[objectID];

    Field<Type> ghostPsi(ctfAddr.size(), pTraits<Type>::zero);

    // Do interpolation
    forAll (ctfAddr, triI)
    {
        const labelList& curAddr = ctfAddr[triI];
        const scalarList& curWeights = ctfWeights[triI];

        forAll (curAddr, cellI)
        {
            ghostPsi[triI] += curWeights[cellI]*masterGhostValues[curAddr[cellI]];
        }
    }

    return ghostPsi;
}

// global will be done totally
template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::mapFromIBHitToDualEdge
(
    const Field<Type>& ibValues,
    label objectID
)const
{
    if (ibValues.size() != ibCellsList()[objectID].size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::Field<Type>\n"
            "immersedBoundaryFvMesh::mapFromIBHitToDualEdge\n"
            "(\n"
            "    const Field<Type>& ibValues,\n"
            "    label objectID\n"
            ") const"
        )   << "Field size does not correspond to size of IB points "
            << "dual surface for object " << objectNames(objectID) << nl
            << "Field size = " << ibValues.size()
            << " IB points size = " << ibCellsList()[objectID].size()
            << abort(FatalError);
    }
    
    // Map values to triFaceCenter/dualFaceVertice
  
    Field<Type> triFacePsi = mapFromIBHitToTriFace(ibValues,objectID);
    
    // find relationship between triFaceCenter/dualFaceVertice and dualmeshface/point
    const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];
    const dualPatch& dPatch = cTSM.newDualPatch();
    const PPatch& patch = dPatch.newPatch();
    const labelList& mFFTDP = dPatch.mapFromFacesToDualPoints();
    
    //dual patch point values, whole patch
    Field<Type> dualPtValues(patch.localPoints().size(), pTraits<Type>::zero);
    
    forAll(mFFTDP,I)
    {
        if(mFFTDP[I]>-1)
        {
            dualPtValues[mFFTDP[I]] = triFacePsi[I];
        }
        
        
    }
    
    // dual patch edge values
    const edgeList& edges = patch.edges();
    
    Field<Type> ibPsi(edges.size(), pTraits<Type>::zero);
    
    forAll(edges,I)
    {
        edge e = edges[I];
        ibPsi[I] = 0.5*(dualPtValues[e[0]]+dualPtValues[e[1]]);
        
    }
    
   

    return ibPsi;
}

// global will be done totally
template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::mapFromGHOSTHitToDualEdge
(
    const Field<Type>& ghostValues,
    label objectID
)const
{
    if (ghostValues.size() != ghostCellsList()[objectID].size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::Field<Type>\n"
            "immersedBoundaryFvMesh::mapFromGHOSTHitToDualEdge\n"
            "(\n"
            "    const Field<Type>& ibValues,\n"
            "    label objectID\n"
            ") const"
        )   << "Field size does not correspond to size of IB points "
            << "dual surface for object " << objectNames(objectID) << nl
            << "Field size = " << ghostValues.size()
            << " IB points size = " << ghostCellsList()[objectID].size()
            << abort(FatalError);
    }

    // Map values to triFaceCenter/dualFaceVertice
    Field<Type> triFacePsi = mapFromGHOSTHitToTriFace(ghostValues,objectID);
    
    // find relationship between triFaceCenter/dualFaceVertice and dualmeshface/point
    const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];
    const dualPatch& dPatch = cTSM.newDualPatch();
    const PPatch& patch = dPatch.newPatch();
    const labelList& mFFTDP = dPatch.mapFromFacesToDualPoints();
    
    //dual patch point values, whole patch
    Field<Type> dualPtValues(patch.localPoints().size(), pTraits<Type>::zero);
    
    forAll(mFFTDP,I)
    {
        if(mFFTDP[I]>-1)
        {
            dualPtValues[mFFTDP[I]] = triFacePsi[I];
        }
    }
    
    // dual patch edge values
    const edgeList& edges = patch.edges();
    
    Field<Type> ghostPsi(edges.size(), pTraits<Type>::zero);
    
    forAll(edges,I)
    {
        edge e = edges[I];
        ghostPsi[I] = 0.5*(dualPtValues[e[0]]+dualPtValues[e[1]]);
    }

    return ghostPsi;
}


// global will be done totally
template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::mapFromIBHitToDualFace
(
    const Field<Type>& ibValues,
    label objectID
)const
{
    if (ibValues.size() != ibCellsList()[objectID].size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::Field<Type>\n"
            "immersedBoundaryFvMesh::mapFromIBHitToDualFace\n"
            "(\n"
            "    const Field<Type>& ibValues,\n"
            "    label objectID\n"
            ") const"
        )   << "Field size does not correspond to size of IB points "
            << "dual surface for object " << objectNames(objectID) << nl
            << "Field size = " << ibValues.size()
            << " IB points size = " << ibCellsList()[objectID].size()
            << abort(FatalError);
    }
    // Map values to triFaceCenter/dualFaceVertice
    Field<Type> triFacePsi = mapFromIBHitToTriFace(ibValues,objectID);
        
    // find relationship between triFaceCenter/dualFaceVertice and dualmeshface/point
    const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];
    const dualPatch& dPatch = cTSM.newDualPatch();
    const PPatch& patch = dPatch.newPatch();
    const labelList& mFFTDP = dPatch.mapFromFacesToDualPoints();
        
    Field<Type> dualPtValues(patch.localPoints().size(), pTraits<Type>::zero);
    
    forAll(mFFTDP,I)
    {
        if(mFFTDP[I]>-1)
        {
            dualPtValues[mFFTDP[I]] = triFacePsi[I];
        }
    }
    
    const PPatchInterpolation& ppI = dPatch.dualPatchInterp();
    
    Field<Type> ibPsi = ppI.pointToFaceInterpolate(dualPtValues);
    
    return ibPsi;
}

// global will be done totally
template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::mapFromGHOSTHitToDualFace
(
    const Field<Type>& ghostValues,
    label objectID
)const
{
    if (ghostValues.size() != ghostCellsList()[objectID].size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::Field<Type>\n"
            "immersedBoundaryFvMesh::mapFromGHOSTHitToDualFace\n"
            "(\n"
            "    const Field<Type>& ibValues,\n"
            "    label objectID\n"
            ") const"
        )   << "Field size does not correspond to size of IB points "
            << "dual surface for object " << objectNames(objectID) << nl
            << "Field size = " << ghostValues.size()
            << " IB points size = " << ghostCellsList()[objectID].size()
            << abort(FatalError);
    }

    // Map values to triFaceCenter/dualFaceVertice
    Field<Type> triFacePsi = mapFromGHOSTHitToTriFace(ghostValues,objectID);
    
    // find relationship between triFaceCenter/dualFaceVertice and dualmeshface/point
    const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];
    const dualPatch& dPatch = cTSM.newDualPatch();
    const PPatch& patch = dPatch.newPatch();
    const labelList& mFFTDP = dPatch.mapFromFacesToDualPoints();
    
    Field<Type> dualPtValues(patch.localPoints().size(), pTraits<Type>::zero);
    
    forAll(mFFTDP,I)
    {
        if(mFFTDP[I]>-1)
        {
            dualPtValues[mFFTDP[I]] = triFacePsi[I];
        }
    }
    
    const PPatchInterpolation& ppI = dPatch.dualPatchInterp();
    
    Field<Type> ghostPsi = ppI.pointToFaceInterpolate(dualPtValues);

    return ghostPsi;
}

// not used any more
template<class T>
void Foam::immersedBoundaryFvMesh::sendToMaster
(
    List<T>& Values
) const
{
    if (Pstream::parRun())
    {
        if (Pstream::myProcNo() != Pstream::masterNo())
        {
            
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
                
                 // create the input stream from processor procI
                 IPstream iStream(Pstream::commsTypes::blocking, procI);
                 iStream >> Values[procI];
                

            }

        }
    }
}
// ************************************************************************* //

