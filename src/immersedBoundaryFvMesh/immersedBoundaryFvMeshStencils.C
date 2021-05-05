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
    This is used to obtain:
    ibHitPointsListPtr_, ghostHitPointsListPtr_
    ibHitFacesListPtr_, ghostHitFacesListPtr_
    samplingPointsListPtr_, imagePointsListPtr_
SourceFiles
    immersedBoundaryFvMeshStencils.C
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvMesh.H"
#include "MeshWave.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::immersedBoundaryFvMesh::makeStencilsInfo() const
{
    if (debug==2)
    {
        InfoIn("void immersedBoundaryFvMesh::makeStencils() const")
            << "make stencils for IB and ghost cells "
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibHitPointsListPtr_ || ghostHitPointsListPtr_ ||
    ibHitFacesListPtr_ || ghostHitFacesListPtr_ ||
    samplingPointsListPtr_ || farSamplingPointsListPtr_ ||imagePointsListPtr_ )
    {
        FatalErrorIn("immersedBoundaryFvMesh::readDict() const")
            << "make stencils for IB and ghost cells "
            << "ibHitPointsListPtr_ || ghostHitPointsListPtr_ ||"
            << "ibHitFacesListPtr_ || ghostHitFacesListPtr_ ||"
            <<"samplingPointsListPtr_ || imagePointsListPtr_ ||"
            << " already exists"
            << abort(FatalError);
    }

    ibHitPointsListPtr_ = new PtrList<pointField> (objectsList().size());
    ibHitFacesListPtr_ = new PtrList<labelList> (objectsList().size());
    ibTriNormalsFlipListPtr_ = new PtrList<bool> (objectsList().size());
    samplingPointsListPtr_ = new PtrList<pointField> (objectsList().size());
    farSamplingPointsListPtr_ = new PtrList<pointField> (objectsList().size());

    ibFaceHitPointsListPtr_ = new PtrList<pointField> (objectsList().size());
    ibFaceHitFacesListPtr_ = new PtrList<labelList> (objectsList().size());

    ibFaceSamplingPointsListPtr_ = new PtrList<pointField> (objectsList().size());



    ghostHitPointsListPtr_ = new PtrList<pointField> (objectsList().size());
    ghostHitFacesListPtr_ = new PtrList<labelList> (objectsList().size());
    imagePointsListPtr_ = new PtrList<pointField> (objectsList().size());

    samplingStencilsListPtr_= new PtrList<immersedBoundaryStencils>
            (objectDictList().size());
    farSamplingStencilsListPtr_= new PtrList<immersedBoundaryStencils>
            (objectDictList().size());
    ibCellStencilsListPtr_= new PtrList<immersedBoundaryStencils>
            (objectDictList().size());

    ibFaceSamplingStencilsListPtr_= new PtrList<immersedBoundaryStencils>
            (objectDictList().size());

    imageStencilsListPtr_= new PtrList<immersedBoundaryStencils>
            (objectDictList().size());

   
}

Foam::labelList Foam::immersedBoundaryFvMesh::findCellID(const pointField& points) const
{
    // make sure parallel works
    labelList cellIDList(points.size(),-1);

    // if the point is out of domain, then find the nearest cell
    if (Pstream::parRun())//Start of mpi run
    {
        // It is a compromise. It can be more accurate but more costy.
        forAll(cellIDList,I)
        {
            cellIDList[I]=queryMeshPtr_->findCell(points[I], -1, true);
            if(cellIDList[I]<0)cellIDList[I]=queryMeshPtr_->findNearestCell(points[I], -1, true);
        }
    }
    else
    {
        forAll(cellIDList,I)
        {
            cellIDList[I]=queryMeshPtr_->findCell(points[I], -1, true);
            if(cellIDList[I]<0)cellIDList[I]=queryMeshPtr_->findNearestCell(points[I], -1, true);
        }
    }

    return cellIDList;
}

Foam::scalar Foam::immersedBoundaryFvMesh::cellSize(label cellID) const
{
    scalar delta;

    if (this->nGeometricD() == 3)
    {
        delta = Foam::pow(this->V().field()[cellID], 1.0/3.0);
    }
    else
    {
        scalar thickness = 0.0;
        const Vector<label>& directions = this->geometricD();
        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = this->bounds().span()[dir];
                break;
            }
        }

        delta = Foam::sqrt(this->V().field()[cellID]/thickness);


    }
    return delta;
}




Foam::scalar Foam::immersedBoundaryFvMesh::cellLength(label cellID) const
{

    const edgeList& edges = this->edges();
    const pointField& pp = this->points();
    const labelList& cEdges = this->cellEdges()[cellID];

    // Make a list with all the edge lenghts
    scalarField eLengths(cEdges.size(), 0.0);

    forAll (cEdges, edgei)
    {
        label edgeID = cEdges[edgei];
        eLengths[edgei] = edges[edgeID].mag(pp);
    }
    
    // Find the minimum per cell
    scalar  minLength = GREAT;

     forAll (cEdges, edgei)
     {
         minLength = Foam::min(minLength, eLengths[edgei]);
     }



    return minLength;
}









void Foam::immersedBoundaryFvMesh::makeIbHitAndSamplingPoints(const label& objectID) const
{
    radiusFactor_ = ibProperties().lookupOrDefault<scalar>("radiusFactor",3.0);
    radiusFactor_ = objectDictList()[objectID].lookupOrDefault<scalar>("radiusFactor",radiusFactor_);


    distFactor_ = ibProperties().lookupOrDefault<scalar>("distFactor",3);
    distFactor_ = objectDictList()[objectID].lookupOrDefault<scalar>("distFactor",distFactor_);
    

    


    const labelList& ibc = ibCellsList()[objectID];

    const triSurfaceSearch& tss = triSurfaceSearchList()[objectID];

    pointField hitPts(ibc.size());
    labelList hitFaces(ibc.size());
    pointField samplingPts(ibc.size());
    pointField farSamplingPts(ibc.size());
    vectorField surfNormals(ibc.size());
    bool ibTriNormalsFlip(false);
    forAll(ibc, I)
    {
        label cellID = ibc[I];

        vector span
        (
            2*radiusFactor_*cellSize(cellID),
            2*radiusFactor_*cellSize(cellID),
            2*radiusFactor_*cellSize(cellID)
        );

        pointIndexHit pih = tss.nearest(this->cellCentres()[cellID], span);

        if (pih.hit())
        {
            hitPts[I] = pih.hitPoint();
            surfNormals[I] =
                triSurfaceTools::surfaceNormal
                (
                    objectsList()[objectID],
                    pih.index(),
                    pih.hitPoint()
                );
            hitFaces[I] = pih.index();
        }
        else
        {
            FatalErrorIn
                (
                    "Foam::bool Foam::immersedBoundaryFvMesh::makeIbHitAndSamplingPoints"
                )   << "cell " << cellID<< " at location " << this->cellCentres()[cellID]
                    <<" does not has hit point on object "<< objectNames(objectID)
                    << abort(FatalError);
        }

          vector normal_= this->cellCentres()[cellID] - hitPts[I];
          normal_=normal_/mag(normal_);

           
         samplingPts[I] = hitPts[I]
                        + distFactor_* cellLength(cellID)*normal_;
       
    }

    if(samplingPts.size()>0)
    {
        scalar tmpR= surfNormals[0]&(samplingPts[0]-hitPts[0]);
        if(tmpR<0)
        {
            ibTriNormalsFlip=true;
        }
    }
    ibHitPointsListPtr_->set
    (
        objectID,
        new pointField(hitPts)
    );

    ibHitFacesListPtr_->set
    (
        objectID,
        new labelList(hitFaces)
    );

    ibTriNormalsFlipListPtr_->set
    (
        objectID,
        new bool(ibTriNormalsFlip)
    );

    samplingPointsListPtr_->set
    (
        objectID,
        new pointField(samplingPts)
    );
    farSamplingPointsListPtr_->set
    (
        objectID,
        new pointField(farSamplingPts)
    );
}

void Foam::immersedBoundaryFvMesh::makeIbFacesHitAndSamplingPoints(const label& objectID) const
{
    radiusFactor_ = ibProperties().lookupOrDefault<scalar>("radiusFactor",3.0);
    radiusFactor_ = objectDictList()[objectID].lookupOrDefault<scalar>("radiusFactor",radiusFactor_);


    distFactor_ = ibProperties().lookupOrDefault<scalar>("distFactor",1.5);
    distFactor_ = objectDictList()[objectID].lookupOrDefault<scalar>("distFactor",distFactor_);




    const labelList& ibf = ibFacesList()[objectID];
    const labelList& ibc = ibCellsList()[objectID];
    const labelList& ibfc = ibFaceCellsList()[objectID];

    const triSurfaceSearch& tss = triSurfaceSearchList()[objectID];

    pointField hitPts(ibf.size());
    labelList hitFaces(ibf.size());
    pointField samplingPts(ibf.size());
    vectorField surfNormals(ibf.size());
//    bool ibTriNormalsFlip(false);
    forAll(ibf, I)
    {
        label faceID = ibf[I];
        label cellID = ibc[max(0,ibfc[I])];


        scalar cellSize =  this->cellSize(cellID);

        vector span
        (
            2*radiusFactor_*cellSize,
            2*radiusFactor_*cellSize,
            2*radiusFactor_*cellSize
        );

        pointIndexHit pih = tss.nearest(this->faceCentres()[faceID], span);

        if (pih.hit())
        {
            hitPts[I] = pih.hitPoint();
            surfNormals[I] =
                triSurfaceTools::surfaceNormal
                (
                    objectsList()[objectID],
                    pih.index(),
                    pih.hitPoint()
                );
            hitFaces[I] = pih.index();
        }
        else
        {
            FatalErrorIn
                (
                    "Foam::bool Foam::immersedBoundaryFvMesh::makeIbFacesHitAndSamplingPoints"
                )   << "face " << faceID<< " at location " << this->faceCentres()[faceID]
                    <<" does not has hit point on object "<< objectNames(objectID)
                    << abort(FatalError);
        }
        vector normal_= this->faceCentres()[faceID] - hitPts[I];
        normal_.z()+=SMALL;
        normal_=normal_/mag(normal_);
    
        
        samplingPts[I] = hitPts[I]
                    + distFactor_* cellLength(cellID)*normal_;
                    

    }

    ibFaceHitPointsListPtr_->set
    (
        objectID,
        new pointField(hitPts)
    );

    ibFaceHitFacesListPtr_->set
    (
        objectID,
        new labelList(hitFaces)
    );

    ibFaceSamplingPointsListPtr_->set
    (
        objectID,
        new pointField(samplingPts)
    );
}

void Foam::immersedBoundaryFvMesh::makeGhostHitAndImagePoints(const label& objectID) const
{
    scalar radiusFactor = ibProperties().lookupOrDefault<scalar>("radiusFactor",3.0);
    radiusFactor = objectDictList()[objectID].lookupOrDefault<scalar>("radiusFactor",radiusFactor);

    const labelList& gc = ghostCellsList()[objectID];

    const triSurfaceSearch& tss = triSurfaceSearchList()[objectID];

    pointField hitPts(gc.size());
    labelList hitFaces(gc.size());
    pointField imagePts(gc.size());

    forAll(gc, I)
    {
        label cellID = gc[I];

        vector span
        (
            2*radiusFactor*cellSize(cellID),
            2*radiusFactor*cellSize(cellID),
            2*radiusFactor*cellSize(cellID)
        );

        pointIndexHit pih = tss.nearest(this->cellCentres()[cellID], span);

        if (pih.hit())
        {
            hitPts[I] = pih.hitPoint();

            hitFaces[I] = pih.index();
        }
        vector normal_= hitPts[I]-this->cellCentres()[cellID];
        normal_=normal_/mag(normal_);
        imagePts[I] = hitPts[I] + distFactor_*cellSize(cellID)*normal_;
    }

    ghostHitPointsListPtr_->set
    (
        objectID,
        new pointField(hitPts)
    );

    ghostHitFacesListPtr_->set
    (
        objectID,
        new labelList(hitFaces)
    );

    imagePointsListPtr_->set
    (
        objectID,
        new pointField(imagePts)
    );
}

void Foam::immersedBoundaryFvMesh::makeSamplingStencils(const label& objectID) const
{
    maxCellCellRows_ =  ibProperties().lookupOrDefault<scalar>("maxCellCellRows",4.0);
    maxCellCellRows_ = objectDictList()[objectID].lookupOrDefault<scalar>("maxCellCellRows",maxCellCellRows_);

    radiusFactor_ = ibProperties().lookupOrDefault<scalar>("radiusFactor",3.0);
    radiusFactor_ = objectDictList()[objectID].lookupOrDefault<scalar>("radiusFactor",radiusFactor_);

    dictionary dict;
    dict.add("maxCellCellRows",maxCellCellRows_);
    dict.add("radiusFactor",radiusFactor_);
    dict.add("gammaInclude_c",-0.5);
    dict.add("gammaWeight_c",0.0);

    samplingStencilsListPtr_->set
    (
        objectID,
        new immersedBoundaryStencils
        (
            "sampling points",
            samplingPointsList()[objectID],
            findCellID(samplingPointsList()[objectID]),
            ibGammaList()[objectID],
            queryMesh(),
            dict
        )
    );
}





void Foam::immersedBoundaryFvMesh::makeFarSamplingStencils(const label& objectID) const
{
    maxCellCellRows_ =  ibProperties().lookupOrDefault<scalar>("maxCellCellRows",4.0);
    maxCellCellRows_ = objectDictList()[objectID].lookupOrDefault<scalar>("maxCellCellRows",maxCellCellRows_);

    radiusFactor_ = ibProperties().lookupOrDefault<scalar>("radiusFactor",3.0);
    radiusFactor_ = objectDictList()[objectID].lookupOrDefault<scalar>("radiusFactor",radiusFactor_);

    dictionary dict;
    dict.add("maxCellCellRows",maxCellCellRows_);
    dict.add("radiusFactor",radiusFactor_);
    dict.add("gammaInclude_c",-0.5);
    dict.add("gammaWeight_c",0.0);

    farSamplingStencilsListPtr_->set
    (
        objectID,
        new immersedBoundaryStencils
        (
            "farSampling points",
            farSamplingPointsList()[objectID],
            findCellID(samplingPointsList()[objectID]),
            ibGammaList()[objectID],
            queryMesh(),
            dict
        )
    );
}




void Foam::immersedBoundaryFvMesh::makeIbCellStencils(const label& objectID) const
{
    maxCellCellRows_ =  ibProperties().lookupOrDefault<scalar>("maxCellCellRows",4.0);
    maxCellCellRows_ = objectDictList()[objectID].lookupOrDefault<scalar>("maxCellCellRows",maxCellCellRows_);

    radiusFactor_ = ibProperties().lookupOrDefault<scalar>("radiusFactor",3.0);
    radiusFactor_ = objectDictList()[objectID].lookupOrDefault<scalar>("radiusFactor",radiusFactor_);

    dictionary dict;
    dict.add("maxCellCellRows",maxCellCellRows_);
    dict.add("radiusFactor",radiusFactor_);
    dict.add("gammaInclude_c",-0.5);
    dict.add("gammaWeight_c",0.0);

    const vectorField& C = this->cellCentres();
    pointField ibc_centers(C,ibCellsList()[objectID]);

    ibCellStencilsListPtr_->set
    (
        objectID,
        new immersedBoundaryStencils
        (
            "IB cells",
            ibc_centers,
            ibCellsList()[objectID],
            ibGammaList()[objectID],
            queryMesh(),
            dict
        )
    );
}

void Foam::immersedBoundaryFvMesh::makeIbFaceSamplingStencils(const label& objectID) const
{
    maxCellCellRows_ =  ibProperties().lookupOrDefault<scalar>("maxCellCellRows",4.0);
    maxCellCellRows_ = objectDictList()[objectID].lookupOrDefault<scalar>("maxCellCellRows",maxCellCellRows_);

    radiusFactor_ = ibProperties().lookupOrDefault<scalar>("radiusFactor",3.0);
    radiusFactor_ = objectDictList()[objectID].lookupOrDefault<scalar>("radiusFactor",radiusFactor_);

    dictionary dict;
    dict.add("maxCellCellRows",maxCellCellRows_);
    dict.add("radiusFactor",radiusFactor_);
    dict.add("gammaInclude_c",-0.5);
    dict.add("gammaWeight_c",0.0);

    ibFaceSamplingStencilsListPtr_->set
    (
        objectID,
        new immersedBoundaryStencils
        (
            "IB face sampling points",
            ibFaceSamplingPointsList()[objectID],
            findCellID(ibFaceSamplingPointsList()[objectID]),
            ibGammaList()[objectID],
            queryMesh(),
            dict
        )
    );
}

void Foam::immersedBoundaryFvMesh::makeImageStencils(const label& objectID) const
{
    maxCellCellRows_ =  ibProperties().lookupOrDefault<scalar>("maxCellCellRows",4.0);
    maxCellCellRows_ = objectDictList()[objectID].lookupOrDefault<scalar>("maxCellCellRows",maxCellCellRows_);

    radiusFactor_ = ibProperties().lookupOrDefault<scalar>("radiusFactor",3.0);
    radiusFactor_ = objectDictList()[objectID].lookupOrDefault<scalar>("radiusFactor",radiusFactor_);

    dictionary dict;
    dict.add("maxCellCellRows",maxCellCellRows_);
    dict.add("radiusFactor",radiusFactor_);
    dict.add("gammaInclude_c",-0.5);
    dict.add("gammaWeight_c",0.0);

    imageStencilsListPtr_->set
    (
        objectID,
        new immersedBoundaryStencils
        (
            "image points",
            imagePointsList()[objectID],
            findCellID(imagePointsList()[objectID]),
            ibGammaList()[objectID],
            queryMesh(),
            dict
        )
    );
}

const Foam::immersedBoundaryStencils& Foam::immersedBoundaryFvMesh::samplingStencils
(
    const label& objectID
) const
{
    if(!samplingStencilsListPtr_->set(objectID))
    {
        makeSamplingStencils(objectID);
    }

    PtrList<immersedBoundaryStencils>& samplingStencilsList = *samplingStencilsListPtr_;

    return samplingStencilsList[objectID];
}

const Foam::immersedBoundaryStencils& Foam::immersedBoundaryFvMesh::farSamplingStencils
(
    const label& objectID
) const
{
    if(!farSamplingStencilsListPtr_->set(objectID))
    {
        makeFarSamplingStencils(objectID);
    }

    PtrList<immersedBoundaryStencils>& farSamplingStencilsList = *farSamplingStencilsListPtr_;

    return farSamplingStencilsList[objectID];
}


const Foam::immersedBoundaryStencils& Foam::immersedBoundaryFvMesh::ibCellStencils
(
    const label& objectID
) const
{
    if(!ibCellStencilsListPtr_->set(objectID))
    {
        makeIbCellStencils(objectID);
    }

    PtrList<immersedBoundaryStencils>& ibCellStencilsList = *ibCellStencilsListPtr_;

    return ibCellStencilsList[objectID];
}

const Foam::immersedBoundaryStencils& Foam::immersedBoundaryFvMesh::ibFaceSamplingStencils
(
    const label& objectID
) const
{
    if(!ibFaceSamplingStencilsListPtr_->set(objectID))
    {
        makeIbFaceSamplingStencils(objectID);
    }

    PtrList<immersedBoundaryStencils>& ibFaceSamplingStencilsList = *ibFaceSamplingStencilsListPtr_;

    return ibFaceSamplingStencilsList[objectID];
}

const Foam::immersedBoundaryStencils& Foam::immersedBoundaryFvMesh::imageStencils
(
    const label& objectID
) const
{
    if(!imageStencilsListPtr_->set(objectID))
    {
        makeImageStencils(objectID);
    }

    PtrList<immersedBoundaryStencils>& imageStencilsList = *imageStencilsListPtr_;

    return imageStencilsList[objectID];
}


// ************************************************************************* //

