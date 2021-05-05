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
    Mark IB cells, ghost cells, live cells, dead cells, IB faces and so forth.

    This is used to obtain
    gammaCellTypeListPtr_, gammaCellTypePtr_, ibCellsListPtr_, ghostCellsListPtr_
    gammaPtr_

SourceFiles
    immersedBoundaryFvMeshMarkCells.C
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvMesh.H"
#include "MeshWave.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::immersedBoundaryFvMesh::markCells()const
{
    if (debug==2)
    {
        InfoIn("void immersedBoundaryFvMesh::markCells() const")
            << "make IB and ghost cells "
            << endl;
    }

    // initialize cell infor list as outside
    List<cellInfo> cellInfoList(this->nCells());
    cellInfoList = cellClassification::OUTSIDE;

    if (gammaCellTypeListPtr_||ibCellsListPtr_||ghostCellsListPtr_||gammaCellTypePtr_||gammaPtr_||oldIbCellsListPtr_||oldIbDeadCellsListPtr_||oldIbLiveCellsListPtr_)
    {
        FatalErrorIn("immersedBoundaryFvMesh::markCells() const")
            << "make information for cut cells, fluid cells, solid cells, IB cells, and ghost cells "
            << "gammaCellTypeListPtr_||ibCellsListPtr_||ghostCellsListPtr_||gammaCellTypePtr_||gammaPtr_||oldIbCellsListPtr_||oldIbDeadCellsListPtr_||oldIbLiveCellsListPtr_"
            << " already exists"
            << abort(FatalError);
    }

    // initialize pointer lists related to IB/Ghost/live/dead cells/faces
    // some pointer lists may not be used any more, such as oldibCellsListPtr_
    gammaCellTypeListPtr_ = new PtrList<volScalarField>(objectsList().size());
    ibCellsListPtr_ = new PtrList<labelList>(objectsList().size());
    oldIbCellsListPtr_ = new PtrList<labelList>(objectsList().size());
    ibDeadCellsListPtr_ = new PtrList<labelList>(objectsList().size());
    oldIbDeadCellsListPtr_ = new PtrList<labelList>(objectsList().size());
    ibLiveCellsListPtr_ = new PtrList<labelList>(objectsList().size());
    oldIbLiveCellsListPtr_ = new PtrList<labelList>(objectsList().size());
    ibGammaListPtr_ = new PtrList<volScalarField>(objectsList().size());
    sGammaListPtr_ = new PtrList<surfaceScalarField>(objectsList().size());
    ibFacesListPtr_ = new PtrList<labelList>(objectsList().size());
    ibFaceCellsListPtr_ = new PtrList<labelList>(objectsList().size());
    ibFaceFlipsListPtr_ = new PtrList<List<bool>>(objectsList().size());
    ghostCellsListPtr_ = new PtrList<labelList>(objectsList().size());
    
    
    // This switch is for mutilple objects moving in the fluid (e.g. sediment bed and a movable object)
    bool multipleObject_ =  ibProperties().lookupOrDefault<bool>("multipleObject",false);
    
    
    
    if (multipleObject_)
    { 
        // mark cells for each object
        forAll(objectsList(), objectID)
        {
            // check if all the PointsInFluid are inside mesh, make cut cells
            const cellClassification& cellType = this->cellType(objectID);

            // make gammaCellTypeList
            // fluid cells (1), cut cells (0), and solid cells(-1)
            makeGammaCellType(cellType,objectID);

            // make ghost cells and IB cells, as well as fluid cells and solid cells
            makeGhostAndIbCells(objectID);
        }
       
            deleteExtraIbCellsAndGhostCells();
        
       forAll(objectsList(), objectID)
       {
            // check if all the PointsInFluid are inside mesh, make cut cells
           const cellClassification& cellType = this->cellType(objectID);
            
            // make ibGamma, ibFace information
            makeIbInfo(objectID);

            // make face gamma information
            makeSGamma(objectID);

            // combine different objects
            forAll(cellInfoList,cellI)
            {
                if(cellType[cellI] != cellClassification::OUTSIDE)
                {
                    cellInfoList[cellI] = cellType[cellI];
                }
            }
        }
    }
    else
    {
         // mark cells for each object
        forAll(objectsList(), objectID)
        {
            // check if all the PointsInFluid are inside mesh, make cut cells
            const double Oldtime1=time().elapsedCpuTime();

            const cellClassification& cellType = this->cellType(objectID);
            
            const double Oldtime2=time().elapsedCpuTime();
            
            // make gammaCellTypeList
            // fluid cells (1), cut cells (0), and solid cells(-1)
            makeGammaCellType(cellType,objectID);

            const double Oldtime3=time().elapsedCpuTime();

            // make ghost cells and IB cells, as well as fluid cells and solid cells
            makeGhostAndIbCells(objectID);
            
            const double Oldtime4=time().elapsedCpuTime();

            // make ibGamma, ibFace information
            makeIbInfo(objectID);

            const double Oldtime5=time().elapsedCpuTime();


            // make face gamma information
            makeSGamma(objectID);
             
            const double Oldtime6=time().elapsedCpuTime();
                        
              
            // combine different objects
            forAll(cellInfoList,cellI)
            {
                if(cellType[cellI] != cellClassification::OUTSIDE)
                {
                    cellInfoList[cellI] = cellType[cellI];
                }
            }
            
            
            if(debug)
            {
                Info<<"cellType Executation Time = "<<Oldtime2-Oldtime1<< " s"<<endl;
                Info<<"makeGhostAndIbCells Executation Time = "<<Oldtime4-Oldtime3<< " s"<<endl;
                Info<<"makeIbInfo Executation Time = "<<Oldtime5-Oldtime4<< " s"<<endl;
                Info<<"makeSGamma Executation Time = "<<Oldtime6-Oldtime5<< " s"<<endl;

            }
        }
    }
  
    // make global gammaCellType, global live/ib/ghost/dead cells
    makeGamma(cellInfoList);
     

}

const Foam::cellClassification Foam::immersedBoundaryFvMesh::cellType
(
    const label& objectID
)const
{
    // read global pointsInFluid
    pointField globalPointsInFluid = ibProperties().lookup("pointsInFluid");

    pointField pointsInFluid(0);

    // read local pointsInFluid for each object
    pointsInFluid.append(objectDictList()[objectID].lookupOrDefault<pointField>
        ("pointsInFluid",globalPointsInFluid));

    const triSurfaceSearch& querySurf = triSurfaceSearch(objectsList()[objectID]);

    // check if pointsInFluid is inside any cell

    forAll(pointsInFluid, I)
    {
        const point& pointInFluid = pointsInFluid[I];

        label cellIndex = queryMeshPtr_->findCell(pointInFluid, -1, true);
        if(debug)
        {
            if (returnReduce(cellIndex,maxOp<label>()) < SMALL)
            {
                FatalErrorIn("immersedBoundaryFvMesh::cellType")
                    << "pointInFluid " << pointInFluid
                    << " is not inside any cell"
                    << exit(FatalError);
            }
        }
    }



    cellClassification cellType
    (
        *this,
        *queryMeshPtr_,
        querySurf,
        pointsInFluid
    );
    // cut cells has to have at least a solid or fluid neighbour cell, otherwise ifInFluid does not work
    return cellType;
}

void Foam::immersedBoundaryFvMesh::makeGammaCellType
(
    const cellClassification& cellType,
    const label& objectID
) const
{
    if (debug==2)
    {
        InfoIn("void immersedBoundaryFvMesh::makeGammaCellType() const")
            << "make information for cut cells, fluid cells, and solid cells "
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (gammaCellTypeListPtr_->set(objectID))
    {
        FatalErrorIn("immersedBoundaryFvMesh::makeGammaCellType() const")
            << "make information for cut cells, fluid cells, and solid cells "
            << "gammaCellTypeListPtr_[objectID]"
            << " already set"
            << abort(FatalError);
    }

    gammaCellTypeListPtr_->set
    (
        objectID,
        new volScalarField
            (
                IOobject
                (
                    "gammaCellType_"+ibProperties().subDict("objects").toc()[objectID],
                    this->time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                *this,
                dimensionedScalar("1", dimless, 1)
            )
    );

    PtrList<volScalarField>& gammaCellTypeList = *gammaCellTypeListPtr_;
    volScalarField& gammaCellType = gammaCellTypeList[objectID];
    scalarField& gammaCellTypeI = gammaCellType.primitiveFieldRef();
    forAll(gammaCellTypeI,cellI)
    {
        if(cellType[cellI]==cellClassification::OUTSIDE) gammaCellTypeI[cellI] = 1;
        else if(cellType[cellI]==cellClassification::CUT) gammaCellTypeI[cellI] = 0;
        else gammaCellTypeI[cellI] = -1;
    }

    // Evaluate the coupled patchField, to replace evaluateCoupled()
    // evaluateCoupled() only exist in OF-ext
    evaluateCoupled(gammaCellType);
    evaluateUnCoupled(gammaCellType);
}


// determine if a cut cell is in fluid region based on set of fluidpoints and solid points
bool Foam::immersedBoundaryFvMesh::ifInFluid
(
    const label& cellID,
    const triSurfaceSearch& querySurf,
    const label& objectID
) const
{
    bool ifInFluid = true;

    vector C = this->cellCentres()[cellID];

    scalar Delta = cellSize(cellID);

    vector span
        (
            30*Delta,
            30*Delta,
            30*Delta
        );

    pointIndexHit pih = querySurf.nearest(C, span);

    if (!pih.hit())
    {

        Delta = this->bounds().mag();
        vector span_New(Delta,Delta,Delta);
        pih = querySurf.nearest(C, span);

    }

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
            ifInFluid = true; // cellID inside
        }
        else
        {
            ifInFluid = false; // cellID outside
        }
    }
    else
    {
        FatalErrorIn
            (
                "Foam::bool Foam::immersedBoundaryFvMesh::ifInFluid"
            )   << "Can't find nearest triSurface point for point "
                << C<< ", "
                << "span = " << span
                << "\nYou could try to increase the search span. "
                << abort(FatalError);
    }

    return ifInFluid;
}


// make sure yIB is large enough.
// this function is controlled by two parameters:
// minNearWallGridSize and cellSizeRatio
// if minNearWallGridSize is not defined or negative,
// only cellSizeRatio works.
// otherwise, use minNearWallGridSize as limit
void Foam::immersedBoundaryFvMesh::correctIbCell
(
    const label& objectID,
    labelHashSet& ibCellsSet,
    labelHashSet& ghostCellsSet, // not changed
    volScalarField& gammaExt
)const
{
    bool correctIbCell = objectDictList()[objectID].lookupOrDefault<bool>
                            ("correctIbCell",false);

    if(!correctIbCell) return;

    labelHashSet oldIbCellsSet(ibCellsSet);

    Info<<"Start correcting IB cells based on IB distance"<<endl;

    const triSurfaceSearch& querySurf = triSurfaceSearchList()[objectID];

    const labelListList& cellCells = this->cellCells();

    scalar minNearWallGridSize = objectDictList()[objectID].lookupOrDefault<scalar>
        ("minNearWallGridSize",-1);// recommend to use roughness height

    scalar cellSizeRatio = objectDictList()[objectID].lookupOrDefault<scalar>
        ("cellSizeRatio",0.3);

    Info<<"  minNearWallGridSize ";
    if(minNearWallGridSize>0)
    {
        Info<<" is "<<minNearWallGridSize<<" m"<<endl;
    }
    else
    {
        Info<<" depends on cellSize, ratio is "<<cellSizeRatio<<endl;
    }

    label counter=0;
    label TaddCellsN=0;
    label TdeleteCellsN=0;
    label addCellsN=0;
    label deleteCellsN=0;

    labelHashSet proved;// already visited and proved okay with wall distance

    // loop until no adding or deleting cells
    do
    {
        const labelList ibCells=ibCellsSet.toc();
        addCellsN=0;
        deleteCellsN=0;

        forAll(ibCells,I)
        {
            label cellID=ibCells[I];
            if(proved.found(cellID))
            {
                continue;
            }
            vector C = this->cellCentres()[cellID];
            scalar Delta = cellSize(cellID);
            vector span
                (
                    30*Delta,
                    30*Delta,
                    30*Delta
                );
            pointIndexHit pih = querySurf.nearest(C, span);

            if (!pih.hit())
            {
                Delta = this->bounds().mag();
                vector span_New(Delta,Delta,Delta);
                pih = querySurf.nearest(C, span_New);
            }

            point nearestPoint = pih.hitPoint();
            scalar cDelta = cellSizeRatio*Delta;
            if(minNearWallGridSize>0) cDelta = minNearWallGridSize;
            if(mag(nearestPoint-C) < cDelta)// trigger the condition to move ib cells
            {
                labelList neiCells = cellCells[cellID];
                ibCellsSet.erase(cellID);
                gammaExt[cellID] = 0; // IB cells to delete
                deleteCellsN++;

                // find its neigh live cells to be added to IB cells
                forAll(neiCells,I)
                {
                    label celli = neiCells[I];
                    if(gammaExt[celli]>0)
                    {
                        if(!ibCellsSet.found(celli) and !oldIbCellsSet.found(celli))
                        {
                            addCellsN++;
                            ibCellsSet.insert(celli);
                            gammaExt[celli] = 0.5;
                        }
                    }
                }
            }
            else
            {
                if(!proved.found(cellID)) proved.insert(cellID);
            }
        }
        evaluateCoupled(gammaExt);
        evaluateUnCoupled(gammaExt);

        forAll (gammaExt.boundaryField(), patchI)
        {
            if (gammaExt.boundaryField()[patchI].coupled())
            {
                scalarField gammaExtOwn =
                    gammaExt.boundaryField()[patchI].patchInternalField();

                scalarField gammaExtNei =
                    gammaExt.boundaryField()[patchI].patchNeighbourField();

                const unallocLabelList& fCells =
                    this->boundary()[patchI].faceCells();

                forAll (gammaExtNei, faceI)
                {
                    if(gammaExtNei[faceI] == 0 and gammaExtOwn[faceI] > 0.5)
                    {
                        label celli = fCells[faceI];
                        if(!ibCellsSet.found(celli) and !oldIbCellsSet.found(celli))
                        {
                            addCellsN++;
                            ibCellsSet.insert(celli);
                            gammaExt[celli] = 0.5;
                        }
                    }
                }
            }
        }

        evaluateCoupled(gammaExt);
        evaluateUnCoupled(gammaExt);
        TaddCellsN +=addCellsN;
        TdeleteCellsN +=deleteCellsN;
        counter++;

    }while(returnReduce(addCellsN, maxOp<label>())>0 or returnReduce(deleteCellsN, maxOp<label>())>0);

    const labelList ibCells = ibCellsSet.toc();
    labelHashSet neighLiveCellsSet;
    forAll(ibCells,I)
    {
        const labelList& curCells = cellCells[ibCells[I]];
        forAll(curCells,II)
        {
            label cellID=curCells[II];
            if(!neighLiveCellsSet.found(cellID) and gammaExt[cellID]>0.5)
            {
                neighLiveCellsSet.insert(cellID);
            }
        }
    }


    if(debug)
    {
        Info<<"After "<<returnReduce(counter, maxOp<label>())<<" time(s) of IB cells correction: "<<endl;
        Info<<"add " << returnReduce(TaddCellsN, sumOp<label>())<<" cells"<<endl;
        Info<<"delete " << returnReduce(TdeleteCellsN, sumOp<label>())<<" cells"<<endl;
    }
}

// ghost cells still have some problem
void Foam::immersedBoundaryFvMesh::makeGhostAndIbCells(const label& objectID) const
{
    const volScalarField& gammaCellType = gammaCellTypeList()[objectID];
    const scalarField& gammaCellTypeI = gammaCellType.internalField();

    // a temporary scalar field
    volScalarField gammaExt(gammaCellType);
    gammaExt.rename("gammaExt_"+objectNames(objectID));
    scalar DeltaRatio_ = ibProperties().lookupOrDefault<scalar>
        ("DeltaRatio",10);
    pointField globalPointsInFluid = ibProperties().lookup("pointsInFluid");
    pointField pointsInFluid(0);
    pointsInFluid.append(objectDictList()[objectID].lookupOrDefault<pointField>
        ("pointsInFluid",globalPointsInFluid));

    const triSurfaceSearch& querySurf = triSurfaceSearchList()[objectID];

    // use first point
    point pIF = pointsInFluid[0];


    scalar Delta = this->bounds().mag();
    vector span_New(DeltaRatio_*Delta,DeltaRatio_*Delta,DeltaRatio_*Delta);
    pointIndexHit   pih = querySurf.nearest(pIF, span_New);


    // determine if triangle normal direction needs flip
    bool flip = false;
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

        scalar Indicator = nearestPointNormal & ( nearestPoint - pIF );

        if (Indicator>0)
        {
            flip = false;
        }
        else
        {
            flip = true;
        }
    }
    else
    {
        FatalErrorIn("immersedBoundaryFvMesh::makeGhostAndIbCells")
            << "pointInFluid " << pIF
            << " cannot find its nearest point to surface "
            << objectNames(objectID)
            << exit(FatalError);
    }

    // insert gammaExt values in cut cells
    // determine if cut cell center is in fluid or solid region
    // fluid(1), solid(-1)
    if(!flip)
    {
        forAll(gammaCellTypeI, cellI)
        {
            if(gammaCellTypeI[cellI] == 0)
            {
                if(ifInFluid(cellI,querySurf,objectID))
                {
                    gammaExt[cellI] = 1;
                }
                else
                {
                    gammaExt[cellI] = -1;
                }
            }
        }
    }
    else
    {
        forAll(gammaCellTypeI, cellI)
        {
            if(gammaCellTypeI[cellI] == 0)
            {
                if(ifInFluid(cellI,querySurf,objectID))
                {
                    gammaExt[cellI] = -1;
                }
                else
                {
                    gammaExt[cellI] = 1;
                }
            }
        }
    }

//    gammaExt.write();
    evaluateCoupled(gammaExt);
    evaluateUnCoupled(gammaExt);

    labelHashSet ibCellsSet;
    labelHashSet ghostCellsSet;

    const fvMesh& mesh_= *this;
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const volScalarField& gE = gammaExt;
    scalarField& gammaExtI = gammaExt.primitiveFieldRef();

    forAll (neighbour, faceI)
    {
        if (mag(gammaExtI[neighbour[faceI]] - gammaExtI[owner[faceI]]) > 1)
        {
            if (gammaExtI[owner[faceI]] > SMALL) // owner is fluid, neighbour is solid
            {
                if (!ibCellsSet.found(owner[faceI]))
                {
                    gammaExtI[owner[faceI]] = 0.5;
                    ibCellsSet.insert(owner[faceI]);
                }
                if (!ghostCellsSet.found(neighbour[faceI]))
                {
                    gammaExtI[neighbour[faceI]] = -0.5;
                    ghostCellsSet.insert(neighbour[faceI]);
                }
            }
            else // owner is solid, neighbour is fluid
            {
                if (!ibCellsSet.found(neighbour[faceI]))
                {
                    gammaExtI[neighbour[faceI]] = 0.5;
                    ibCellsSet.insert(neighbour[faceI]);
                }
                if (!ghostCellsSet.found(owner[faceI]))
                {
                    gammaExtI[owner[faceI]] = -0.5;
                    ghostCellsSet.insert(owner[faceI]);
                }
            }
        }
    }

    // check cells next to the processor interface
    forAll (gE.boundaryField(), patchI)
    {
        if (gE.boundaryField()[patchI].coupled())
        {
            scalarField gammaExtOwn =
                gE.boundaryField()[patchI].patchInternalField();

            scalarField gammaExtNei =
                gE.boundaryField()[patchI].patchNeighbourField();

            const unallocLabelList& fCells =
                mesh_.boundary()[patchI].faceCells();
 

            forAll (gammaExtOwn, faceI)
            {
                if( mag(gammaExtNei[faceI] - gammaExtOwn[faceI]) > 1 )
                {
                    if (gammaExtOwn[faceI] > SMALL) // owner is fluid, neighbour is solid
                    {
                        if (!ibCellsSet.found(fCells[faceI]))
                        {
                            gammaExtI[fCells[faceI]]=0.5;
                            ibCellsSet.insert(fCells[faceI]);
                        }
                    }
                    else if (2*gammaExtOwn.size() == fCells.size())
                    {
                        if
                        (
                           !ibCellsSet.found
                            (
                                fCells[gammaExtOwn.size() + faceI]
                            )
                        )
                        {
                            ibCellsSet.insert
                            (
                                fCells[gammaExtOwn.size() + faceI]
                            );
                            gammaExtI[fCells[gammaExtOwn.size() + faceI]] = 0.5;
                        }
                    }
                    else // owner is solid, neighbour is fluid
                    {
                        if (!ghostCellsSet.found(fCells[faceI]))
                        {
                            gammaExtI[fCells[faceI]] = -0.5;
                            ghostCellsSet.insert(fCells[faceI]);
                        }
                    }
                }
            }
        }
    }

    // make sure yIB is large enough.
    // this function is controlled by two parameters:
    // minNearWallGridSize and cellSizeRatio
    // if minNearWallGridSize is not defined or negative,
    // only cellSizeRatio works.
    // otherwise, use minNearWallGridSize as limit
    correctIbCell(objectID,ibCellsSet,ghostCellsSet,gammaExt);
    bool multipleObject_ =  ibProperties().lookupOrDefault<bool>("multipleObject",false);
    if(multipleObject_)
    {
        oldIbCellsListPtr_->set
            (
                objectID,
                new labelList(ibCellsSet.toc())
            );
    }
    else
    {
        ibCellsListPtr_->set
            (
                objectID,
                new labelList(ibCellsSet.toc())
            );
    }
    ghostCellsListPtr_->set
        (
            objectID,
            new labelList(ghostCellsSet.toc())
        );

    if(debug)
    {
        Info << "Number of IB cells: " << returnReduce(ibCellsSet.size(), sumOp<label>())<< endl;
        Info << "Number of ghost cells: " << returnReduce(ghostCellsSet.size(), sumOp<label>())<< endl;
    }

    labelHashSet ibDeadCellsSet; // include include ghost cells
    labelHashSet ibLiveCellsSet;
    forAll(gammaExt,cellI)
    {
        if(gammaExt[cellI]>SMALL)
        {
            if(!ibCellsSet.found(cellI) and !ibLiveCellsSet.found(cellI) )
            {
                ibLiveCellsSet.insert(cellI);
            }
        }
        else
        {
            if(!ibCellsSet.found(cellI) and !ibDeadCellsSet.found(cellI) )
            {
                ibDeadCellsSet.insert(cellI);
            }
        }
    }
    if(multipleObject_)
    {
        oldIbDeadCellsListPtr_->set
            (
                objectID,
                new labelList(ibDeadCellsSet.toc())
            );
    }
    else
    {
        ibDeadCellsListPtr_->set
            (
                objectID,
                new labelList(ibDeadCellsSet.toc())
            );
    }
    
    
    ibLiveCellsListPtr_->set
        (
            objectID,
            new labelList(ibLiveCellsSet.toc())
        );
    if(debug)
    { 
        
        Info<<"ibLiveCells number: "<<returnReduce(ibLiveCellsSet.size(), sumOp<label>())<<endl;
        Info<<"ibDeadCells number: "<<returnReduce(ibDeadCellsSet.size(), sumOp<label>())<<endl;
    }
    gammaExt.clear();
}



//Delete repete ibCells and ghostCells
void Foam::immersedBoundaryFvMesh::deleteExtraIbCellsAndGhostCells() const
{

    forAll(objectsList(), objectID)
    {
        
        const labelList& oldIbCells=oldIbCellsList()[objectID];
        
        const labelList& oldIbDeadCells=oldIbDeadCellsList()[objectID];
         
        
        labelHashSet ibCellsSet;
        labelHashSet oldIbDeadCellsSet0(oldIbDeadCells);
        forAll(objectsList(), I)
        {
            if (I != objectID)
            {
                labelHashSet oldIbDeadCellsSet(oldIbDeadCellsList()[I]); 
                forAll(oldIbCells, cellID)
                {
                    if(oldIbDeadCellsSet.found(oldIbCells[cellID]))
                    {
                        oldIbDeadCellsSet0.insert(oldIbCells[cellID]);
                    }
                    else
                    {
                        ibCellsSet.insert(oldIbCells[cellID]);
                    }
                }
            }
              
        }
        if(debug)
        {
            Info<<"ibCells number: "<<returnReduce(ibCellsSet.size(), sumOp<label>())<<endl;
            Info<<"ibDeadCells number: "<<returnReduce(oldIbDeadCellsSet0.size(), sumOp<label>())<<endl;
        }
        ibCellsListPtr_->set
        (
            objectID,
            new labelList(ibCellsSet.toc())
        );
        ibDeadCellsListPtr_->set
        (
            objectID,
            new labelList(oldIbDeadCellsSet0.toc())
        );
        
        
        
    }

    
}





// make ibGamma, ibFace information
void Foam::immersedBoundaryFvMesh::makeIbInfo(const label& objectID) const
{
    const labelList& ibCells=ibCellsList()[objectID];
    labelHashSet ibCellsSet(ibCells);
    labelHashSet ibDeadCellsSet0(ibDeadCellsList()[objectID]);
    volScalarField ibGamma
        (
            IOobject
            (
                "ibGamma_"+objectNames(objectID),
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimensionedScalar("1", dimless, 1)
        );

    scalarField& ibGammaI=ibGamma.primitiveFieldRef();

    forAll(ibGammaI,cellID)
    {
        if(ibCellsSet.found(cellID))
        {
            ibGammaI[cellID]=0;
        }
        else if(ibDeadCellsSet0.found(cellID))
        {
            ibGammaI[cellID]=-1;
        }
    }

    // Mark all live cell whose has at least one ib cell neighour ibGamma=0.5


    forAll (ibGamma.boundaryFieldRef(), patchI)
    {
        // Note: take faceCells from fvPatch (because of empty)
        if (ibGamma.boundaryFieldRef()[patchI].coupled())
        {
            scalarField gammaOwn =
                ibGamma.boundaryField()[patchI].patchInternalField();

            scalarField gammaNei =
                ibGamma.boundaryField()[patchI].patchNeighbourField();

            
         }
         
    }
    
    // Evaluate the coupled patchField, to replace evaluateCoupled()
    evaluateCoupled(ibGamma);
    evaluateUnCoupled(ibGamma);

    // make ibFaces information
    DynamicList<label> ibF;
    DynamicList<label> ibFC;
    DynamicList<bool> ibFF;

    // reverse of ibCells
    labelList ibCellIndicator(this->nCells(), -1);

    forAll (ibCells, ibCellID)
    {
        ibCellIndicator[ibCells[ibCellID]] = ibCellID;
    }

    const unallocLabelList& owner = this->owner();
    const unallocLabelList& neighbour = this->neighbour();

    forAll (neighbour, faceI)
    {
        if
        (
            ibGammaI[owner[faceI]] >SMALL
         && ibGammaI[neighbour[faceI]] ==0
        )
        {
            // Owner is live, neighbour IB.  Its IB index is in
            // ibCellIndicator
            ibF.append(faceI);
            ibFC.append(ibCellIndicator[neighbour[faceI]]); // ibFC related to ibCell ID
            ibFF.append(false);
        }
        else if
        (
             ibGammaI[owner[faceI]] ==0
          && ibGammaI[neighbour[faceI]] >SMALL
        )
        {
            // Neighbour is live, owner IB.  Its IB index is in
            // ibCellIndicator
            ibF.append(faceI);
            ibFC.append(ibCellIndicator[owner[faceI]]);
            ibFF.append(true);
        }

    }

    forAll (ibGamma.boundaryFieldRef(), patchI)
    {
        // Note: take faceCells from fvPatch (because of empty)
        const labelList& fc = this->boundary()[patchI].faceCells();
        const label start = this->boundaryMesh()[patchI].start();

        if (ibGamma.boundaryFieldRef()[patchI].coupled())
        {
            ibGamma.boundaryFieldRef()[patchI].initEvaluate(Pstream::commsTypes::blocking);
            ibGamma.boundaryFieldRef()[patchI].evaluate(Pstream::commsTypes::blocking);

            scalarField gammaOwn =
                ibGamma.boundaryField()[patchI].patchInternalField();

            scalarField gammaNei =
                ibGamma.boundaryField()[patchI].patchNeighbourField();

            forAll (gammaOwn, patchFaceI)
            {
                if
                (
                    mag(gammaNei[patchFaceI] - gammaOwn[patchFaceI]) > SMALL
                &&    gammaNei[patchFaceI]>-1 && gammaOwn[patchFaceI]>-1
                )
                {
                    if (ibCellIndicator[fc[patchFaceI]] > -1)
                    {
                        // Owner cell is IB
                        ibF.append(start + patchFaceI);
                        ibFC.append(ibCellIndicator[fc[patchFaceI]]);
                        ibFF.append(false);
                    }
                    else
                    {
                        // Neighbour cell is IB
                        ibF.append(start + patchFaceI);
                        ibFC.append(-1);
                        ibFF.append(true);
                    }
                }

            }
        }
    }

    // Pack the data
    ibF.shrink();
    ibFC.shrink();
    ibFF.shrink();

    if(debug)
    {
        Info<<"ibFaces number: "<<returnReduce(ibF.size(), sumOp<label>())<<endl;
    }
    ibFacesListPtr_->set
    (
        objectID,
        new labelList(ibF.xfer())
    );

    ibFaceCellsListPtr_->set
    (
        objectID,
        new labelList(ibFC.xfer())
    );

    ibFaceFlipsListPtr_->set
    (
        objectID,
        new boolList(ibFF.xfer())
    );

    ibGammaListPtr_->set
    (
        objectID,
        ibGamma
    );

    ibGamma.clear();
}



// make global gammaCellType, global live/ib/ghost/dead cells
void Foam::immersedBoundaryFvMesh::makeGamma
(
    const List<cellInfo>& cellInfoList
) const
{
    if (debug==2)
    {
        InfoIn("void immersedBoundaryFvMesh::makeGamma() const")
            << "make information for cut cells, fluid cells, and solid cells "
            << endl;
    }

//    // It is an error to attempt to recalculate
//    // if the pointer is already set
    if (gammaCellTypePtr_||gammaPtr_||globalLiveCellsPtr_||globalIbCellsPtr_||
        globalGhostCellsPtr_||globalDeadCellsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvMesh::makeGamma() const")
            << "make information for cut cells, fluid cells, and solid cells "
            << "gammaCellTypePtr_||gammaPtr_||globalLiveCellsPtr_||globalIbCellsPtr_||"
            << "globalGhostCellsPtr_||globalDeadCellsPtr_"
            << " already set"
            << abort(FatalError);
    }

    gammaCellTypePtr_=
    (
        new volScalarField
            (
                IOobject
                (
                    "gamma_Cut",
                    this->time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                *this,
                dimensionedScalar("1", dimless, 1)
            )
    );

    volScalarField& gammaCellType = *gammaCellTypePtr_;

    scalarField& gammaCellTypeI = gammaCellType.primitiveFieldRef();
    forAll(gammaCellTypeI,cellI)
    {
        if(cellInfoList[cellI]==cellClassification::OUTSIDE) gammaCellTypeI[cellI] = 1;
        else if(cellInfoList[cellI]==cellClassification::CUT) gammaCellTypeI[cellI] = 0;
        else gammaCellTypeI[cellI] = -1;
    }

    // Evaluate the coupled patchField, to replace evaluateCoupled()
    evaluateCoupled(gammaCellType);
    evaluateUnCoupled(gammaCellType);

    // calculate global gamma live cells (1), IB cells (0.5), ghost cells (-0.5), and dead cells(-1)

    labelHashSet globalIbCellsSet;
    labelHashSet globalGhostCellsSet;
    labelHashSet globalLiveCellsSet;
    labelHashSet globalDeadCellsSet;

    List<labelHashSet> ibCellsSetList(objectsList().size());
    List<labelHashSet> ghostCellsSetList(objectsList().size());

    forAll(objectsList(),objectID)
    {
        labelHashSet ibCellsSet(ibCellsList()[objectID]);
        ibCellsSetList[objectID]=ibCellsSet;
        labelHashSet ghostCellsSet(ghostCellsList()[objectID]);
        ghostCellsSetList[objectID]=ghostCellsSet;
    }
    
    // Find shared IB/ghost cells between different IB objects
    // The current code does not take care of the overlaps
    // between each IB object.
    forAll(objectsList(),objectID)
    {
        forAll(ibCellsList()[objectID], I)
        {
            label cellID = ibCellsList()[objectID][I];
            if(!globalIbCellsSet.found(cellID))
            {
                globalIbCellsSet.insert(cellID);
            }
            else
            {
                labelList sharedID(0);
                forAll(objectsList(),I)
                {
                    if(I!=objectID and ibCellsSetList[I].found(cellID)) sharedID.append(I);
                }
       
            }
        }

        forAll(ghostCellsList()[objectID], I)
        {
            label cellID = ghostCellsList()[objectID][I];
            if(!globalGhostCellsSet.found(cellID))
            {
                globalGhostCellsSet.insert(cellID);
            }
            else
            {
                labelList sharedID(0);
                forAll(objectsList(),I)
                {
                    if(I!=objectID and ghostCellsSetList[I].found(cellID)) sharedID.append(I);
                }
            
            }
        }
    }
    
    forAll(objectsList(),objectID)
    {
        forAll(ibDeadCellsList()[objectID], I)
        {
            label cellID = ibDeadCellsList()[objectID][I];
            if(!globalDeadCellsSet.found(cellID))// and !globalIbCellsSet.found(cellID))
            {
                globalDeadCellsSet.insert(cellID);
            }
        }
    }
    
    forAll(objectsList(),objectID)
    {
        forAll(ibLiveCellsList()[objectID], I)
        {
            label cellID = ibLiveCellsList()[objectID][I];
            if(!globalLiveCellsSet.found(cellID) 
                and !globalIbCellsSet.found(cellID)
                and !globalDeadCellsSet.found(cellID))
            {
                globalLiveCellsSet.insert(cellID);
            }
        }
    }

    globalLiveCellsPtr_ = new labelList(globalLiveCellsSet.toc());
    globalIbCellsPtr_ = new labelList(globalIbCellsSet.toc());
    globalGhostCellsPtr_ = new labelList(globalGhostCellsSet.toc());
    globalDeadCellsPtr_ = new labelList(globalDeadCellsSet.toc());

    gammaPtr_=
    (
        new volScalarField
            (
                IOobject
                (
                    "gamma",
                    this->time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                *this,
                dimensionedScalar("0", dimless, 0)
            )
    );

    volScalarField& gamma = *gammaPtr_;

    scalarField& gammaI = gamma.primitiveFieldRef();

    forAll(gammaI,cellI)
    {
        if(globalLiveCellsSet.found(cellI)) gammaI[cellI] = 1;
        else if(globalIbCellsSet.found(cellI)) gammaI[cellI] = 0.5;
        else if(globalDeadCellsSet.found(cellI)) gammaI[cellI] = -1;
        else gammaI[cellI] = -1;
    }

    // Evaluate the coupled patchField, to replace evaluateCoupled()
    evaluateCoupled(gamma);
    evaluateUnCoupled(gamma);

    volScalarField& Gamma = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("Gamma"));
    Gamma == gamma;
}

// make face gamma information
void Foam::immersedBoundaryFvMesh::makeSGamma
(
    const label& objectID
) const
{
    sGammaListPtr_->set
    (
        objectID,
        new surfaceScalarField
        (
            IOobject
            (
                "sGamma_"  + objectNames(objectID),
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            *this,
            dimensionedScalar("zero", dimless, 0),
            calculatedFvsPatchScalarField::typeName
         )
    );
       PtrList<surfaceScalarField>& sGammaList = *sGammaListPtr_;
       surfaceScalarField& sGamma = sGammaList[objectID];
    // Get access to components of sGamma
    scalarField& sGammaI = sGamma.primitiveFieldRef();
    surfaceScalarField::Boundary& sGammaPatches =
        sGamma.boundaryFieldRef();

    const unallocLabelList& owner = this->owner();
    const unallocLabelList& neighbour = this->neighbour();

    // Live cells indicator
    const volScalarField& gExt = ibGammaList()[objectID];
    const scalarField& gExtIn = gExt.primitiveField();

    const volScalarField::Boundary& gExtPatches = gExt.boundaryField();
    // Internal faces: flux is live between all active and IB cells
    forAll (sGammaI, faceI)
    {
        // If both cells are live, flux is live
        if
        (
            gExtIn[owner[faceI]] >= 0
         && gExtIn[neighbour[faceI]] >=0
        )
        {
            sGammaI[faceI] = 1.0;
        }
    }
    // Kill fluxes between two IB cells
    forAll (sGammaI, faceI)
    {
        if
        (
            gExtIn[owner[faceI]] == 0
         && gExtIn[neighbour[faceI]] == 0
        )
        {
            sGammaI[faceI] = 0.0;
        }
    }
    forAll (gExtPatches, patchI)
    {
        if (gExtPatches[patchI].coupled())
        {
            scalarField& gP = sGammaPatches[patchI];

            // For coupled patches, check gamma
            scalarField gammaOwn = gExtPatches[patchI].patchInternalField();

            scalarField gammaNei = gExtPatches[patchI].patchNeighbourField();

            forAll (gammaOwn, faceI)
            {
                if
                (
                    gammaOwn[faceI] >= 0
                 && gammaNei[faceI] >=0
                )
                {
                    gP[faceI] = 1.0;
                }
            }

            // For coupled patches, kill IB
            scalarField gammaIbOwn = gExtPatches[patchI].patchInternalField();

            scalarField gammaIbNei = gExtPatches[patchI].patchNeighbourField();

            forAll (gammaOwn, faceI)
            {
                if
                (
                    gammaIbOwn[faceI] == 0
                 && gammaIbNei[faceI] == 0
                )
                {
                    gP[faceI] = 0.0;
                }
            }
        }
        else
        {
            // For regular patches, check live cells only to achieve
            // correct global mass adjustment.
            // HJ, 21/May/2012
            scalarField gammaFc =
                gExtPatches[patchI].patchInternalField();

            scalarField& gP = sGammaPatches[patchI];
            
            const labelList& fc = this->boundary()[patchI].faceCells();
            
            forAll (gammaFc, faceI)
            {
                label cellID = fc[faceI];
             
               
                if (gExtIn[cellID] > 0) // cut cell flux should be zero everywhere on boundary
                {
                    gP[faceI] = 1.0;
                }
                else
                {
                    gP[faceI] = 0.0;
                }
            }
        }
    }
}

void Foam::immersedBoundaryFvMesh::writeibCell() const
{
    forAll(objectsList(), objectID)
    {
        volScalarField ibCell
        (
            IOobject
            (
                "ibCell_"+ibProperties().subDict("objects").toc()[objectID],
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            *this,
            dimensionedScalar("1", dimless, 1)
        );
        ibCell=gammaCellTypeList()[objectID];

        scalarField& ibCellI=ibCell.primitiveFieldRef();

        const labelList& ibCells=ibCellsList()[objectID];

        forAll(ibCellI, cellID)
        {
            if(ibCellI[cellID]==0)ibCellI[cellID]=-1;
        }

        forAll(ibCells, I)
        {
            label cellID = ibCells[I];
            ibCellI[cellID]=0;
        }
        ibCell.write();
        ibCell.clear();
    }
}


void Foam::immersedBoundaryFvMesh::writeGhostCell() const
{
    forAll(objectsList(), objectID)
    {
        volScalarField ghostCell
        (
            IOobject
            (
                "ghostCell_"+ibProperties().subDict("objects").toc()[objectID],
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            *this,
            dimensionedScalar("1", dimless, 1)
        );
        ghostCell=gammaCellTypeList()[objectID];

        scalarField& ghostCellI=ghostCell.primitiveFieldRef();

        const labelList& ghostCells=ghostCellsList()[objectID];

        forAll(ghostCellI, cellID)
        {
            if(ghostCellI[cellID]==0)ghostCellI[cellID]=1;
        }

        forAll(ghostCells, I)
        {
            label cellID = ghostCells[I];
            ghostCellI[cellID]=0;
        }
        ghostCell.write();
        ghostCell.clear();
    }
}

template<class Type>
void Foam::immersedBoundaryFvMesh::evaluateCoupled
(
    GeometricField<Type, fvPatchField, volMesh>& volValues
)const
{
    forAll (volValues.boundaryFieldRef(), patchI)
    {

        if (volValues.boundaryFieldRef()[patchI].coupled())
        {
            volValues.boundaryFieldRef()[patchI].initEvaluate
                (
                    Pstream::commsTypes::blocking
                );
        }

        if (volValues.boundaryFieldRef()[patchI].coupled())
        {
            volValues.boundaryFieldRef()[patchI].evaluate
                (
                    Pstream::commsTypes::blocking
                );
        }
    }
}

template<class Type>
void Foam::immersedBoundaryFvMesh::evaluateUnCoupled
(
    GeometricField<Type, fvPatchField, volMesh>& volValues
)const
{
    forAll (volValues.boundaryFieldRef(), patchI)
    {
        if (!volValues.boundaryFieldRef()[patchI].coupled())
        {
            Field<Type> PIF = volValues.boundaryField()[patchI].patchInternalField();
            forAll(volValues.boundaryFieldRef()[patchI],I)
            {
                volValues.boundaryFieldRef()[patchI][I]=PIF[I];
            }
        }
    }
}
// ************************************************************************* //

