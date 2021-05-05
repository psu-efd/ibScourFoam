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

#include "immersedBoundaryStencils.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "vectorTools.H"
#include "simpleObjectRegistry.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
 defineTypeNameAndDebug(immersedBoundaryStencils, 0);
 //addToRunTimeSelectionTable(searchableSurface, immersedBoundaryStencils, dict);
 //word immersedBoundaryStencils::meshSubDir = "triSurface";
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryStencils::immersedBoundaryStencils
(
    const word& name,
    const pointField& points,
    const labelList& cellIDList,
    const volScalarField& gamma,
    const meshSearch& queryMesh,
    const dictionary& dict
)
:
mesh_(gamma.mesh()),
queryMesh_(queryMesh),
gamma_(gamma),
points_(points),
cellIDList_(cellIDList),
namePtr_(new word(name)),
cellCellsPtr_(nullptr),
procCentresPtr_(nullptr),
procGammaPtr_(nullptr),
cellProcCellsPtr_(nullptr),
cellWeightsPtr_(nullptr),
cellProcWeightsPtr_(nullptr),
dictPtr_(new dictionary(dict)),
maxCellCellRows_(dict.lookupOrDefault<scalar>("maxCellCellRows",4.0)),
radiusFactor_(dict.lookupOrDefault<scalar>("radiusFactor",3.0)),
gammaInclude_c_(dict.lookupOrDefault<scalar>("gammaInclude_c",-0.5)),
gammaWeight_c_(dict.lookupOrDefault<scalar>("gammaWeight_c",0.))
{
Info <<" set immersedBoundaryStencils "<<endl;
    if(cellIDList.size()!=points.size())
    {
        FatalErrorIn("immersedBoundaryStencils::immersedBoundaryStencils() const")
            << "size of cellIDList "<<cellIDList.size()
            <<" does not equal to "
            << "size of points "<<points.size()
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::immersedBoundaryStencils::~immersedBoundaryStencils()
{
    clearOut();
}


void Foam::immersedBoundaryStencils::clearOut()
{
    deleteDemandDrivenData(namePtr_);
    deleteDemandDrivenData(cellCellsPtr_);
    deleteDemandDrivenData(procCentresPtr_);
    deleteDemandDrivenData(procGammaPtr_);
    deleteDemandDrivenData(cellProcCellsPtr_);
    deleteDemandDrivenData(procCellsPtr_);
    deleteDemandDrivenData(cellWeightsPtr_);
    deleteDemandDrivenData(cellProcWeightsPtr_);
    deleteDemandDrivenData(dictPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::immersedBoundaryStencils::makeStencils() const
{
    cellCellsPtr_ = new labelListList(points().size());
    labelListList& cellCells = *cellCellsPtr_;

    procCentresPtr_ = new vectorListList(Pstream::nProcs());
    vectorListList& procCentres = *procCentresPtr_;

    procGammaPtr_ = new scalarListList(Pstream::nProcs());
    scalarListList& procGamma = *procGammaPtr_;

    cellProcCellsPtr_ = new List<List<labelPair> > (points().size());
    List<List<labelPair> >& cellProcCells = * cellProcCellsPtr_;

    procCellsPtr_ = new labelListList(Pstream::nProcs());
    labelListList& procCells = * procCellsPtr_;

    const cellList& meshCells = mesh().cells();
    const vectorField& C = mesh().cellCentres();
    scalarField rM(points().size(),0);
    const labelList& ibc=cellIDList();
    forAll(rM,I)
    {
        rM[I]=cellSize(ibc[I]);
    }
    rM *= radiusFactor_;

    forAll (cellCells, cellI)
    {
        labelList curCells;

        findCellCells
        (
            C[ibc[cellI]],
            ibc[cellI],
            curCells
        );
        cellCells[cellI] = labelList(curCells.size(), -1);

        label cI = 0;

        for (label i = 0; i < curCells.size(); i++)
        {
            label curCell = curCells[i];
            scalar r = mag(C[curCell] - C[ibc[cellI]]);

            if (r <= rM[cellI])
            {

                cellCells[cellI][cI++] = curCell;
            }
        }

        cellCells[cellI].setSize(cI);

        if(cI<1 and debug==2)
        {
           Info<<ibc[cellI]<<endl;
        }

        //cellCells[cellI].setSize(min(cI,maxCellCellNum_)); // to limit the number of cellCells
    }
    if (Pstream::parRun())
    {
        // Find immersed boundary cells whose cellCells next to processor boundaries
        labelHashSet nextProcIbCellsSet;
        forAll (ibc, cellI)
        {
            labelList curCellCells = cellCells[cellI];

            //curCellCells.append(ibc[cellI]);// need to include ibc cell itself

            if (curCellCells.size())
            {
                forAll (curCellCells, cI)
                {
                    const labelList& faces = meshCells[curCellCells[cI]];

                    bool foundProcessorFace = false;

                    forAll (faces, faceI)
                    {
                        label patchID =
                            mesh().boundaryMesh().whichPatch(faces[faceI]);

                        if (patchID != -1)
                        {
                            if
                            (
                                isA<processorPolyPatch>
                                (
                                    mesh().boundaryMesh()[patchID]
                                )
                            )
                            {
                                foundProcessorFace = true;
                            }
                        }
                    }

                    if (foundProcessorFace)
                    {
                        nextProcIbCellsSet.insert(cellI);
                        break;
                    }
                }
            }
            else
            {
                const labelList& faces = meshCells[ibc[cellI]];

                bool foundProcessorFace = false;

                forAll (faces, faceI)
                {
                    label patchID =
                        mesh().boundaryMesh().whichPatch(faces[faceI]);

                    if (patchID != -1)
                    {
                        if
                        (
                            isA<processorPolyPatch>
                            (
                                mesh().boundaryMesh()[patchID]
                            )
                        )
                        {
                            foundProcessorFace = true;
                        }
                    }
                }

                if (foundProcessorFace)
                {
                    nextProcIbCellsSet.insert(cellI);
                }
            }
        }
        labelList nextProcIbCells = nextProcIbCellsSet.toc();

        sort(nextProcIbCells);

        // Send and receive ibc centres and radii
        vectorListList ctrs(Pstream::nProcs());

        ctrs[Pstream::myProcNo()].setSize(nextProcIbCells.size());
        vectorList& centres = ctrs[Pstream::myProcNo()];

        forAll (centres, cellI)
        {
            centres[cellI] = C[ibc[nextProcIbCells[cellI]]];
        }
        //centres = vectorList(C,ibc);

        Pstream::gatherList(ctrs);
        Pstream::scatterList(ctrs);

        scalarListList rMax(Pstream::nProcs());

        rMax[Pstream::myProcNo()] = scalarField(rM);

        Pstream::gatherList(rMax);
        Pstream::scatterList(rMax);

        // Find cells needed by other processors

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {

            labelHashSet procCellSet;

            if (procI != Pstream::myProcNo())
            {

                forAll (ctrs[procI], cellI)
                {
                  label nearestCellID =
                        //findNearestCell(ctrs[procI][cellI]);
                        queryMesh().findNearestCell(ctrs[procI][cellI], -1, true);

                    if (nearestCellID == -1)
                    {
                        FatalErrorIn
                        (
                            "immersedBoundaryFvMesh::makeIbCellCells() const"
                        ) << "Can't find nearest cell."
                            << abort(FatalError);
                    }

                    scalar R = mag(C[nearestCellID] - ctrs[procI][cellI]);

                    if (R < rMax[procI][cellI])
                    {
                        if (!procCellSet.found(nearestCellID))
                        {
                            procCellSet.insert(nearestCellID);
                        }

                        labelList tmpCellList;

                        findCellCells
                        (
                            ctrs[procI][cellI],
                            nearestCellID,
                            tmpCellList
                        );

                        forAll (tmpCellList, cI)
                        {
                            scalar r =
                                mag
                                (
                                    C[tmpCellList[cI]]
                                  - ctrs[procI][cellI]
                                );

                            if (r <= rMax[procI][cellI])
                            {
                                if (!procCellSet.found(tmpCellList[cI]))
                                {
                                    procCellSet.insert(tmpCellList[cI]);
                                }
                            }
                        }

                    }
                }
            }
            procCells[procI] = procCellSet.toc();
        }

        // Send and receive sizes
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::commsTypes::blocking,
                        procI,
                        sizeof(label)
                    );

                    toProc << procCells[procI].size();
                }
            }
        }

        labelList procSizes(Pstream::nProcs(), 0);
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::commsTypes::blocking,
                        procI,
                        sizeof(label)
                    );

                    fromProc >> procSizes[procI];
                }
            }
        }

        // Send cell centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                vectorField centres(C, procCells[procI]);

                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::commsTypes::blocking,
                        procI,
                        centres.size()*sizeof(vector)
                    );

                    toProc << centres;
                }
            }
        }

        // Receive cell centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                procCentres[procI].setSize(procSizes[procI]);

                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::commsTypes::blocking,
                        procI,
                        procSizes[procI]*sizeof(vector)
                    );

                    fromProc >> procCentres[procI];
                }
            }
            // else: already set to zero-size field
        }
        scalarField gammaI = gamma().primitiveField();
        // Send cell gamma
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                scalarField gamma(gammaI, procCells[procI]);

                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::commsTypes::blocking,
                        procI,
                        gamma.size()*sizeof(scalar)
                    );

                    toProc << gamma;
                }
            }
        }

        // Receive cell gamma
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                procGamma[procI].setSize(procSizes[procI]);

                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::commsTypes::blocking,
                        procI,
                        procSizes[procI]*sizeof(scalar)
                    );

                    fromProc >> procGamma[procI];
                }
            }
        }

        // Cell-procCells addressing
        forAll (cellProcCells, cellI)
        {
            scalar rMax = rM[cellI];

            cellProcCells[cellI].setSize(1000);
            //cellProcCells[cellI].clear();

            label index = 0;
            forAll (procCentres, procI)
            {
                if(procI!=Pstream::myProcNo())
                {
                    forAll (procCentres[procI], pointI)
                    {
                        scalar r =
                            mag
                            (
                                procCentres[procI][pointI]
                              - C[ibc[cellI]]
                            );
                        if (r <= rMax)
                        {
                            cellProcCells[cellI][index].first() = procI;
                            cellProcCells[cellI][index].second() = pointI;
                            index++;
                        }
                    }
                }
            }
            cellProcCells[cellI].setSize(index);
        }
//        // to sort all cellProcCells and cellCells
//        scalarField maxDistance(ibc.size(),0);

//        forAll (cellProcCells, cellI)
//        {
//            scalarField distances(cellProcCells[cellI].size()+cellCells[cellI].size(),0);

//            label totalIndex = 0;

//            forAll (cellCells[cellI], index)
//            {

//                scalar r =
//                    mag
//                    (
//                        C[cellCells[cellI][index]]
//                      - C[ibc[cellI]]
//                    );
//                distances[totalIndex] = r;
//                totalIndex++;
//            }

//            forAll (cellProcCells[cellI], index)
//            {
//                label procI = cellProcCells[cellI][index].first();
//                label pointI = cellProcCells[cellI][index].second();
//                scalar r =
//                    mag
//                    (
//                        procCentres[procI][pointI]
//                      - C[ibc[cellI]]
//                    );
//                distances[totalIndex] = r;
//                totalIndex++;
//            }

//            SortableList<scalar> sortedDistances(distances);
//            //maxDistance[cellI] = sortedDistances[min(maxCellCellNum_,distances.size())-1];
//            maxDistance[cellI] = sortedDistances[distances.size()-1];
//            if(distances.size()<maxCellCellNum_)
//            {
//                //Pout<<"maxCellCellRows and radiusFactor needs to be larger"<<endl;
//            }
//        }
//
//        // make sure every cellCell and cellProcCells is smaller than maxDistance
//        forAll (ibc, cellI)
//        {
//            label newIndex = 0;

//            forAll (cellCells[cellI], index)
//            {

//                scalar r =
//                    mag
//                    (
//                        C[cellCells[cellI][index]]
//                      - C[ibc[cellI]]
//                    );

//                if (r <= maxDistance[cellI])
//                {
//                    newIndex++;
//                }
//            }
//            cellCells[cellI].setSize(newIndex);
//            newIndex = 0;
//            forAll (cellProcCells[cellI], index)
//            {
//                label procI = cellProcCells[cellI][index].first();
//                label pointI = cellProcCells[cellI][index].second();
//                scalar r =
//                    mag
//                    (
//                        procCentres[procI][pointI]
//                      - C[ibc[cellI]]
//                    );
//                if (r <= maxDistance[cellI])
//                {
//                    newIndex++;
//                }
//            }
//            cellProcCells[cellI].setSize(newIndex);
//        }
    }
}



void Foam::immersedBoundaryStencils::findCellCells
(
    const vector& pt,
    const label cellID,
    labelList& cellCells
)const
{
    const labelListList& cellPoints = mesh().cellPoints();
    const labelListList& pointCells = mesh().pointCells();

    //const scalarField& gammaI = ibGammaList()[objectID].internalField();
    const scalarField& gammaI = gamma().internalField();
    labelHashSet cellSet;
    cellSet.insert(cellID);

    // First row

    const labelList& curCellPoints = cellPoints[cellID];

    forAll (curCellPoints, pointI)
    {
        label curPoint = curCellPoints[pointI];
        const labelList& curPointCells = pointCells[curPoint];

        forAll (curPointCells, cI)
        {
            // stencil has to be live cells //or IB cells
            if (gammaI[curPointCells[cI]] > gammaInclude_c())
            {
                if (!cellSet.found(curPointCells[cI]))
                {
                    cellSet.insert(curPointCells[cI]);
                }
            }
        }
    }

    labelList curCells = cellSet.toc();

     // Second and other rows
    for (label nRows = 1; nRows < maxCellCellRows_; nRows++)
    {
        curCells = cellSet.toc();

        forAll (curCells, cellI)
        {
            label curCell = curCells[cellI];
            const labelList& curCellPoints = cellPoints[curCell];

            forAll (curCellPoints, pointI)
            {
                label curPoint = curCellPoints[pointI];
                const labelList& curPointCells = pointCells[curPoint];

                forAll (curPointCells, cI)
                {
                    if (gammaI[curPointCells[cI]] > gammaInclude_c())
                    {

                        if (!cellSet.found(curPointCells[cI]))
                        {
                            cellSet.insert(curPointCells[cI]);
                        }
                    }
                }
            }
        }
    }
    // Erase current cell
    cellSet.erase(cellID);

    // Sorting cells
    const vectorField& C = mesh().cellCentres();

    curCells = cellSet.toc();
    scalarField distances(curCells.size(), 0);

    forAll (distances, cI)
    {
        distances[cI] =
            mag(C[curCells[cI]] - pt);
    }

    SortableList<scalar> sortedDistances(distances);

    labelList sortedCells(curCells.size(), -1);

    //labelList sortedCells(min(maxSortedCellCellNum_,curCells.size()), -1);

    for (label i = 0; i < sortedCells.size(); i++)
    {
        sortedCells[i] =
            curCells[sortedDistances.indices()[i]];
    }

    cellCells = sortedCells;
}

Foam::scalar Foam::immersedBoundaryStencils::cellSize(label cellID) const
{
    scalar delta;

    if (mesh().nGeometricD() == 3)
    {
        delta = Foam::pow(mesh().V().field()[cellID], 1.0/3.0);
    }
    else
    {
        scalar thickness = 0.0;
        const Vector<label>& directions = mesh().geometricD();
        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh().bounds().span()[dir];
                break;
            }
        }

        delta = Foam::sqrt(mesh().V().field()[cellID]/thickness);
//        delta = Foam::pow(mesh().V().field()[cellID], 1.0/3.0);

    }
    return delta;
}

void Foam::immersedBoundaryStencils::makeStencilsWeights() const
{
    cellWeightsPtr_ = new scalarListList(points().size());
    scalarListList& cellWeights = * cellWeightsPtr_;

    cellProcWeightsPtr_ = new scalarListList(points().size());
    scalarListList& cellProcWeights = *cellProcWeightsPtr_;

    const labelListList& ibcc = cellCells();
    const List<List<labelPair> >& ibcProcC = cellProcCells();
    const vectorField& points = this->points();
    const vectorListList& CProc = procCentres();
    const scalarListList& gammaProc = procGamma();

    if(debug)
    {

        scalarField nibcc(ibcc.size());
        scalarField nibcc1(ibcc.size());
        forAll(ibcc,I)
        {
            scalar n1(ibcc[I].size()+ibcProcC[I].size());
            scalar n2(ibcProcC[I].size());

            nibcc[I]=n1;
            nibcc1[I]=n2;
        }
        Info<<"Num. of stencil cells for each "<<name()<<endl;
        Info<<"Min Max Average"<<endl;
        Info<<gMin(nibcc)<<" "<<gMax(nibcc)<<" "<<gAverage(nibcc)<<endl;

        if (Pstream::parRun())//Start of mpi run
        {
            Info<<"Num. of stencil cells for each "<<name()<<endl;
            Info<<"Min Max Average"<<endl;
            Info<<gMin(nibcc1)<<" "<<gMax(nibcc1)<<" "<<gAverage(nibcc1)<<endl;
        }
    }
    forAll (cellWeights, cellI)
    {
        cellWeights[cellI].setSize(ibcc[cellI].size(), 0);
    }

    forAll (cellProcWeights, cellI)
    {
        cellProcWeights[cellI].setSize(ibcProcC[cellI].size(), 0);
    }

    //const scalarField& gammaIn = ibGammaList()[objectID].internalField();
    const scalarField& gammaIn = gamma().internalField();

    const vectorField& CIn = mesh().C().internalField();

    // count Insufficient live neighbourhood points
    scalar counter=0;
    // Go through all cellCells and calculate inverse distance for
    // all live points
    forAll (points, cellI)
    {
        const vector& curP = points[cellI];

        scalar sumW = 0;

        // Local weights
        scalarList& curCW = cellWeights[cellI];

        const labelList& curCells = ibcc[cellI];
        forAll (curCells, ccI)
        {
            // Only pick live cells
            if (gammaIn[curCells[ccI]] >gammaWeight_c())
            {
                curCW[ccI] = 1.0/(SMALL+mag(CIn[curCells[ccI]] - curP));
                curCW[ccI] = pow(curCW[ccI],1);

                sumW += curCW[ccI];
            }
            else
            {
                curCW[ccI] = 0;
            }
        }
        // Processor weights
        const List<labelPair>& interpProcCells = ibcProcC[cellI];

        scalarList& curProcCW = cellProcWeights[cellI];

        forAll (interpProcCells, cProcI)
        {
            if
            (
                gammaProc
                [
                    interpProcCells[cProcI].first()
                ]
                [
                    interpProcCells[cProcI].second()
                ] > gammaWeight_c()
            )
            {
                curProcCW[cProcI] =
                    1/(mag
                        (
                            CProc
                            [
                                interpProcCells[cProcI].first()
                            ]
                            [
                                interpProcCells[cProcI].second()
                            ] - curP
                        )+SMALL);

                    curProcCW[cProcI]=pow(curProcCW[cProcI],1);

                sumW += curProcCW[cProcI];
            }
            else
            {
                curProcCW[cProcI] = 0;
            }
        }
        // Divide through by the sum
        if (sumW < SMALL and debug==2)
        {
            counter++;
            InfoIn
            (
                "void immersedBoundaryStencils::makeStencilsWeights()"
            )   << "Insufficient live neighbourhood for "<<name()<<" point "
                << points[cellI] << "." << nl
                << ibcc[cellI].size()<<tab
                << "Please adjust radiusFactor, distFactor or "
                << "maxCellCellRows "
                << "in immersedBoundaryProperties."
                << endl;

            // Reset sum and weights and use all points
            sumW = 0;
            curCW = 0;

            forAll (curCells, ccI)
            {
                // Use all cells
                curCW[ccI] = 1.0/(SMALL+mag(CIn[curCells[ccI]] - curP));

                curCW[ccI] = pow(curCW[ccI],1);
                sumW += curCW[ccI];
            }
        }

        if(sumW<SMALL)
        {
            curCW=0.0;
            curProcCW=0.0;
        }
        else
        {
            forAll (curCells, ccI)
            {
                curCW[ccI] /= sumW;
            }

            forAll (curProcCW, cProcI)
            {
                curProcCW[cProcI] /= sumW;
            }
        }

    }

    if(debug)
    {
        Info<<"Number of points with insufficient live neighbourhood: "<<returnReduce(counter, sumOp<label>())<<endl;
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
