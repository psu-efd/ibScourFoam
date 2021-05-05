/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include "meshSearch.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "cellClassification.H"
#include "wallPolyPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        immersedBoundaryFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::immersedBoundaryFvMesh::readDict()
{
    Info<<"Read constant/immersedBoundaryProperties"<<endl;

    if (debug==2)
    {
        InfoIn("void immersedBoundaryFvMesh::readDict() const")
            << "read immersedBoundaryProperties "
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (objectDictListPtr_ || ibPropertiesPtr_ ||IBtypeListPtr_)
    {
        FatalErrorIn("immersedBoundaryFvMesh::readDict() const")
            << "read immersedBoundaryProperties  "
            << "objectDictListPtr_ || ibPropertiesPtr_||IBtypeListPtr_"
            << abort(FatalError);
    }

    ibPropertiesPtr_=
        new IOdictionary(
            IOobject
            (
                "immersedBoundaryProperties",
                time().constant(),
                time(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

    // generate list of object names from immersedBoundaryProperties
    const wordList objectNames(ibProperties().subDict("objects").toc());
   
    objectDictListPtr_ = new PtrList<dictionary>(objectNames.size());
    PtrList<dictionary>& objectDictList = *objectDictListPtr_;
    
    IBtypeListPtr_ = new PtrList<word>(objectNames.size());
    
    forAll(objectNames,objectID)
    {

        // initialize sub-dict for each object
        Info<<"Object "<<objectID<<": "<<objectNames[objectID];
        objectDictList.set
        (
            objectID,
            new dictionary
            (
                ibProperties().subDict("objects").subDict(objectNames[objectID])
            )
        );
        
            // initialize IBtype for each object
        // classic or mix or ghost-cell
        word IBtype = ibProperties()
                        .subDict("objects")
                        .subDict(objectNames[objectID])
                        .lookupOrDefault<word>("IBtype","classic");

        IBtypeListPtr_->set
        (
            objectID,
            new word
            (
                IBtype
            )
        );

        if(IBtype == "ghost-cell")
        {

        }
        else if(IBtype == "classic")
        {

        }
        else if(IBtype == "mix")
        {

        }
        else
        {
            FatalErrorIn("immersedBoundaryFvMesh::markCells")
                << "IBtype " << IBtype
                << " is not ghost-cell or classic"
                << exit(FatalError);
        }

        Info<<", Uses "<<IBtype<<" method"<<endl;
        
    }
    
    


}






bool Foam::immersedBoundaryFvMesh::readStl()
{
    if(!readStl_)
    {
        return false;
    }
    
    
    Info<<"Read constant/triSurfaces"<<endl;

    if (debug==2)
    {
        InfoIn("void immersedBoundaryFvMesh::readStl() const")
            << "read triSurfaces "
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if ( objectsListPtr_ ||addObjectsListPtr_)
    {
        FatalErrorIn("immersedBoundaryFvMesh::readStl() const")
            << "read triSurfaces  "
            << "objectsListPtr_||addObjectsListPtr_"
            << abort(FatalError);
    }


    const wordList objectNames(ibProperties().subDict("objects").toc());
    objectsListPtr_ = new PtrList<triSurface>(objectDictList().size());
    PtrList<triSurface>& objectsList = *objectsListPtr_;
    
    addObjectsListPtr_ = new PtrList<triSurface>(objectDictList().size()-1);
    PtrList<triSurface>& addObjectsList = *addObjectsListPtr_;
    
    addObjectsNameListPtr_ = new PtrList<word>(objectDictList().size()-1);
    PtrList<word>& addObjectsNameList = *addObjectsNameListPtr_;
    
    addObjectsIDListPtr_ = new PtrList<label>(objectDictList().size()-1);
    PtrList<label>& addObjectsIDList = *addObjectsIDListPtr_;
    
    


    Info<<"Read objects surface mesh: "<<nl;
    
    label nExtraSurface = 0;
    forAll(objectDictList(),objectID)
    {

        // initialize triSurface for each object
        objectsList.set
        (
            objectID,
            new triSurfaceMesh
            (
                IOobject
                (
                    objectNames[objectID]  + ".stl",
                    mesh().time().constant(), // instance
                    "triSurface",                // local
                    mesh(),                   // registry
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
       
        bool multipleObject_ =  ibProperties().lookupOrDefault<bool>("multipleObject",false);

 	    if(multipleObject_)
        {       
            if( !objectDictList()[objectID].found("sediment"))
            {
                
                addObjectsList.set
                (
                    nExtraSurface,
                    new triSurface(objectsList[objectID])
                );
                addObjectsNameList.set
                (
                    nExtraSurface,
                    new word(objectNames[objectID])
                );
                addObjectsIDList.set
                (
                    nExtraSurface,
                    new label(objectID)
                );
                nExtraSurface++;
                
                Info<<"add   Object "<<objectID<<": "<<objectNames[objectID];
            }
        }

         
    }   
        
    readStl_ = false;
    return true;
}

void Foam::immersedBoundaryFvMesh::readDualDict()
{
    Info<<"Read constant/immersedBoundaryProperties"<<endl;

    if (debug==2)
    {
        InfoIn("void immersedBoundaryFvMesh::readDualDict() const")
            << "read immersedBoundaryProperties "
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if ( triSurfaceSearchListPtr_ ||dualListPtr_)
    {
        FatalErrorIn("immersedBoundaryFvMesh::readStl() const")
            << "read triSurfaces  "
            << "triSurfaceSearchListPtr_  ||dualListPtr_"
            << abort(FatalError);
    }

    
    
    const wordList objectNames(ibProperties().subDict("objects").toc());
    triSurfaceSearchListPtr_ = new PtrList<triSurfaceSearch>(objectDictList().size());
    PtrList<triSurfaceSearch>& triSurfaceSearchList = *triSurfaceSearchListPtr_; 
          
    
    dualListPtr_ = new PtrList<bool>(objectDictList().size());
    PtrList<bool>& dualList = *dualListPtr_;    
    

    forAll(objectDictList(),objectID)
    {
        triSurfaceSearchList.set
        (
            objectID,
            new triSurfaceSearch(objectsList()[objectID])
        );  
       

        // check if use dual mesh
        dualList.set 
        (
            objectID,
            new bool(
                        ibProperties()
                        .subDict("objects")
                        .subDict(objectNames[objectID])
                        .lookupOrDefault<bool>("dualMeshSwicth",false)
                    )
        );        

        

    }

    // initialize calculated values
    // variables transferred from hydro to IB (hit points based)
    wallShearStressListPtr_ =  new PtrList<vectorField>(objectNames.size());
    nutListPtr_ =  new PtrList<scalarField>(objectNames.size());
    depositionListPtr_ =  new PtrList<scalarField>(objectNames.size());
    entrainmentListPtr_ =  new PtrList<scalarField>(objectNames.size());
    accumDItaListPtr_ =  new PtrList<scalarField>(objectNames.size());
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryFvMesh::immersedBoundaryFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    updateIB_(true),
    readStl_(true),
    globalMeshPtr_(nullptr),
    queryMeshPtr_(new meshSearch(*this)),
    ibPropertiesPtr_(nullptr),
    objectDictListPtr_(nullptr),
    IBtypeListPtr_(nullptr),
    triSurfaceMeshListPtr_(nullptr),
    objectsListPtr_(nullptr),
    addObjectsListPtr_(nullptr),
    addObjectsNameListPtr_(nullptr),
    addObjectsIDListPtr_(nullptr),
    sampledTriSurfaceMeshListPtr_(nullptr),
    cutTriSurfaceMeshListPtr_(nullptr),
    triSurfaceSearchListPtr_(nullptr),
    gammaCellTypeListPtr_(nullptr),
    gammaCellTypePtr_(nullptr),
    gammaPtr_(nullptr),
    globalLiveCellsPtr_(nullptr),
    globalIbCellsPtr_(nullptr),
    globalGhostCellsPtr_(nullptr),
    globalDeadCellsPtr_(nullptr),
    ibCellsListPtr_(nullptr),
    oldIbCellsListPtr_(nullptr),

    ibDeadCellsListPtr_(nullptr),
    oldIbDeadCellsListPtr_(nullptr),
    
    ibLiveCellsListPtr_(nullptr),
    oldIbLiveCellsListPtr_(nullptr),
    
    ibGammaListPtr_(nullptr),
    sGammaListPtr_(nullptr),
    ibFacesListPtr_(nullptr),
    ibFaceCellsListPtr_(nullptr),
    ibFaceFlipsListPtr_(nullptr),
    ibFaceHitPointsListPtr_(nullptr),
    ibFaceHitFacesListPtr_(nullptr),
    ibFaceSamplingPointsListPtr_(nullptr),
    ibHitPointsListPtr_(nullptr),
    ibHitFacesListPtr_(nullptr),
    ibTriNormalsFlipListPtr_(nullptr),
    samplingPointsListPtr_(nullptr),
    farSamplingPointsListPtr_(nullptr),

    ghostCellsListPtr_(nullptr),
    ghostHitPointsListPtr_(nullptr),
    ghostHitFacesListPtr_(nullptr),
    imagePointsListPtr_(nullptr),

    samplingStencilsListPtr_(nullptr),
    farSamplingStencilsListPtr_(nullptr),
    ibCellStencilsListPtr_(nullptr),
    ibFaceSamplingStencilsListPtr_(nullptr),
    imageStencilsListPtr_(nullptr),

    ibCellsToTriAddrListPtr_(nullptr),
    ibCellsToTriWeightsListPtr_(nullptr),
    ibCellsToTriMatrixListPtr_(nullptr),
    ghostCellsToTriAddrListPtr_(nullptr),
    ghostCellsToTriWeightsListPtr_(nullptr),
    ghostCellsToTriMatrixListPtr_(nullptr),
    ibCellsToTriEdgeAddrListPtr_(nullptr),
    ibCellsToTriEdgeWeightsListPtr_(nullptr),
    ghostCellsToTriEdgeAddrListPtr_(nullptr),
    ghostCellsToTriEdgeWeightsListPtr_(nullptr),
    proTriFacesInMeshListPtr_(nullptr),   
 

    dualMeshAddressingListPtr_(nullptr),

    oldToNewFcMapListPrt_(nullptr),
    oldToNewEdgeMapListPrt_(nullptr),
    newToOldEdgeMapListPrt_(nullptr),
    oldToNewDualEdgeMapListPrt_(nullptr),
    newToOldDualEdgeMapListPrt_(nullptr),

    wallShearStressListPtr_(nullptr),
    nutListPtr_(nullptr),
    depositionListPtr_(nullptr),
    entrainmentListPtr_(nullptr),
    accumDItaListPtr_(nullptr),
    dualMeshListPtr_(nullptr),
    dualListPtr_(nullptr)
{


}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundaryFvMesh::~immersedBoundaryFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::immersedBoundaryFvMesh::update()
{
    // if
    if(!updateIB())
    {
        return false;
    }

    // clear out pointers associated with previous IB information
    clearOut();

 

    Info<<endl<<"Begin updating IB information"<<endl;
    const double Oldtime1=time().elapsedCpuTime();
    // read sub-dict and surfaceMesh
    readDict();
    readStl();
    readDualDict();
    const double Oldtime2=time().elapsedCpuTime();

    // marking ghost cells and IB cells.
    // making gammaCellTypePtr_ and gammaCellTypeListPtr_
    markCells();
    const double Oldtime3=time().elapsedCpuTime();

    
    const double Oldtime4=time().elapsedCpuTime();

    // Make stencil info for IB cells, ghost cells.
    // In the current version, makeStencilsInfo only initializes
    // all pointer lists associated to interpolation stencil.
    // Corresponding interpolation stencil is stored in
    // samplingStencilsListPtr_, ibCellStencils, ibFaceSamplingStencils, imageStencils
    // which are based on sampling points, ib cell centers, ib face centers,
    // and image points (image of ghost cell), respectively.
    //
    makeStencilsInfo();
    const double Oldtime5=time().elapsedCpuTime();

    // make addressing from hydrodynamics domain to IB
    // hit points -> triangle cell center/edge
    // located in immersedBoundaryFvMeshMapping.C
    // contains both makeTriAddressing and makeTriEdgeAddressing
    makeTriAddressing();
    const double Oldtime6=time().elapsedCpuTime();

    // Generate 2D mesh for immersed boundary.
    // It includes cutTriSurfaceMesh to provide necessary mesh data

    makeDualMesh();
    const double Oldtime8=time().elapsedCpuTime();

    if(debug)
    {
        Info<<"readDict Executation Time = "<<Oldtime2-Oldtime1<< " s"<<endl;
        Info<<"markCellsExecutation Time = "<<Oldtime3-Oldtime2<< " s"<<endl;
        Info<<"makeStencilsInfo Executation Time = "<<Oldtime5-Oldtime4<< " s"<<endl;
        Info<<"makeTriAddressing Executation Time = "<<Oldtime6-Oldtime5<< " s"<<endl;
        
        Info<<"makeDualMesh Executation Time = "<<Oldtime8-Oldtime6<< " s"<<endl;

        Info<<"Total update Executation Time = "<<Oldtime8-Oldtime1<< " s"<<endl;
    }

    Info<<endl<<"****** Finish updating surface mesh information ******"<<endl<<endl;


    IBHasUpdated();

    return true;
}

void Foam::immersedBoundaryFvMesh::finalClearOut()
{
    Info<<"final clearOut"<<endl;
    
    deleteDemandDrivenData(objectsListPtr_);
    deleteDemandDrivenData(addObjectsListPtr_);
    deleteDemandDrivenData(addObjectsNameListPtr_);
    deleteDemandDrivenData(addObjectsIDListPtr_);
    deleteDemandDrivenData(globalMeshPtr_);

    deleteDemandDrivenData(queryMeshPtr_);

    clearOut();
}

void Foam::immersedBoundaryFvMesh::clearOut()
{
    Info<<"clearOut"<<endl;

    deleteDemandDrivenData(ibPropertiesPtr_);
    deleteDemandDrivenData(objectDictListPtr_);
    deleteDemandDrivenData(triSurfaceSearchListPtr_);
    deleteDemandDrivenData(IBtypeListPtr_);

    deleteDemandDrivenData(sampledTriSurfaceMeshListPtr_);
    deleteDemandDrivenData(cutTriSurfaceMeshListPtr_);
    deleteDemandDrivenData(triSurfaceMeshListPtr_);

    deleteDemandDrivenData(gammaCellTypeListPtr_);
    deleteDemandDrivenData(gammaCellTypePtr_);
    deleteDemandDrivenData(gammaPtr_);
    deleteDemandDrivenData(globalLiveCellsPtr_);
    deleteDemandDrivenData(globalIbCellsPtr_);
    deleteDemandDrivenData(globalGhostCellsPtr_);
    deleteDemandDrivenData(globalDeadCellsPtr_);
    // IB cells
    deleteDemandDrivenData(ibCellsListPtr_);
    deleteDemandDrivenData(oldIbCellsListPtr_);

    
    deleteDemandDrivenData(ibDeadCellsListPtr_);
    deleteDemandDrivenData(oldIbDeadCellsListPtr_);
    deleteDemandDrivenData(ibLiveCellsListPtr_);
    deleteDemandDrivenData(oldIbLiveCellsListPtr_);
    deleteDemandDrivenData(ibGammaListPtr_);
    deleteDemandDrivenData(sGammaListPtr_);
    deleteDemandDrivenData(ibFacesListPtr_);
    deleteDemandDrivenData(ibFaceCellsListPtr_);
    deleteDemandDrivenData(ibFaceFlipsListPtr_);
    deleteDemandDrivenData(ibFaceHitPointsListPtr_);
    deleteDemandDrivenData(ibFaceHitFacesListPtr_);
    deleteDemandDrivenData(ibFaceSamplingPointsListPtr_);
    deleteDemandDrivenData(ibHitPointsListPtr_);
    deleteDemandDrivenData(ibHitFacesListPtr_);
    deleteDemandDrivenData(ibTriNormalsFlipListPtr_);
    deleteDemandDrivenData(samplingPointsListPtr_);
    deleteDemandDrivenData(farSamplingPointsListPtr_);


    // ghost cells
    deleteDemandDrivenData(ghostCellsListPtr_);

    deleteDemandDrivenData(ghostHitPointsListPtr_);
    deleteDemandDrivenData(ghostHitFacesListPtr_);
    deleteDemandDrivenData(imagePointsListPtr_);

    deleteDemandDrivenData(samplingStencilsListPtr_);
    deleteDemandDrivenData(farSamplingStencilsListPtr_);
    deleteDemandDrivenData(ibCellStencilsListPtr_);
    deleteDemandDrivenData(ibFaceSamplingStencilsListPtr_);
    deleteDemandDrivenData(imageStencilsListPtr_);

    deleteDemandDrivenData(ibCellsToTriAddrListPtr_);
    deleteDemandDrivenData(ibCellsToTriWeightsListPtr_);
    deleteDemandDrivenData(ibCellsToTriMatrixListPtr_);
    deleteDemandDrivenData(ghostCellsToTriAddrListPtr_);
    deleteDemandDrivenData(ghostCellsToTriWeightsListPtr_);
    deleteDemandDrivenData(ghostCellsToTriMatrixListPtr_);
    deleteDemandDrivenData(ibCellsToTriEdgeAddrListPtr_);
    deleteDemandDrivenData(ibCellsToTriEdgeWeightsListPtr_);
    deleteDemandDrivenData(ghostCellsToTriEdgeAddrListPtr_);
    deleteDemandDrivenData(ghostCellsToTriEdgeWeightsListPtr_);
    deleteDemandDrivenData(proTriFacesInMeshListPtr_);
    deleteDemandDrivenData(oldToNewFcMapListPrt_);
    deleteDemandDrivenData(oldToNewEdgeMapListPrt_);
    deleteDemandDrivenData(newToOldEdgeMapListPrt_);
    deleteDemandDrivenData(oldToNewDualEdgeMapListPrt_);
    deleteDemandDrivenData(newToOldDualEdgeMapListPrt_);    

    // calculated values
    deleteDemandDrivenData(wallShearStressListPtr_);
    deleteDemandDrivenData(nutListPtr_);
    deleteDemandDrivenData(depositionListPtr_);
    deleteDemandDrivenData(entrainmentListPtr_);
    deleteDemandDrivenData(accumDItaListPtr_);      
    deleteDemandDrivenData(dualMeshListPtr_);
    deleteDemandDrivenData(dualListPtr_); 
    deleteDemandDrivenData(dualMeshAddressingListPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "immersedBoundaryFvMeshMarkCells.C"
#include "immersedBoundaryFvMeshStencils.C"
#include "immersedBoundaryFvMeshCorrectPhi.C"
#include "immersedBoundaryFvMeshUpdateCellValues.C"
#include "immersedBoundaryFvMeshPressureCorrector.C"
#include "immersedBoundaryFvMeshTurbulence.C"
#include "immersedBoundaryFvMeshMapping.C"
#include "immersedBoundaryFvMeshTriSurfaceTools.C"
#include "immersedBoundaryFvMeshDualMesh.C"
#include "immersedBoundaryFvMeshPostEvaluation.C"
#include "immersedBoundaryFvMeshEvaluateC.C"
#include "immersedBoundaryFvMeshDualCalcs.C"


// ************************************************************************* //
