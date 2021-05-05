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

Class
    Foam::immersedBoundaryFvMesh

Description
    Calculate IB information. It is based on dynamicFvMesh.

    It also provides functions for local scour modelling.




SourceFiles
    immersedBoundaryFvMesh.C
    immersedBoundaryFvMeshMarkCells.C
    immersedBoundaryFvMeshStencils.C
    immersedBoundaryFvMeshCorrectPhi.C
    immersedBoundaryFvMeshUpdateCellValues.C
    immersedBoundaryFvMeshPressureCorrector.C
    immersedBoundaryFvMeshTurbulence.C
    immersedBoundaryFvMeshMapping.C
    immersedBoundaryFvMeshTriSurfaceTools.C
    immersedBoundaryFvMeshDualMesh.C
    immersedBoundaryFvMeshPostEvaluation.C
    immersedBoundaryFvMeshEvaluateC.C

  
    
\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryFvMesh_H
#define immersedBoundaryFvMesh_H

#include "dynamicFvMesh.H"
#include "triSurfaceMesh.H"
#include "meshSearch.H"
#include "cellClassification.H"
#include "cellInfo.H"
#include "vectorList.H"
#include "FieldFields.H"
#include "fvMatrix.H"
#include "turbulenceModel.H"

#include "edgeIntersections.H"
#include "surfaceIntersection.H"
#include "booleanSurface.H"
#include "intersectedSurface.H"
#include "sampledTriSurfaceMesh.H"
#include "isoSurfaceCell.H"
#include "isoSurface.H"
#include "volPointInterpolation.H"
#include "cutTriSurfaceMesh.H"
#include "CompactIOList.H"
#include "PrimitivePatchInterpolation.H"
#include "mappedPolyPatch.H"
#include "globalIndex.H"
#include "IOmanip.H"
//#include "immersedBoundaryInterpolation.H"
//#include "immersedBoundaryFvMeshLduAddressing.H"
#include "immersedBoundaryStencils.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class immersedBoundaryFvMesh Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryFvMesh
:
    public dynamicFvMesh
{
public:

    // Public data types

        //- Type of cell.
        enum cellType
        {
            FLUID,
            IB,     // Inside of surface
            GHOST,    // Outside ,,
            SOLID         // cut by surface
        };

        typedef List<vectorList> vectorListList;
        typedef PrimitivePatch<face, ::Foam::List, pointField, point> PPatch;
        typedef PrimitivePatchInterpolation< PrimitivePatch<face, ::Foam::List, pointField, point> >  PPatchInterpolation;


private:

    // Private data

        //- if need to update surface mesh information clearOut() & update()
        //- need to change for multiple objects
        mutable bool updateIB_;
        
        //- if need to read stil
        mutable bool readStl_;

        //- global mesh used in MPI
        mutable fvMesh* globalMeshPtr_;

        //- mesh search
        mutable meshSearch* queryMeshPtr_;

        //- immersed boundary dictionary
        mutable IOdictionary* ibPropertiesPtr_;

        //- sub-dict for objects
        mutable PtrList<dictionary>* objectDictListPtr_;

        //- IB type  classic or ghost cell method
        mutable PtrList<word>* IBtypeListPtr_;

        //- surface mesh for objects
        mutable PtrList<triSurfaceMesh>* triSurfaceMeshListPtr_;

        //- stl for objects
        mutable PtrList<triSurface>* objectsListPtr_;
        
        mutable PtrList<triSurface>* addObjectsListPtr_;
        mutable PtrList<word>* addObjectsNameListPtr_;
        mutable PtrList<label>* addObjectsIDListPtr_;
        //- cut triSurface related
        //- new cut surface
        mutable PtrList<sampledTriSurfaceMesh>* sampledTriSurfaceMeshListPtr_;

        mutable PtrList<cutTriSurfaceMesh>* cutTriSurfaceMeshListPtr_;

        //- surface mesh search for objects
        mutable PtrList<triSurfaceSearch>* triSurfaceSearchListPtr_;

        //- indicator for fluid cells (1), cut cells (0), and solid cells(-1)
        mutable PtrList<volScalarField>* gammaCellTypeListPtr_;

        //- global gammaCellType fluid cells (1), cut cells (0), and solid cells(-1)
        mutable volScalarField* gammaCellTypePtr_;

        //- global gamma live cells (1), IB cells (0.5), ghost cells (-0.5), and dead cells(-1)
        mutable volScalarField* gammaPtr_;

        mutable labelList* globalLiveCellsPtr_;

        mutable labelList* globalIbCellsPtr_;

        mutable labelList* globalGhostCellsPtr_;

        mutable labelList* globalDeadCellsPtr_;

        //- IB cell list
        mutable PtrList<labelList>* ibCellsListPtr_;
        mutable PtrList<labelList>* oldIbCellsListPtr_;

        mutable PtrList<labelList>* ibDeadCellsListPtr_;
        mutable PtrList<labelList>* oldIbDeadCellsListPtr_;
        
        mutable PtrList<labelList>* ibLiveCellsListPtr_;
        mutable PtrList<labelList>* oldIbLiveCellsListPtr_;
        //- local gamma live cells (1), IB cells (0), and dead cells(-1)
        mutable  PtrList<volScalarField>* ibGammaListPtr_;
        mutable  PtrList<surfaceScalarField>* sGammaListPtr_;

        mutable PtrList<labelList>* ibFacesListPtr_;

        mutable PtrList<labelList>* ibFaceCellsListPtr_;

        mutable PtrList<List<bool>>* ibFaceFlipsListPtr_;

        mutable PtrList<pointField>* ibFaceHitPointsListPtr_;

        mutable PtrList<labelList>* ibFaceHitFacesListPtr_;

        mutable PtrList<pointField>* ibFaceSamplingPointsListPtr_;

        mutable PtrList<pointField>* ibHitPointsListPtr_;

        mutable PtrList<labelList>* ibHitFacesListPtr_;

        mutable PtrList<bool>* ibTriNormalsFlipListPtr_;

        mutable PtrList<pointField>* samplingPointsListPtr_;
        mutable PtrList<pointField>* farSamplingPointsListPtr_;
      

        //- ghost cell list
        mutable PtrList<labelList>* ghostCellsListPtr_;

        mutable PtrList<pointField>* ghostHitPointsListPtr_;

        mutable PtrList<labelList>* ghostHitFacesListPtr_;

        mutable PtrList<pointField>* imagePointsListPtr_;

        mutable PtrList<immersedBoundaryStencils>* samplingStencilsListPtr_;
        mutable PtrList<immersedBoundaryStencils>* farSamplingStencilsListPtr_;
        mutable PtrList<immersedBoundaryStencils>* ibCellStencilsListPtr_;

        mutable PtrList<immersedBoundaryStencils>* ibFaceSamplingStencilsListPtr_;

        mutable PtrList<immersedBoundaryStencils>* imageStencilsListPtr_;

        // tri-addressing
        mutable PtrList<labelListList>* ibCellsToTriAddrListPtr_;

        mutable PtrList<scalarListList>* ibCellsToTriWeightsListPtr_;

        mutable PtrList< PtrList<scalarSquareMatrix>>* ibCellsToTriMatrixListPtr_;

        mutable PtrList<labelListList>* ghostCellsToTriAddrListPtr_;

        mutable PtrList<scalarListList>* ghostCellsToTriWeightsListPtr_;

        mutable PtrList< PtrList<scalarSquareMatrix>>* ghostCellsToTriMatrixListPtr_;

        //tri-edge-addressing
        mutable PtrList<labelListList>* ibCellsToTriEdgeAddrListPtr_;

        mutable PtrList<scalarListList>* ibCellsToTriEdgeWeightsListPtr_;

        mutable PtrList<labelListList>* ghostCellsToTriEdgeAddrListPtr_;

        mutable PtrList<scalarListList>* ghostCellsToTriEdgeWeightsListPtr_;

        mutable PtrList<labelListList>* proTriFacesInMeshListPtr_;

      
        mutable PtrList<labelListList>* dualMeshAddressingListPtr_;
        


        mutable PtrList<labelList>* oldToNewFcMapListPrt_;

        mutable PtrList<labelList>* oldToNewEdgeMapListPrt_;

        mutable PtrList<labelList>* newToOldEdgeMapListPrt_;
        
        mutable PtrList<labelList>* oldToNewDualEdgeMapListPrt_;
        
        mutable PtrList<labelList>* newToOldDualEdgeMapListPrt_;

        // read from dictionary
        mutable scalar radiusFactor_;
        mutable scalar distFactor_;
        mutable scalar maxCellCellRows_;
        mutable scalar maxFindCellCellNum_;
        mutable scalar maxSortedCellCellNum_;
        mutable scalar maxCellCellNum_;
        mutable word BCtype_;
        mutable scalar gammI_c;


        mutable PtrList<vectorField>* wallShearStressListPtr_;
        mutable PtrList<scalarField>* nutListPtr_;

        // sediment related
        // deposition
        mutable PtrList<scalarField>* depositionListPtr_;

        // entrainment
        mutable PtrList<scalarField>* entrainmentListPtr_;

        // accumulated dIta
        mutable PtrList<scalarField>* accumDItaListPtr_;
       
        //
        mutable PtrList<fvMesh>* dualMeshListPtr_;
        
        mutable PtrList<bool>* dualListPtr_;

    // Private Member Functions
        //- if need to update surface mesh information
        bool updateIB() const
        {
            return updateIB_;
        }
        void IBNeedUpdated() const
        {
            updateIB_ = true;
        }
        void IBHasUpdated() const
        {
            updateIB_ = false;
        }

        //- read dictionary, surface mesh
        void readDict();
        
        void readDualDict();
        
        bool readStl();  
        //- mark IB/Ghost cells
        void markCells() const;

        //- check if all the PointsInFluid are inside mesh, make cut cells
        const cellClassification cellType(const label& objectID) const;

        //- make gammaCellTypeListPtr_
        void makeGammaCellType
        (
            const cellClassification& cellType,
            const label& objectID
        ) const;

        //- check if cellID center in fluid domain
        bool ifInFluid
        (
            const label& cellID,
            const triSurfaceSearch& querySurf,
            const label& objectID
        ) const;

        //- make ghost cells and IB cells, as well as fluid cells and solid cells
        void makeGhostAndIbCells(const label& objectID) const;

        //Delete repete ibCells and ghostCells
        void deleteExtraIbCellsAndGhostCells() const;

        //- make ibGamma, ibFace information
        void makeIbInfo(const label& objectID) const;


        //- make face gamma information
        void makeSGamma(const label& objectID) const;

        //- make gammaCellTypePtr_
        void makeGamma
        (
            const List<cellInfo>& cellInfoList
        ) const;

        //- make sGammaCellTypePtr_

        template<class Type>
        void evaluateCoupled
        (
            GeometricField<Type, fvPatchField, volMesh>& volValues
        )const;

        template<class Type>
        void evaluateUnCoupled
        (
            GeometricField<Type, fvPatchField, volMesh>& volValues
        )const;
        scalar cellSize(label cellID) const;
        scalar cellLength(label cellID) const;

        //-
        void makeStencilsInfo() const;

        void makeIbHitAndSamplingPoints(const label& objectID) const;

        void makeGhostHitAndImagePoints(const label& objectID) const;

        void makeIbFacesHitAndSamplingPoints(const label& objectID) const;

        void makeSamplingStencils(const label& objectID) const;
        void makeFarSamplingStencils(const label& objectID) const;
        void makeIbCellStencils(const label& objectID) const;

        void makeIbFaceSamplingStencils(const label& objectID) const;

        void makeImageStencils(const label& objectID) const;

        void correctIbCell
        (
            const label& objectID,
            labelHashSet& ibCellsSet,
            labelHashSet& ghostCellsSet,
            volScalarField& gammaExt
        ) const;

        labelList findCellID(const pointField& points) const;

        void makeSamplingStencilsWeights(const label& objectID) const;

        void makeIbCellStencilsWeights(const label& objectID) const;

        void makeIbFaceSamplingStencilsWeights(const label& objectID) const;

        void makeImageStencilsWeights(const label& objectID) const;

        void makeInvDirichletMatrices(const label& objectID) const;
        void makeInvNeumannMatrices(const label& objectID) const;

     
        //- send and receive parallel
        template<class Type>
        Foam::FieldField<Foam::Field, Type>
        sendAndReceive
        (
            const Field<Type>& psi,
            const labelListList& procCells,
            const vectorListList& CProc
        ) const;

        // reconstruction
        template<class Type>
        void fromSPointReconstruction
        (
            GeometricField<Type,fvPatchField,volMesh>& psi,
            label objectID
        ) const;

        void fromIbFaceSPointReconstruction
        (
            surfaceScalarField& phi,
            const volVectorField& U,
            label objectID
        ) const;

        template<class Type>
        void fromSPointReconstructionNeumann
        (
            GeometricField<Type,fvPatchField,volMesh>& psi,
            label objectID
        ) const;



        template<class Type>
        void ghostNeumannReconstruction
        (
            GeometricField<Type,fvPatchField,volMesh>& psi,
            label objectID
            ) const;

        template<class Type>
        void correctDiag
        (
            fvMatrix<Type>& eqn,
            label objectID
        ) const;
        
       
        void manipulateMatrix
        (
            fvScalarMatrix& eqn,
            label objectID
        ) const;
        
        void manipulateMatrix
        (
            fvVectorMatrix& eqn,
            label objectID
        ) const;
        
        
        template<class Type>
        void ghostCorrectDiag
        (
            fvMatrix<Type>& eqn,
            label objectID
        ) const;

        void correctOffDiag
        (
            fvScalarMatrix& eqn,
            label objectID
        ) const;



        // - turbulence model
        void kEpsilonCorrection
        (
            const turbulenceModel& turbulence,
            label objectID
        ) const;

        void kOmegaIbCorrection
        (
            const turbulenceModel& turbulence,
            label objectID
        ) const;
        void nutIbCorrection
        (
            label objectID
        ) const;



        void SmagorinskyCorrection
        (
            const turbulenceModel& turbulence,
            const volScalarField& k,
            label objectID
        ) const;

        void SmagorinskyGhostCorrection
        (
            const turbulenceModel& turbulence,
            const volScalarField& k,
            label objectID
        ) const;


        // - mapping
        //- main program
        void makeTriAddressing() const;

        void makeTriAddressing(const label& objectID, bool ifGhost) const;

        void makeTriEdgeAddressing(const label& objectID, bool ifGhost) const;


        
        // - extrude mesh
        void makeDualMesh() const;
        void makeDualMesh(const label& objectID) const;
    
   
        
        template<class Type>
        void setDualEdgeValueToDualMesh
        (
            const label& objectID,
            GeometricField<Type, fvsPatchField, surfaceMesh>& volValues,
            const Field<Type>& values,
            const bool setBC = false
        )const;
             
        template<class Type>
        Field<Type>  mapFromDualPatchToTriPoints
        (
            const Field<Type>& dualValues,
            label objectID
        )const;
        
        template<class Type>
        Field<Type>  mapFromTriPointsToDualPatch
        (
            const Field<Type>& triPoints,
            label objectID
        )const;        
        
    
        template<class Type>
        void dualMeshToPatch
        (
            const label& objectID,
            const GeometricField<Type, fvPatchField, volMesh>& volValues,
            Field<Type>& values,
            const bool setBC = false
        )const;        

        template<class Type>
        void setValueToDualMesh
        (
            const label& objectID,
            GeometricField<Type, fvPatchField, volMesh>& volValues,
            const Field<Type>& values,
            const bool setBC = false
        )const;

        template<class T>
        void sendToMaster
        (
            List<T>& Values
        ) const;

        template<class Type>
        void hitPointExportToMesh
        (
            const Field<Type>& psi,
            GeometricField<Type, fvPatchField, volMesh>& PSI,
            label objectID,
            bool write = true,
            bool dual = false
        )const;

        template<class Type>
        void hitPointExportToMesh
        (
            const word& psiName,
            const Field<Type>& psi,
            label objectID,
            bool write = true
        )const;

        
        template<class Type>
        void hitPointExportToDualMesh
        (
            const Field<Type>& psi,
            GeometricField<Type, fvsPatchField, surfaceMesh>& PSI,
            label objectID,
            bool write = true
        )const;

        template<class Type>
        void hitPointExportToDualMesh
        (
            const word& psiName,
            const Field<Type>& psi,
            label objectID,
            bool write = true
        )const;        

        template<class Type>
        void fieldExportToMesh
        (
            const word& psiName,
            const Field<Type>& psi,
            label objectID
        )const;

        // pressure integral
        void pressureCoeff
        (
            label objectID
        )const;

        // used in MPI to write extrudedMesh results to serial

        template<class Type, template<class> class PatchField, class GeoMesh>
        void parallelWrite
        (
            const GeometricField<Type, PatchField, GeoMesh>& volValues
        ) const;

        template<class Type>
        void parallelWrite
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& surfValues
        ) const;

        template<class Type>
        void parallelWriteMesh
        (
            const IOList<Type>& IOField
        ) const;
        template<class Type>
        void parallelWriteMesh
        (
            const IOField<Type>& IOField
        ) const;

        void parallelWriteMesh
        (
            const faceCompactIOList& IOField
        ) const;

        void parallelWriteMesh
        (
            const polyBoundaryMesh& IOField
        ) const;
        template<class ZoneType>
        void parallelWriteMesh
        (
            const ZoneMesh< ZoneType, polyMesh >& IOField
        )const;


        void smoothOutsidePoints
        (
            pointField& newPoints,
            const label& objectID
        )const ;

        void makeFunctionObjectList
        (
            const label& objectID
        ) const;


protected:


    // Protected Member Functions

    //- Disallow default bitwise copy construct
    immersedBoundaryFvMesh(const immersedBoundaryFvMesh&);

    //- Disallow default bitwise assignment
    void operator=(const immersedBoundaryFvMesh&);

    // Protected Member Functions



public:

    //- Runtime type information
    TypeName("immersedBoundaryFvMesh");


    // Constructors

    //- Construct from IOobject
    immersedBoundaryFvMesh(const IOobject& io);


    //- Destructor
    ~immersedBoundaryFvMesh();


    // Member Functions



    //- Update the mesh for initialization
    bool update();

    //- clear out all pointers
    void clearOut();
    void finalClearOut();

    //- return fvMesh
    const fvMesh& mesh() const
    {
        return *this;
    }

    const fvMesh& globalMesh() const
    {
        return *globalMeshPtr_;
    }

    //- Return pointer information

    //- return mesh search
    const PtrList<word>& IBtypeList() const
    {
        return *IBtypeListPtr_;
    }

    //- return mesh search
    const meshSearch& queryMesh() const
    {
        return *queryMeshPtr_;
    }

    //- return main dictionary
    const IOdictionary& ibProperties() const
    {
        return *ibPropertiesPtr_;
    }

    //- return object sub-dict list
    const PtrList<dictionary>& objectDictList() const
    {
        return *objectDictListPtr_;
    }

    //- return surface mesh list for objects
    const PtrList<triSurfaceMesh>& triSurfaceMeshList() const
    {
        return *triSurfaceMeshListPtr_;
    }

    //- return surface mesh triSurfaceSearch list for objects
    PtrList<triSurface>& objectsList() const
    {
        return *objectsListPtr_;
    }

    PtrList<triSurface>& addObjectsList() const
    {
        return *addObjectsListPtr_;
    }
            
            
    const PtrList<word>& addObjectsNameList() const
    {
        return *addObjectsNameListPtr_;
    }
    const PtrList<label>& addObjectsIDList() const
    {
        return *addObjectsIDListPtr_;
    }
    //- return
    const PtrList<sampledTriSurfaceMesh>& sampledTriSurfaceMeshList() const
    {
        return *sampledTriSurfaceMeshListPtr_;
    }

    const PtrList<cutTriSurfaceMesh>& cutTriSurfaceMeshList() const
    {
        return *cutTriSurfaceMeshListPtr_;
    }

    //- return surface mesh triSurfaceSearch list for objects
     PtrList<triSurfaceSearch>& triSurfaceSearchList() const
    {
        return *triSurfaceSearchListPtr_;
    }

    //- indicator for fluid cells (1), cut cells (0), and solid cells(-1)
    const PtrList<volScalarField>& gammaCellTypeList() const
    {
        return *gammaCellTypeListPtr_;
    }

    const volScalarField& gammaCellType() const
    {
        return *gammaCellTypePtr_;
    }

       //- global gamma live cells (1), IB cells (0.5), ghost cells (-0.5), and dead cells(-1)
    const volScalarField& gamma() const
    {
        return *gammaPtr_;
    }
    const labelList& globalLiveCells() const
    {
        return *globalLiveCellsPtr_;
    }
    const labelList& globalIbCells() const
    {
        return *globalIbCellsPtr_;
    }
    const labelList& globalGhostCells() const
    {
        return *globalGhostCellsPtr_;
    }
    const labelList& globalDeadCells() const
    {
        return *globalDeadCellsPtr_;
    }

    //- IB cell stencil information
    const PtrList<labelList>& ibCellsList() const
    {
        return *ibCellsListPtr_;
    }
    const PtrList<labelList>& oldIbCellsList() const
    {
        return *oldIbCellsListPtr_;
    }

    const PtrList<labelList>& ibDeadCellsList() const
    {
        return *ibDeadCellsListPtr_;
    }
    
    const PtrList<labelList>& oldIbDeadCellsList() const
    {
        return *oldIbDeadCellsListPtr_;
    }
    
    const PtrList<labelList>& ibLiveCellsList() const
    {
        return *ibLiveCellsListPtr_;
    }
    
    const PtrList<labelList>& oldIbLiveCellsList() const
    {
        return *oldIbLiveCellsListPtr_;
    }
    
    const PtrList<volScalarField>& ibGammaList() const
    {
        return *ibGammaListPtr_;
    }

    const PtrList<surfaceScalarField>& sGammaList() const
    {
        return *sGammaListPtr_;
    }

    const PtrList<labelList>& ibFacesList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!ibFacesListPtr_->set(objectID))
            {
                makeIbFacesHitAndSamplingPoints(objectID);
            }
        }
        return *ibFacesListPtr_;
    }

    const PtrList<labelList>& ibFaceCellsList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!ibFaceCellsListPtr_->set(objectID))
            {
                makeIbFacesHitAndSamplingPoints(objectID);
            }
        }
        return *ibFaceCellsListPtr_;
    }

    const PtrList<List<bool>>& ibFaceFlipsList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!ibFaceFlipsListPtr_->set(objectID))
            {
                makeIbFacesHitAndSamplingPoints(objectID);
            }
        }
        return *ibFaceFlipsListPtr_;
    }

    const PtrList<pointField>& ibFaceHitPointsList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!ibFaceHitPointsListPtr_->set(objectID))
            {
                makeIbFacesHitAndSamplingPoints(objectID);
            }
        }
        return *ibFaceHitPointsListPtr_;
    }

    const PtrList<labelList>& ibFaceHitFacesList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!ibFaceHitFacesListPtr_->set(objectID))
            {
                makeIbFacesHitAndSamplingPoints(objectID);
            }
        }
        return *ibFaceHitFacesListPtr_;
    }

    const PtrList<pointField>& ibFaceSamplingPointsList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!ibFaceSamplingPointsListPtr_->set(objectID))
            {
                makeIbFacesHitAndSamplingPoints(objectID);
            }
        }
        return *ibFaceSamplingPointsListPtr_;
    }

    const PtrList<pointField>& ibHitPointsList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!ibHitPointsListPtr_->set(objectID))
            {
                makeIbHitAndSamplingPoints(objectID);
            }
        }
        return *ibHitPointsListPtr_;
    }

    const PtrList<labelList>& ibHitFacesList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!ibHitFacesListPtr_->set(objectID))
            {
                makeIbHitAndSamplingPoints(objectID);
            }
        }
        return *ibHitFacesListPtr_;
    }

    const PtrList<bool>& ibTriNormalsFlipList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!ibTriNormalsFlipListPtr_->set(objectID))
            {
                makeIbHitAndSamplingPoints(objectID);
            }
        }
        return *ibTriNormalsFlipListPtr_;
    }

    const PtrList<pointField>& samplingPointsList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!samplingPointsListPtr_->set(objectID))
            {
                makeIbHitAndSamplingPoints(objectID);
            }
        }
        return *samplingPointsListPtr_;
    }
    const PtrList<pointField>& farSamplingPointsList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!farSamplingPointsListPtr_->set(objectID))
            {
                makeIbHitAndSamplingPoints(objectID);
            }
        }
        return *farSamplingPointsListPtr_;
    }


    //- ghost cell stencil information
    const PtrList<labelList>& ghostCellsList() const
    {
        return *ghostCellsListPtr_;
    }

    const PtrList<pointField>& ghostHitPointsList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!ghostHitPointsListPtr_->set(objectID))
            {
                makeGhostHitAndImagePoints(objectID);
            }
        }
        return *ghostHitPointsListPtr_;
    }

    const PtrList<labelList>& ghostHitFacesList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!ghostHitFacesListPtr_->set(objectID))
            {
                makeGhostHitAndImagePoints(objectID);
            }
        }
        return *ghostHitFacesListPtr_;
    }

    const PtrList<pointField>& imagePointsList() const
    {
        forAll(objectDictList(),objectID)
        {
            if(!imagePointsListPtr_->set(objectID))
            {
                makeGhostHitAndImagePoints(objectID);
            }
        }
        return *imagePointsListPtr_;
    }
    const immersedBoundaryStencils& samplingStencils(const label& objectID) const;
    const immersedBoundaryStencils& farSamplingStencils(const label& objectID) const;
    const immersedBoundaryStencils& ibCellStencils(const label& objectID) const;

    const immersedBoundaryStencils& ibFaceSamplingStencils(const label& objectID) const;

    const immersedBoundaryStencils& imageStencils(const label& objectID) const;


    //- names of objects
    const word objectNames(const label& objectID) const
    {
        return ibProperties().subDict("objects").toc()[objectID];
    }

    const PtrList<labelListList>& ibCellsToTriAddrList() const
    {
        return *ibCellsToTriAddrListPtr_;
    }

    const PtrList<scalarListList>& ibCellsToTriWeightsList() const
    {
        return *ibCellsToTriWeightsListPtr_;
    }

    const PtrList<PtrList<scalarSquareMatrix>>& ibCellsToTriMatrixList() const
    {
        return *ibCellsToTriMatrixListPtr_;
    }

    const PtrList<labelListList>& ghostCellsToTriAddrList() const
    {
        return *ghostCellsToTriAddrListPtr_;
    }

    const PtrList<scalarListList>& ghostCellsToTriWeightsList() const
    {
        return *ghostCellsToTriWeightsListPtr_;
    }

    const PtrList<PtrList<scalarSquareMatrix>>& ghostCellsToTriMatrixList() const
    {
        return *ghostCellsToTriMatrixListPtr_;
    }

    const PtrList<labelListList>& ibCellsToTriEdgeAddrList() const
    {
        return *ibCellsToTriEdgeAddrListPtr_;
    }

    const PtrList<scalarListList>& ibCellsToTriEdgeWeightsList() const
    {
        return *ibCellsToTriEdgeWeightsListPtr_;
    }

    const PtrList<labelListList>& ghostCellsToTriEdgeAddrList() const
    {
        return *ghostCellsToTriEdgeAddrListPtr_;
    }

    const PtrList<scalarListList>& ghostCellsToTriEdgeWeightsList() const
    {
        return *ghostCellsToTriEdgeWeightsListPtr_;
    }
    
    const PtrList<labelListList>& proTriFacesInMeshList() const
    {
        return *proTriFacesInMeshListPtr_;
    }



    
    const labelListList& dualMeshAddressing(const label& objectID) const
    {
        if(!dualMeshAddressingListPtr_)
        {
            dualMeshAddressingListPtr_= new PtrList<labelListList>(objectDictList().size());
        }
        if(!dualMeshAddressingListPtr_->set(objectID))
        {
            makeDualMeshAddressing(objectID);
        }

        PtrList<labelListList>& dualMeshAddressingList = *dualMeshAddressingListPtr_;

        return dualMeshAddressingList[objectID];
    }
    void makeDualMeshAddressing(label objectID) const;

   
       
    const labelList& oldToNewFcMap(const label& objectID) const
    {
        PtrList<labelList>& oldToNewFcMapList = *oldToNewFcMapListPrt_;

        return oldToNewFcMapList[objectID];
    }

    const labelList& oldToNewEdgeMap(const label& objectID) const
    {
        PtrList<labelList>& oldToNewEdgeMapList = *oldToNewEdgeMapListPrt_;

        return oldToNewEdgeMapList[objectID];
    }

    const labelList& newToOldEdgeMap(const label& objectID) const
    {
        PtrList<labelList>& newToOldEdgeMapList = *newToOldEdgeMapListPrt_;

        return newToOldEdgeMapList[objectID];
    }
    
    const labelList& oldToNewDualEdgeMap(const label& objectID) const
    {
        PtrList<labelList>& oldToNewDualEdgeMapList = *oldToNewDualEdgeMapListPrt_;

        return oldToNewDualEdgeMapList[objectID];
    }            
    
    const labelList& newToOldDualEdgeMap(const label& objectID) const
    {
        PtrList<labelList>& newToOldDualEdgeMapList = *newToOldDualEdgeMapListPrt_;

        return newToOldDualEdgeMapList[objectID];
    }            

    //- write IB-cell field
    void writeibCell() const;

    //- write ghost-cell field
    void writeGhostCell() const;

    //- update cell values, called in solvers
    //-update U before momentum predictor
    void evaluateU() const;
    // set initial condition
    void initialIB() const;

    // set Inlet boundary condition
    void setInlet() const;

    //-update p
    void evaluateP() const;
    
    //-update k
    void evaluatek() const;
    
    //-update C before C predictor
    void evaluateC() const;

    //-update UEqn
    void manipulateUEqn
    (
        fvVectorMatrix& UEqn
    ) const;

    void ibCellForcing
    (
        fvVectorMatrix& UEqn,
        label objectID
    ) const;


    void ibCellReconstruction
    (
        fvVectorMatrix& UEqn,
        label objectID
    ) const;


    void ghostCellReconstruction
    (
        fvVectorMatrix& UEqn,
        label objectID
    ) const;

    template<class Type>
    Foam::Field<Type> stencilInterpolation
    (
        const Field<Type>& psi,
        label objectID,
        const labelList& ibc,
        const labelListList& ibcc,
        const List<List<labelPair> >& ibcProcC,
        const labelListList& procCells,
        const vectorListList& CProc,
        const scalarListList& cellWeights,
        const scalarListList& cellProcWeights
    )const;

    template<class Type>
    Foam::Field<Type> samplingPointsValues
    (
        const Field<Type>& psi,
        label objectID
    )const;
    
    template<class Type>
    Foam::Field<Type> farSamplingPointsValues
    (
        const Field<Type>& psi,
        label objectID
    )const;
    
    template<class Type>
    Foam::Field<Type> ibFaceSamplingPointsValues
    (
        const Field<Type>& psi,
        label objectID
    )const;

    template<class Type>
    Foam::Field<Type> ibCellCentersValues
    (
        const Field<Type>& psi,
        label objectID
    )const;

    template<class Type>
    Foam::Field<Type> imagePointsValues
    (
        const Field<Type>& psi,
        label objectID
    )const;

    void correctDiag
    (
        fvVectorMatrix& eqn
    ) const;
    
    void correctDiag
    (
        fvScalarMatrix& eqn
    ) const;
                
    void manipulateMatrix
    (
        fvVectorMatrix& eqn
    ) const;
    void manipulateMatrix
    (
        fvScalarMatrix& eqn
    ) const;


    void correctOffDiag
    (
        fvScalarMatrix& eqn
    ) const;

    

    //- Adjust the fluxes on immersed boundary to obey continuity.

    void ibCorrectPhi
    (
        surfaceScalarField& phi
    ) const;

    void immersedBoundaryAdjustPhi
    (
        surfaceScalarField& phi,
        const volVectorField& U
    ) const;

    void immersedBoundaryAdjustPhi
    (
        surfaceScalarField& phi,
        const volVectorField& U,
        const label& objectID
    ) const;

    void ibFacesAdjust
    (
        surfaceScalarField& phi,
        const volVectorField& U
    ) const;

    void ibFacesAdjust
    (
        surfaceScalarField& phi,
        const volVectorField& U,
        const label& objectID
    ) const;

    void adjustDivPhiHbyA
    (
        volScalarField& divPhiHbyA
    ) const;

    void adjustDivPhiHbyA
    (
        volScalarField& divPhiHbyA,
        const label& objectID
    ) const;

    void evaluateGradP
    (
        volVectorField& psi
    ) const;

// - Turbulence model
    //- k-Epsilon
    void kEpsilonCorrection(const turbulenceModel& turbulence) const;

    //- k-Omega
    void kOmegaCorrection(const turbulenceModel& turbulence) const;
    
    void nutCorrection() const;
    
    // Smagorinsky
    void SmagorinskyCorrection
    (
        const turbulenceModel& turbulence,
        const volScalarField& k
    ) const;
// - mapping
    template<class Type>
    Foam::Field<Type> mapFromIBHitToTriFace
    (
        const Field<Type>& ibValues,
        label objectID
    )const;

    template<class Type>
    Foam::Field<Type> mapFromGHOSTHitToTriFace
    (
        const Field<Type>& ghostValues,
        label objectID
    )const;

    template<class Type>
    Foam::Field<Type> mapFromIBHitToTriEdge
    (
        const Field<Type>& ibValues,
        label objectID
    )const;

    template<class Type>
    Foam::Field<Type> mapFromGHOSTHitToTriEdge
    (
        const Field<Type>& ghostValues,
        label objectID
    )const;
    
    template<class Type>
    Foam::Field<Type> mapFromIBHitToDualEdge
    (
        const Field<Type>& ibValues,
        label objectID
    )const;

    template<class Type>
    Foam::Field<Type> mapFromGHOSTHitToDualEdge
    (
        const Field<Type>& ghostValues,
        label objectID
    )const;            

    template<class Type>
    Foam::Field<Type> mapFromIBHitToDualFace
    (
        const Field<Type>& ibValues,
        label objectID
    )const;

    template<class Type>
    Foam::Field<Type> mapFromGHOSTHitToDualFace
    (
        const Field<Type>& ghostValues,
        label objectID
    )const; 
    
    template<class Type>
    void solveMatrix
    (
        Type& psi,
        const scalarSquareMatrix& matrix,
        const Field<Type>& source,
        const point& point,
        const scalarList& weight

    )const;

    // - vtk writer
    template<class Type>
    void writeVTK
    (
        const fileName &surfaceName,
        const word &fieldName,
        const Field<Type>& psi,
        label objectID
    ) const;

    // - write out calculated values
    const vectorField& wallShearStress(const label& objectID)const
    {
        const PtrList<vectorField>& tauWall= *wallShearStressListPtr_;
        return tauWall[objectID];
    }
    
    const scalarField& nutNew(const label& objectID)const
    {
        const PtrList<scalarField>& nutNew= *nutListPtr_;
        return nutNew[objectID];
    }

    // sediment related
    // deposition
    const scalarField& deposition(const label& objectID)const
    {
        const PtrList<scalarField>& depositionList= *depositionListPtr_;
        return depositionList[objectID];
    }

    // entrainment
    const scalarField& entrainment(const label& objectID)const
    {
        const PtrList<scalarField>& entrainmentList= *entrainmentListPtr_;
        return entrainmentList[objectID];
    }

    // accumulated dIta
    const scalarField& accumDItaList(const label& objectID)const
    {
        const PtrList<scalarField>& accumDItaList= *accumDItaListPtr_;
        return accumDItaList[objectID];
    }


    // -
    const PtrList<fvMesh>& dualMeshList()const
    {

        return *dualMeshListPtr_;
    }
    
    //- if using dual mesh 
    const PtrList<bool>& dualList()const
    {
        return *dualListPtr_;
    }
    
    const bool& dual(const label& objectID)const
    {
        return dualList()[objectID];
    } 
        
             
    
    // - Post evalulation

    void postEvaulation() const;


    
    //- using Roulund et al 2005, Jacobsen et al 2014, Batkal et al 2015
    void sediment_dual(const label& objectID) const;            

    void updateVs(const label& objectID) const;

    // - all kinds of different calc tools, located in Calcs.C
    void pressureOutput() const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
