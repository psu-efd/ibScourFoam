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
#include "PatchTools.H"
#include "fvcSmooth.H"
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::immersedBoundaryFvMesh::postEvaulation() const
{
    
    const double Oldtime1 = time().elapsedCpuTime();

    forAll(objectsList(),objectID)
    {
        
        if(objectDictList()[objectID].found("sediment"))
        {

            IOdictionary newDict
            (
                IOobject
                (
                    "immersedBoundaryProperties",
                    mesh().time().constant(),
                    mesh(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            );
            bool changeSTL = newDict.subDict("objects")
                            .subDict(objectNames(objectID))
                            .subDict("sediment")
                            .lookupOrDefault<bool>("changeSTL",true);

            if(changeSTL)
            {
                if(!updateIB())
                {
                    Info<<endl<<"****** Change to mobile bed ******"<<endl<<endl;
                }
                IBNeedUpdated();
            }
            else
            {
                if(updateIB())
                {
                    Info<<endl<<"****** Change to fixed bed ******"<<endl<<endl;
                }
                IBHasUpdated();
            }

            // It provides 4 options.
            // Seen in Foam::immersedBoundaryFvMesh::updateVs
            updateVs(objectID);
                      
            sediment_dual(objectID);

        }

    }
    const double Oldtime2=time().elapsedCpuTime();

    if(debug)
    {
        Info<<"postEvaulation Executation Time = "<<Oldtime2-Oldtime1<< " s"<<endl;
    }
}




template<class Type>
void Foam::immersedBoundaryFvMesh::dualMeshToPatch
(
    const label& objectID,
    const GeometricField<Type, fvPatchField, volMesh>& volValues, // dualMesh
    Field<Type>& values, // triFace values
    const bool setBC  // default = false, let boundary faces has same as corresponding boundary edge
)const
{
    const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];
    const labelList& includeFaces = cTSM.includeDualFaces();

    // fill includeFaces
    const Field<Type>& volValuesI = volValues.primitiveField();

    forAll(includeFaces,I)
    {
        values[includeFaces[I]]=volValuesI[I];
    }

    if(setBC)
    {
        labelHashSet includeFacesSet(includeFaces);
        const labelListList& edgeFaces = cTSM.newDualPatch().newPatch().edgeFaces();
        const labelList& boundaryEdges = cTSM.boundaryDualEdges();
        forAll(boundaryEdges,I)
        {
            label edgeID = boundaryEdges[I];
            const labelList& edgeFace = edgeFaces[edgeID];
            if(edgeFace.size()>1)
            {
                if(includeFacesSet.found(edgeFace[0]))
                {
                    values[edgeFace[1]]=values[edgeFace[0]];
                }
                else
                {
                    values[edgeFace[0]]=values[edgeFace[1]];
                }
            }
        }
    }
}

// also make sure wall patchField type
template<class Type>
void Foam::immersedBoundaryFvMesh::setValueToDualMesh
(
    const label& objectID,
    GeometricField<Type, fvPatchField, volMesh>& volValues, // dualMesh
    const Field<Type>& values, // patch values
    const bool setBC  // default = false   only for flux
)const
{
    const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];
    const labelList& includeFaces = cTSM.includeDualFaces();

    // fill includeFaces
    Field<Type>& volValuesI = volValues.primitiveFieldRef();

    forAll(includeFaces,I)
    {
        volValuesI[I]=values[includeFaces[I]];
    }

    const fvMesh& mesh=volValues.mesh();
    if(setBC)
    {
        forAll(volValues.boundaryFieldRef(),patchI)
        {
            if (isA<wallPolyPatch>(mesh.boundaryMesh()[patchI]))
            {

                volValues.boundaryFieldRef().set
                (
                    patchI,
                    fvPatchField<Type>::New("zeroGradient", mesh.boundary()[patchI], volValues)
                );
                forAll(volValues.boundaryFieldRef()[patchI],I)
                {
                    volValues.boundaryFieldRef()[patchI][I] = pTraits<Type>::zero; //only for flux
                }
            }
            else if (!isA<emptyPolyPatch>(mesh.boundaryMesh()[patchI]))
            {

                volValues.boundaryFieldRef().set
                (
                    patchI,
                    fvPatchField<Type>::New("zeroGradient", mesh.boundary()[patchI], volValues)
                );
            }
        }
    }

    // Evaluate the coupled patchField
    evaluateCoupled(volValues);
    evaluateUnCoupled(volValues);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::immersedBoundaryFvMesh::parallelWrite
(
    const GeometricField<Type, PatchField, GeoMesh>& volValues

) const
{
    if(Pstream::parRun())
    {
        fileName casePath = time().rootPath()/time().globalCaseName()/time().timeName()/volValues.db().dbDir();
        mkDir(casePath);

        fileName filePath =  casePath/volValues.name();

        autoPtr<Ostream> osPtr
        (
            fileHandler().NewOFstream
            (
                filePath,
                time().writeFormat(),
                IOstream::currentVersion,
                time().writeCompression()
            )
        );
        volValues.writeHeader(osPtr());
        volValues.writeData(osPtr());
        IOobject::writeEndDivider(osPtr());
    }
    else
    {
        volValues.write();
    }
}

template<class Type>
void Foam::immersedBoundaryFvMesh::parallelWriteMesh
(
    const IOList<Type>& IOField
) const
{
    fileName casePath = time().rootPath()/time().globalCaseName()/IOField.instance()/IOField.db().dbDir()/IOField.local();
    mkDir(casePath);
    fileName filePath =  casePath/IOField.name();
    autoPtr<Ostream> osPtr
    (
        fileHandler().NewOFstream
        (
            filePath,
            time().writeFormat(),
            IOstream::currentVersion,
            time().writeCompression()
        )
    );
    IOField.writeHeader(osPtr());
    IOField.writeData(osPtr());
    IOobject::writeEndDivider(osPtr());
}
template<class Type>
void Foam::immersedBoundaryFvMesh::parallelWriteMesh
(
    const IOField<Type>& IOField
) const
{
    fileName casePath = time().rootPath()/time().globalCaseName()/IOField.instance()/IOField.db().dbDir()/IOField.local();
    mkDir(casePath);
    fileName filePath =  casePath/IOField.name();
    autoPtr<Ostream> osPtr
    (
        fileHandler().NewOFstream
        (
            filePath,
            time().writeFormat(),
            IOstream::currentVersion,
            time().writeCompression()
        )
    );
    IOField.writeHeader(osPtr());
    IOField.writeData(osPtr());
    IOobject::writeEndDivider(osPtr());
}

void Foam::immersedBoundaryFvMesh::parallelWriteMesh
(
    const faceCompactIOList& IOField
) const
{

    fileName casePath = time().rootPath()/time().globalCaseName()/IOField.instance()/IOField.db().dbDir()/IOField.local();
    mkDir(casePath);
    fileName filePath =  casePath/IOField.name();
    autoPtr<Ostream> osPtr
    (
        fileHandler().NewOFstream
        (
            filePath,
            time().writeFormat(),
            IOstream::currentVersion,
            time().writeCompression()
        )
    );
    IOField.writeHeader(osPtr(),"faceList");
    IOField.writeData(osPtr());
    IOobject::writeEndDivider(osPtr());
}

void Foam::immersedBoundaryFvMesh::parallelWriteMesh
(
    const polyBoundaryMesh& IOField
) const
{
    fileName casePath = time().rootPath()/time().globalCaseName()/IOField.instance()/IOField.db().dbDir()/IOField.local();
    mkDir(casePath);
    fileName filePath =  casePath/IOField.name();
    autoPtr<Ostream> osPtr
    (
        fileHandler().NewOFstream
        (
            filePath,
            time().writeFormat(),
            IOstream::currentVersion,
            time().writeCompression()
        )
    );
    IOField.writeHeader(osPtr());
    IOField.writeData(osPtr());
    IOobject::writeEndDivider(osPtr());
}

template<class ZoneType>
void Foam::immersedBoundaryFvMesh::parallelWriteMesh
(
    const ZoneMesh< ZoneType, polyMesh >& IOField
) const
{
    fileName casePath = time().rootPath()/time().globalCaseName()/IOField.instance()/IOField.db().dbDir()/IOField.local();
    mkDir(casePath);
    fileName filePath =  casePath/IOField.name();
    autoPtr<Ostream> osPtr
    (
        fileHandler().NewOFstream
        (
            filePath,
            time().writeFormat(),
            IOstream::currentVersion,
            time().writeCompression()
        )
    );
    IOField.writeHeader(osPtr());
    IOField.writeData(osPtr());
    IOobject::writeEndDivider(osPtr());
}


template<class Type>
void Foam::immersedBoundaryFvMesh::parallelWrite
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& surfValues
) const
{
    GeometricField<Type, fvPatchField, volMesh> volValues = fvc::average(surfValues);
    volValues.rename(surfValues.name());
    parallelWrite(volValues);
}


void Foam::immersedBoundaryFvMesh::smoothOutsidePoints
(
    pointField& newPoints, // global point ID
    const label& objectID
)const
{
    // read information from dictionary
    const dictionary& dict = objectDictList()[objectID].subDict("sediment");
    vector gravity(dict.lookup("gravity"));

    pointField oldPoints = newPoints;

    const triSurface& surf(objectsList()[objectID]);

    const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];
    const labelList& outsidePoints = cTSM.outsidePoints(); // local point ID
    const labelListList& outsidePointsStencil = cTSM.outsidePointsStencil(); // local point ID

    const scalarListList& outsidePointsStencilWeight = cTSM.outsidePointsStencilWeight(); // local point ID
    const labelList& meshPoints = surf.meshPoints();
    forAll(outsidePoints,I)
    {
        label localPointID = outsidePoints[I];
        label pointID = meshPoints[localPointID]; // global ID

        newPoints[pointID].z()= 0;

        forAll(outsidePointsStencil[I],II)
        {
            label localNearestPointID=outsidePointsStencil[I][II];

            label nearestPointID = meshPoints[localNearestPointID]; // global ID
            scalar weight = outsidePointsStencilWeight[I][II];

            newPoints[pointID].z()=newPoints[pointID].z()+weight*newPoints[nearestPointID].z();
            
        }

    }
}










void Foam::immersedBoundaryFvMesh::pressureOutput() const
{
    forAll(objectsList(),objectID)
    {
        pressureCoeff(objectID);
    }
}


void Foam::immersedBoundaryFvMesh::pressureCoeff(label objectID )const
{
    volScalarField& p = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("p"));

    scalarField triPressure;
    if(IBtypeList()[objectID]=="classic")
    {
        const scalarField ibPressure(p,ibCellsList()[objectID]);
        triPressure = mapFromIBHitToTriFace(ibPressure,objectID);
    }
    else if(IBtypeList()[objectID]=="mix")
    {
        const scalarField ghostPressure(p,ghostCellsList()[objectID]);
        triPressure = mapFromGHOSTHitToTriFace(ghostPressure,objectID);
    }
    const triSurface& surf = objectsList()[objectID];
    vector pForce(pTraits<vector>::zero);

    const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];
    labelHashSet includeFacesSet(cTSM.includeFaces());
    labelHashSet boundaryFacesSet(cTSM.boundaryFaces());
    forAll(triPressure,I)
    {
        if(includeFacesSet.found(I) or boundaryFacesSet.found(I))
        {
            pForce +=triPressure[I]*surf.localFaces()[I].normal(surf.localPoints());
        }
    }
    if(!ibTriNormalsFlipList()[objectID])
    {
        pForce*=-1.0;
    }

    reduce(pForce,sumOp<vector>());

    Info<<"Total force of "<<objectNames(objectID)<<": "<<pForce<<endl;

    if (Pstream::master())
    {
        fileName forcesDir;
        word startTimeName =
         time().timeName(time().startTime().value());

        if (Pstream::parRun())
        {
         // Put in undecomposed case (Note: gives problems for
         // distributed data running)
         forcesDir =
            time().path()/".."/"pressure"/startTimeName;
        }
        else
        {
         forcesDir =
            time().path()/"pressure"/startTimeName;
        }

        // Create directory if does not exist.
        mkDir(forcesDir);

        // Open new file at start up
        OFstream fileWriter(forcesDir/("obj_" + objectNames(objectID) + ".dat"),IOstream::ASCII,IOstream::currentVersion,IOstream::UNCOMPRESSED,true);// append true

        fileWriter<<time().value()<<tab<<pForce[0]<<tab<<pForce[1]<<tab<<pForce[2]<<endl;
    }//end of master
}    
    
 
// ************************************************************************* //

