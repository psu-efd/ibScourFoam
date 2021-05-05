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
    immersedBoundaryFvMeshDualCalcs.C
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvMesh.H"
#include "MeshWave.H"
#include "SortableList.H"
#include "fvCFD.H"
#include "turbulenceModel.H"
#include "fvcSmooth.H"
#include "sedModules.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Main program for scour modelling
// bedLoadFlux is calculated on centers of extrudeMesh / vertices of dualMesh
void Foam::immersedBoundaryFvMesh::sediment_dual(const label& objectID)const
{

    const dictionary& dict = objectDictList()[objectID].subDict("sediment");
    bool changeSTL(dict.lookupOrDefault<bool>("changeSTL",true));
    PtrList<fvMesh>& dualMeshList= *dualMeshListPtr_;
    fvMesh& mesh = dualMeshList[objectID];
    
    
    
    if(changeSTL)
    {
        vector gravity(dict.lookup("gravity"));
        scalar shieldsNumberC0(dict.lookupOrDefault<scalar>("shieldsNumberC0",-1));
        scalar mus(dict.lookupOrDefault<scalar>("staticFrictionCoefficient",0.63));
        scalar mud(dict.lookupOrDefault<scalar>("dynamicFrictionCoefficient",0.51));
        scalar scale_(dict.lookupOrDefault<scalar>("deformScale",1));
        scalar d50(readScalar(dict.lookup("grainDiameter")));
        scalar ssG(dict.lookupOrDefault<scalar>("specificGravity",2.69));
        bool sandSlide_(dict.lookupOrDefault<bool>("sandSlide",true));
        
        scalar clipDIta(dict.lookupOrDefault<scalar>("clipDIta",0.1));

        bool suspendedLoad_ = dict.lookupOrDefault<bool>("suspendedLoad",false);

        const dictionary& transportProperties =
              this->lookupObject<dictionary>
              (
                 "transportProperties"
              );
        dimensionedScalar nu(transportProperties.lookup("nu"));

        PtrList<triSurface>& addSurfs= addObjectsList();
        
        const PtrList<word>& addObjectsName = addObjectsNameList();
        const PtrList<label>& addObjectsID = addObjectsIDList();
        
       
        
            
            
        
        #include "createDualFieldsInSediment.H"
        
         
     
        const vectorField& ibWallShearStress = this->wallShearStress(objectID); // ib value or ghost value

        vectorField wallShearStress; // edge value not flux
     
       
        // map hit point values onto face edges
        if(IBtypeList()[objectID]=="classic")
        {
            wallShearStress = mapFromIBHitToDualEdge(ibWallShearStress,objectID);
        }
        else if(IBtypeList()[objectID]=="mix")
        {
            wallShearStress = mapFromGHOSTHitToDualEdge(ibWallShearStress,objectID);
        }
          
        // wss is surfaceVectorField, 
        // but ita(should be eta) is still volVectorField

        setDualEdgeValueToDualMesh
        (
            objectID,
            wss,// dualMesh edge
            wallShearStress // dual patch Edge value
        );
        

        // sediment transport bed-load

        // read information from dictionary
        

        // evaluate shieldsNumberC0 using modified shields diagram Parker 2005
        scalar Rep=sqrt((ssG-1)*mag(gravity)*d50)*d50/nu.value();
        if(shieldsNumberC0<0)
        {
            shieldsNumberC0=0.22*pow(Rep,-0.6)+0.06*exp(-17.77*pow(Rep,-0.6)); // Brownlie 1981
            if(Rep<4.22)
            {
                shieldsNumberC0=min(shieldsNumberC0,0.135*pow(Rep,-0.261)); // Mantz 1977
            }
            shieldsNumberC0*=0.5;
        }
        
        surfaceScalarField wssMag = mag(wss);

        dimensionedScalar gsd("gsd",sqr(dimLength)/sqr(dimTime),mag(gravity)*(ssG-1)*d50);
        sN = wssMag/gsd;

        // information of the surface
        const PPatch& surf = cutTriSurfaceMeshList()[objectID].newDualPatch().newPatch();
        triSurface& triSurf = objectsList()[objectID];
          
        

        labelList includeDualFaces;
        if(Pstream::parRun())
        {
            includeDualFaces = cutTriSurfaceMeshList()[objectID].includeDualFaces();
        }
        else
        {
            includeDualFaces = cutTriSurfaceMeshList()[objectID].includeDualFaces();
        }
           
        
        // Calculate areas of include faces, corresponding to cutTriSurfaceMesh/dualMesh
        scalarField triAreasCTSM(includeDualFaces.size(),0.0);
        vector gravityNormal = gravity/mag(gravity);
        forAll(includeDualFaces,I)
        {
            vector normal = surf[includeDualFaces[I]].normal(surf.points());
            triAreasCTSM[I]=mag(normal&gravityNormal);
        }

        scalarField zTriPtCenter(triSurf.localPoints().size()); // dual patch

        forAll(zTriPtCenter,ptI)
        {
            zTriPtCenter[ptI] = -(triSurf.localPoints()[ptI]&gravity)/mag(gravity);
        }

        scalarField zCenter = mapFromTriPointsToDualPatch(zTriPtCenter,objectID);
        

            
        setValueToDualMesh
        (
            objectID,
            ita,// dualMesh, will change
            zCenter,// dual patch value
            true
        );

        pointField edgeNormals(surf.nEdges(), Zero);

        const vectorField& faceNormals = surf.faceNormals();
        
        surfaceVectorField eN = q0;
        eN.rename("edgeNormals");
        
        const PPatchInterpolation& PPI = cutTriSurfaceMeshList()[objectID].newDualPatch().dualPatchInterp();
        edgeNormals = PPI.faceToEdgeInterpolate(faceNormals);

        edgeNormals /= mag(edgeNormals)+VSMALL;
        


        setDualEdgeValueToDualMesh
        (
            objectID,
            eN,// dualMesh edge
            edgeNormals // dual patch Edge value
        ); 

        // sedModules is to do sediment transport related calculation
        // seen in src/sedModules
        // It has multiple functions:
        // 1. bedLoadUpdateValues: update bedLoadFlux, shieldsNumber,
        //      shieldsNumberC, entrainment rate, deposition rate on dualMesh
        // 2. ibEDUpdateValues: update entrainment rate, deposition rate on hit points,
        //      so that no interpolation is required from dualMesh to hit points
        // 3. sandSlide: it has 3 options:
        //      3.1 diffusion-based
        //      3.2 velocity-based (not ready)
        //      3.3 geometric-based
        // 4. fixWallFlux: make sure that sum of sed flux accross wall equals to zero
        // 5. feedInlet: make sure inlet sed flux equals to outlet sed flux, not well implemented
        sedModules sediment("Roulund2005");
 
        sediment.bedLoadUpdateValues
        (
            d50,ssG,mus,mud,eN,wss,shieldsNumberC0,gravity,
            BETA,PHI,sN,
            sNC,q0
            
        );
        


        if(suspendedLoad_)
        {
            volScalarField& C = const_cast<volScalarField&>
            (this->lookupObject<volScalarField>("C"));

            volVectorField& Vs = const_cast<volVectorField&>
            (this->lookupObject<volVectorField>("Vs"));
            // evaluate Vs Fredse and Deigaard 1992 pp 201
            scalar Rp = average(mag(Vs.primitiveField()))*d50/nu.value();
            
            scalar NN = 2.39;
            if(Rp>1 and Rp<500)
            {
                NN = 4.45*Foam::pow(Rp,0.10);
            }
            else
            {
                NN = 4.35*Foam::pow(Rp,0.03);
            }    
            forAll(Vs,cellID)
            {
                Vs[cellID] = Vs[cellID]*Foam::pow(1-C[cellID],NN);
            }
        
            
            scalarField Zb(sN.size(),0.0);

            // update deposition using sampled C and Vs, required by ibEDUpdateValues
            scalarField ibCbSample;
            vectorField ibVsSample;

            scalarField ibZb(ibWallShearStress.size(),SMALL);

            scalarField ibEntrainment(ibWallShearStress.size(),0);
            scalarField ibDeposition(ibWallShearStress.size(),0);

            if(IBtypeList()[objectID]=="classic")
            {
                ibCbSample = samplingPointsValues(C,objectID);
            }
            else if(IBtypeList()[objectID]=="mix")
            {
                ibCbSample = imagePointsValues(C,objectID);
                ibVsSample = imagePointsValues(Vs,objectID);
            }

            // update VsSample and CbSample, required by bedLoadUpdateValues
            hitPointExportToDualMesh(ibVsSample,VsSample,objectID,false);
            hitPointExportToDualMesh(ibCbSample,CbSample,objectID,false);
            vectorField ibNormals = samplingPointsList()[objectID]-ibHitPointsList()[objectID];
            ibNormals /= mag(ibNormals)+SMALL;
            sediment.ibEDUpdateValues
            (
                d50,ssG,mus,mud,ibNormals,ibWallShearStress,ibVsSample,ibCbSample,shieldsNumberC0,gravity,ibZb,
                ibEntrainment,ibDeposition
            );

            scalarField Cbstar=ibEntrainment/(mag(ibVsSample&gravity)+SMALL);
        


            // depositionListPtr_ and entrainmentListPtr_ will be used in CEqn.H
            depositionListPtr_->set
            (
                objectID,
                ibDeposition
            );

            entrainmentListPtr_->set
            (
                objectID,
                ibEntrainment
            );

            
            Info<< "SEDIMENT INFOMATION"<<endl;
            Info<<tab << "                        " <<tab << "min"<<tab<<"max"<<tab<<"average"<< endl;
            Info<<tab << "Ref. Height(mm)  " <<tab <<  min(ibZb*1000) <<tab<< max(ibZb*1000) << tab << average(ibZb*1000) << endl;
            Info<<tab << "Cbstar" <<tab <<  min(Cbstar) << tab << max(Cbstar) << tab << average(Cbstar) << endl;
            Info<<tab << "Vs   " <<tab <<  min(mag(Vs.primitiveFieldRef())) << tab << max(mag(Vs.primitiveFieldRef())) << tab << average(mag(Vs.primitiveFieldRef())) << endl;
            Info<<tab << "D50              = "<<d50*1000<<" mm"<<endl;
          
            
            
        }
        
        surfaceScalarField::Boundary& qFluxPatches = qFlux.boundaryFieldRef();
        const surfaceVectorField::Boundary& q0Patches = q0.boundaryField();
           
        qFlux = q0 & mesh.Sf();
        


        
        // Boundary value
        forAll(qFluxPatches,patchID)
        {
            const vectorField& Sf = qFluxPatches[patchID].patch().Sf();
            forAll(qFluxPatches[patchID],I)
            {
                qFluxPatches[patchID][I]=Sf[I]&q0Patches[patchID][I];
                
    
            }
            
           
        }
        // fix zero flux on each wall
        sediment.fixWallFlux(qFlux);

        bool feedInlet(dict.lookupOrDefault<bool>("feedInlet",false));
        if(feedInlet)
        {
            sediment.feedInlet(qFlux);
        }
        
       
            volScalarField divQ = fvc::div(qFlux);
          
        divQ.rename("divQ");
        
        labelHashSet inletBCFacesSet;
        bool fixInlet(dict.lookupOrDefault<bool>("fixInlet",false));
        if(fixInlet)
        {
            Info<<"fix inlet BC"<<endl;
            // find all faces of inlet  // dual mesh inlet
            forAll(mesh.boundaryMesh(),patchi)
            {
                if(mesh.boundaryMesh()[patchi].name()=="inlet")
                {
                    const labelList& faceCells = mesh.boundaryMesh()[patchi].faceCells();
                    const labelListList& cellCells = mesh.cellCells();

                    forAll(faceCells,I)
                    {
                        divQ.boundaryFieldRef()[patchi][I] = 0.0;
                        label cellID = faceCells[I];
                        if(!inletBCFacesSet.found(cellID)) inletBCFacesSet.insert(cellID);
                        divQ.primitiveFieldRef()[cellID] = 0.0;
                        forAll(cellCells[cellID],II)
                        {
                            if(!inletBCFacesSet.found(cellCells[cellID][II]))
                            {
                                inletBCFacesSet.insert(cellCells[cellID][II]);
                            }
                            divQ.primitiveFieldRef()[cellCells[cellID][II]]=0.0;
                        }
                    }
                }
            }
        }

        labelList inletBCFaces = inletBCFacesSet.toc();

        
        scalar porosity = dict.lookupOrDefault<scalar>("porosity",0.4);


        fvScalarMatrix itaEqn
        (
            fvm::ddt(ita)
            +1.0/(1.0-porosity)*divQ
        );
        
       

        itaEqn.solve();

        volScalarField dIta = ita-ita.oldTime();

        dIta.rename("elevationChange");

        scalar updateInterval = ibProperties().lookupOrDefault<scalar>("updateInterval",1);
        if(updateInterval<1) updateInterval = 1;

        scalarField accumDIta = dIta.primitiveField();

       
        if(!accumDItaListPtr_->set(objectID))
        {
            accumDItaListPtr_->set
            (
                objectID,
                accumDIta
            );
        }
        else
        {
            accumDIta +=accumDItaList(objectID);
            accumDItaListPtr_->set
            (
                objectID,
                accumDIta
            );
        }

  

            
        dIta.primitiveFieldRef() = accumDIta;
        dIta = dIta*scale_;

        ita = ita.oldTime()+dIta;
        
        
        Info<< "SEDIMENT INFOMATION"<<endl;

        Info<<tab << "                        " <<tab << "min"<<tab<<"max"<<tab<<"average"<< endl;
   
        Info << tab << "dIta(m)       " << tab << min(dIta).value() << tab << max(dIta).value() << tab << average(dIta).value() << endl;
        Info << tab << "ita(m)       " << tab << min(ita).value() << tab << max(ita).value() << tab << average(ita).value() << endl;
        Info<<tab << "shieldsNumber    " <<tab <<  min(sN.primitiveField()) << tab << max(sN.primitiveField()) << tab << average(sN.primitiveField()) << endl;
        Info<<tab << "shieldsNumberC   " <<tab <<  min(sNC.primitiveField()) << tab << max(sNC.primitiveField()) << tab << average(sNC.primitiveField()) << endl;

        Info<<tab << "wallShearStressMag" <<tab <<  min(wssMag.primitiveField()) << tab << max(wssMag.primitiveField()) << tab << average(wssMag.primitiveField()) << endl;
        
        Info<<tab << "D50              = "<<d50*1000<<" mm"<<endl;
        
        
        
        //Clipping of dIta in case it goes wild
        if(clipDIta>0)
        {
            dIta.min(clipDIta/scale_);
            dIta.max(-clipDIta/scale_);
            if(debug)
            {
                Info<< "After clip, scaled bed elevation change = "
                    << dIta.weightedAverage(mesh.V()).value()*scale_
                    << "  Min(accumDIta) = " << min(dIta).value()*scale_
                    << "  Max(accumDIta) = " << max(dIta).value()*scale_
                    << endl;
            }
        }



       
        if(sandSlide_ and changeSTL and !(time().timeIndex()%label(updateInterval)))
        {
            Info<<"Start sand-slide algorithm"<<endl;

            scalar err1 = dict.lookupOrDefault<scalar>("sand_slide_err1",1e-4);
            scalar err2 = dict.lookupOrDefault<scalar>("sand_slide_err2",1e-5);

            
            scalar sandSlide_maxItr = dict.lookupOrDefault<scalar>("sandSlide_maxItr",200);
            
            scalar sandSlide_diffusivity = dict.lookupOrDefault<scalar>("sandSlide_diffusivity",1);

            scalar timeStep = dict.lookupOrDefault<scalar>("sandSlide_timeStep",0.001);

            scalar oldMass=0.0;
            
            forAll(ita,I)
            {
                oldMass +=(1.0+ita[I])*triAreasCTSM[I];
            }
            
            

            sediment.sandSlide
            (
                ita, // grain size
                mus, // static friction coefficient
                sandSlide_diffusivity, // diffusivity coefficient
                err1,
                err2,
                timeStep,
                sandSlide_maxItr // default 200
            );
            
            dIta = ita-ita.oldTime();
            scalar mass = 0.0;
            forAll(ita,I)
            {
                mass +=(1.0+ita[I])*triAreasCTSM[I];
            }
            Info<<"Mass change after sandSlide: "<<(oldMass-mass)/(oldMass+SMALL)*100.0<<"%, "<<oldMass-mass<<endl;
        }


        // change to degrees
        BETA.primitiveFieldRef()=radToDeg(BETA);
        // PHI has already been changed to Degree
    //    PHI.primitiveFieldRef()=radToDeg(PHI);

        // map ita to dual patch face values
        scalarField dzCenter(surf.size(),0.0);

        scalarField newZCenter(zCenter);
        
        dualMeshToPatch
        (
            objectID,
            ita,// dualMesh
            newZCenter,// dual patch value, will change
            true // default = false
        );

        dualMeshToPatch
        (
            objectID,
            dIta,// dualMesh
            dzCenter,// triFace value, will change
            true // default = false
        );
	
	

        scalarField dzPoint = mapFromDualPatchToTriPoints(dzCenter,objectID);
        scalarField zPoint = mapFromDualPatchToTriPoints(newZCenter,objectID);
        
       // if(debug)
         {
             Info<< "dzPoint " << gMin(dzPoint) << " " << gMax(dzPoint) << " " << gAverage(dzPoint) << endl;
             Info<< "zPoint " << gMin(zPoint) << " " << gMax(zPoint) << " " << gAverage(zPoint) << endl;
         }

        scalar bottomLimit_ = dict.lookupOrDefault<scalar>("bottomLimit",-10);

        scalar upperLimit_ = dict.lookupOrDefault<scalar>("upperLimit",-10); 
        Info << "upperLimit  "<<upperLimit_<<endl;
       scalar upperLimitValue_ = dict.lookupOrDefault<scalar>("upperLimitValue",0);
       Info << "upperLimitValue  "<<upperLimitValue_<<endl;
       pointField newPoints = triSurf.points();
        forAll(newPoints,pointI)
        {
            {
                newPoints[pointI] = newPoints[pointI]-dzPoint[triSurf.whichPoint(pointI)]*gravity/mag(gravity);
                
                if(newPoints[pointI].z() < bottomLimit_)  newPoints[pointI].z() = bottomLimit_;
                if(newPoints[pointI].x() <upperLimit_)  newPoints[pointI].z() = upperLimitValue_;
                
                
            }
            //points --- global index to local index
        }
        
        if(changeSTL  and !(time().timeIndex()%label(updateInterval)))
        {
            smoothOutsidePoints(newPoints,objectID);
            const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];
            const dualPatch& dPatch = cTSM.newDualPatch();
            const labelList& mFDFTP = dPatch.mapFromDualFacesToPoints();
            const labelList& mp = triSurf.meshPoints();
  
            forAll(inletBCFaces,I)
            {
                label faceID = includeDualFaces[inletBCFaces[I]];
                label pointID = mFDFTP[faceID]; // this is local 

                newPoints[mp[pointID]].z()=0; // this should be global
            }
            
            triSurf.movePoints(newPoints); // in old OF version, it may use localPoints instead
            
            forAll(addSurfs, movingObjectID)
            {
                const dictionary& movingObjectdict = objectDictList()[addObjectsID[movingObjectID]] ;
                scalar dropVelocity = movingObjectdict.lookupOrDefault<scalar>("dropVelocity",0);
                triSurface& objectSurf = addSurfs[movingObjectID];
                pointField objectPoints = objectSurf.points();
                scalar dt =time().deltaTValue();
                forAll(objectPoints,I)
                {
                    objectPoints[I].z() += dropVelocity*dt;
                }
                objectSurf.movePoints(objectPoints);
                
                Info<<" dropVelocity "<<dropVelocity<<endl;
               
            }
        }
        
     
           
        scalar timeValue = time().timeOutputValue()+1e-6;
        scalar timeValue0 = time().timeOutputValue()-time().deltaT0Value()+1e-6;

        scalar extrudeMeshOutputTimeStep = ibProperties().lookupOrDefault<scalar>("extrudeMeshOutputTimeStep",-1);

        // write out
       if (time().outputTime() 
            or floor(timeValue/extrudeMeshOutputTimeStep)>floor(timeValue0/extrudeMeshOutputTimeStep) )
        {
            volScalarField gradIta=mag(fvc::grad(ita));
            gradIta.rename("beta");
            forAll(gradIta,I)
            {
                gradIta[I]=atan(gradIta[I]);
                gradIta[I]=radToDeg(gradIta[I]);
            }
       
            if(Pstream::parRun())
            {
                if(Pstream::myProcNo()==Pstream::masterNo())
                {
                    
                      parallelWrite(gradIta); // this is different from BETA
                      parallelWrite(qFlux);
                      parallelWrite(dIta);
                    
                      parallelWrite(sN);
                      parallelWrite(sNC);
                      parallelWrite(wss);
                      parallelWrite(ita);
                      parallelWrite(q0);
                      parallelWrite(PHI);
                      parallelWrite(divQ);

                   
                     mesh.setInstance(time().timeName());
                    if(changeSTL)
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
                }
            }
            else
            {

                    
                parallelWrite(gradIta); // this is different from BETA
                parallelWrite(dIta);
            
                parallelWrite(sN);
                parallelWrite(sNC);
                parallelWrite(wss);
                parallelWrite(ita);
                parallelWrite(q0);
                parallelWrite(PHI);
                parallelWrite(divQ);
        
                mesh.setInstance(time().timeName());
                if(changeSTL) mesh.write();
            }
            
            if(changeSTL and Pstream::myProcNo()==Pstream::masterNo())
            {
                fileName path0=time().rootPath()
                    /time().globalCaseName()/time().timeName()/mesh.dbDir()/"triSurface";
                mkDir(path0);
                fileName path=path0/objectNames(objectID)+".stl";

                triSurf.write(path);
                Info<<"write stl to "<<path<<endl;  
                
                
                fileName path1(mesh.time().constant()/"triSurface"/objectNames(objectID)+".stl");
                triSurf.write(path1);
                Info<<"write stl to "<<path1<<endl;
                
                
                 
                forAll(addObjectsName, movingObjectID)
                {
                    
                    fileName movingPath = path0/addObjectsName[movingObjectID]+".stl";
                    
        
                    addSurfs[movingObjectID].write(movingPath);
                    Info<<"write moving stl to "<<movingPath<<endl;
                    
                    fileName movingPath1(mesh.time().constant()/"triSurface"/addObjectsName[movingObjectID]+".stl");
                    
                    
                    addSurfs[movingObjectID].write(movingPath1);
                    Info<<"write moving stl to "<<movingPath1<<endl;
                    
                } 
            }

        }  
        
        // If changeSTL is false, it cannot be changed back to true.
        if(changeSTL and !(time().timeIndex()%label(updateInterval)))
        {
            IBNeedUpdated();
        }
        else
        {
            IBHasUpdated();
        }
    }
    else
    {
        scalar timeValue = time().timeOutputValue()+1e-6;
        scalar timeValue0 = time().timeOutputValue()-time().deltaT0Value()+1e-6;

        scalar extrudeMeshOutputTimeStep = ibProperties().lookupOrDefault<scalar>("extrudeMeshOutputTimeStep",-1);

        // write out
       if (time().outputTime() 
            or floor(timeValue/extrudeMeshOutputTimeStep)>floor(timeValue0/extrudeMeshOutputTimeStep) )
        {
            
            if(Pstream::parRun())
            {
                if(Pstream::myProcNo()==Pstream::masterNo())
                {
                    
                    mesh.setInstance(time().timeName());
                    
                    parallelWriteMesh(mesh.lookupObject<pointIOField>("points"));
                    parallelWriteMesh(mesh.lookupObject<faceCompactIOList>("faces"));
                    parallelWriteMesh(mesh.lookupObject<labelIOList>("owner"));
                    parallelWriteMesh(mesh.lookupObject<labelIOList>("neighbour"));
                    parallelWriteMesh(mesh.lookupObject<polyBoundaryMesh>("boundary"));
                    parallelWriteMesh(mesh.lookupObject<pointZoneMesh>("pointZones"));
                    parallelWriteMesh(mesh.lookupObject<faceZoneMesh>("faceZones"));
                    parallelWriteMesh(mesh.lookupObject<cellZoneMesh>("cellZones"));
                   
                }
            }
            else
            {        
                mesh.setInstance(time().timeName());
                mesh.write();
            }
        }
        
        IBHasUpdated();
    }
}

// also make sure wall patchField type
template<class Type>
void Foam::immersedBoundaryFvMesh::setDualEdgeValueToDualMesh
(
    const label& objectID,
    GeometricField<Type, fvsPatchField, surfaceMesh>& surfValues, // dual Mesh surf
    const Field<Type>& values, // dual patch edge values
    const bool setBC  // default = false   only for flux
)const
{
    const fvMesh& mesh = surfValues.mesh();
    const labelList& newToOldEdgeMap = this->newToOldDualEdgeMap(objectID);

    // fill oldToNewEdgeMap
    Field<Type>& surfValuesI = surfValues.primitiveFieldRef();
   
    forAll(newToOldEdgeMap,I)
    {
        label patchID = mesh.boundaryMesh().whichPatch(I);
        if(patchID==-1)
        {
            surfValuesI[I] = values[newToOldEdgeMap[I]];

        }
        else if(newToOldEdgeMap[I]>-1)
        {
            surfValues.boundaryFieldRef()[patchID][I-mesh.boundary()[patchID].start()] = values[newToOldEdgeMap[I]];
   
        }
    }

    if(setBC)
    {
        forAll(surfValues.boundaryFieldRef(),patchI)
        {
            if (!isA<emptyPolyPatch>(mesh.boundaryMesh()[patchI]))
            {
                surfValues.boundaryFieldRef().set
                (
                    patchI,
                    fvsPatchField<Type>::New("symmetry", mesh.boundary()[patchI], surfValues)
                );
            }
        }
    }


}

template<class Type>
void Foam::immersedBoundaryFvMesh::hitPointExportToDualMesh
(
    const word& psiName,
    const Field<Type>& psi, //hit pt values it can be ib hit or ghost hit
    label objectID,
    bool write //default=true
)const
{

    GeometricField< Type, fvsPatchField, surfaceMesh > PSI
    (
        IOobject
        (
            psiName,
            time().timeName(),
            dualMeshList()[objectID],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMeshList()[objectID],
        dimensioned<Type>
        (
            "0",
            dimless,
            pTraits<Type>::zero
        )
    );

    hitPointExportToDualMesh(psi,PSI,objectID,write);
}

template<class Type>
void Foam::immersedBoundaryFvMesh::hitPointExportToDualMesh
(
    const Field<Type>& psi,
    GeometricField<Type, fvsPatchField, surfaceMesh>& PSI,
    label objectID,
    bool write //default=true
)const
{

    if(IBtypeList()[objectID]=="classic")
    {
        setDualEdgeValueToDualMesh
        (
            objectID,
            PSI,// dualMesh
            mapFromIBHitToDualEdge(psi,objectID)// triFace value
        );
    }
    else if (IBtypeList()[objectID]=="mix")
    {
        setDualEdgeValueToDualMesh
        (
            objectID,
            PSI,// dualMesh
            mapFromGHOSTHitToDualEdge(psi,objectID)// triFace value
        );
    }
    scalar timeValue = time().timeOutputValue()+1e-6;
    scalar timeValue0 = time().timeOutputValue()-time().deltaT0Value()+1e-6;

    scalar extrudeMeshOutputTimeStep = ibProperties().lookupOrDefault<scalar>("extrudeMeshOutputTimeStep",-1);

    // write out
    if (time().outputTime() or floor(timeValue/extrudeMeshOutputTimeStep)>floor(timeValue0/extrudeMeshOutputTimeStep))
    {
         if(write and Pstream::myProcNo()==Pstream::masterNo()) parallelWrite(PSI);
    }
}


// global will be done totally
template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::mapFromDualPatchToTriPoints
(
    const Field<Type>& dualValues,
    label objectID
)const
{
    const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];
    const dualPatch& dPatch = cTSM.newDualPatch();
    const PPatch& oldPatch = dPatch.oldPatch();
    const labelList& mFDFTP = dPatch.mapFromDualFacesToPoints();
    
    //dual patch point values, whole patch
    Field<Type> triPoints(oldPatch.localPoints().size(), pTraits<Type>::zero);
    
    forAll(mFDFTP,I)
    {
        triPoints[mFDFTP[I]] = dualValues[I];
    }


    return triPoints;
}

// global will be done totally
template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::mapFromTriPointsToDualPatch
(
    const Field<Type>& triPoints,
    label objectID
)const
{
    const cutTriSurfaceMesh& cTSM = cutTriSurfaceMeshList()[objectID];
    const dualPatch& dPatch = cTSM.newDualPatch();
    const labelList& mFDFTP = dPatch.mapFromDualFacesToPoints();
    
    //dual patch point values, whole patch
    Field<Type> dualValues(mFDFTP.size(), pTraits<Type>::zero);
    
    forAll(mFDFTP,I)
    {
        dualValues[I] = triPoints[mFDFTP[I]];
    }

    return dualValues;
}
