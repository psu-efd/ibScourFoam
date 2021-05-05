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
    immersedBoundaryFvMeshUpdateCellValues.C
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvMesh.H"
#include "MeshWave.H"
#include "SortableList.H"
#include "fvCFD.H"
#include "turbulenceModel.H"
#include "wallFuncModules.H"
#include "turbulentTransportModel.H"
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::immersedBoundaryFvMesh::initialIB() const
{
    forAll(objectsList(),objectID)
    {
        
            const triSurfaceSearch& tss = triSurfaceSearchList()[objectID];

            word Inlet_ = objectDictList()[objectID].lookupOrDefault<word>("InletType","notfixedValue");

            word wallFunc_ = ibProperties().lookupOrDefault<word>("wallFunction","vanDriest");
            wallFunc_ = objectDictList()[objectID].lookupOrDefault<word>("wallFunction",wallFunc_);
            if(Inlet_=="fixedValue")
            {
                wallFuncModules wFM(wallFunc_,"Menter2001");

                volVectorField& U = const_cast<volVectorField&>
                (this->lookupObject<volVectorField>("U"));

                volScalarField& omega = const_cast<volScalarField&>
                (this->lookupObject<volScalarField>("omega"));

                volScalarField& k = const_cast<volScalarField&>
                (this->lookupObject<volScalarField>("k"));

                volScalarField& nut = const_cast<volScalarField&>
                (this->lookupObject<volScalarField>("nut"));


                volScalarField& Gamma = const_cast<volScalarField&>
                    (this->lookupObject<volScalarField>("Gamma"));

                scalar Ks_ = ibProperties().lookupOrDefault<scalar>("roughnessHeight",0);
                Ks_ = objectDictList()[objectID].lookupOrDefault<scalar>("roughnessHeight",Ks_);
                
                
                 
                scalar Intensity_ = ibProperties().lookupOrDefault<scalar>
                    ("turbulentIntensity",0.05);
                Intensity_ = objectDictList()[objectID].lookupOrDefault<scalar>
                    ("turbulentIntensity",Intensity_);

                const dictionary& transportProperties =
                      this->lookupObject<dictionary>
                      (
                         "transportProperties"
                      );
                dimensionedScalar nu(transportProperties.lookup("nu"));

                scalar nuLam=nu.value();

                boundBox* meshBB = new boundBox(this->points(), false);
                vector    span = meshBB->span();
                delete meshBB;
                scalar waterDepth = meshBB->max()[2];
                waterDepth = ibProperties().lookupOrDefault<scalar>
                    ("waterDepth",waterDepth);
                waterDepth = objectDictList()[objectID].lookupOrDefault<scalar>
                    ("waterDepth",waterDepth);

                scalar boundaryLayerThickness = ibProperties().lookupOrDefault<scalar>
                    ("boundaryLayerThickness",waterDepth);
                boundaryLayerThickness = objectDictList()[objectID].lookupOrDefault<scalar>
                    ("boundaryLayerThickness",boundaryLayerThickness);
                vector meanVelocity=objectDictList()[objectID].lookupOrDefault<vector>
                    ("meanVelocity",pTraits<vector>::zero);
                if(debug)
                {
                    Info<<"  Given"<<endl
                        <<tab<<"water depth: "<<waterDepth<<" m"<<endl
                        <<tab<<"boundary layer thickness: "<<boundaryLayerThickness<<" m"<<endl
                        <<tab<<"mean velocity: "<<meanVelocity<<" m/s"<<endl;
                    if(Ks_>0)
                    {
                        Info<<tab<<"Roughness Height: "<<Ks_*100<<" cm" <<endl;
                    }
                    else
                    {
                        Info<<tab<<"Bed is smooth"<<endl;
                    }
                }
                // estimate utau for give water depth and mean velocity
               
                scalar U_edge=meanVelocity[0];
                scalar uTau = U_edge/( 5.75 * log10(12.2*waterDepth/Ks_));
		        Info<<"Estimate shear velocity: "<<uTau*100<<" cm/s"<<endl;

                uTau = objectDictList()[objectID].lookupOrDefault<scalar>
                    ("uTau",uTau);
                Info<<"Estimate shear velocity: "<<uTau*100<<" cm/s"<<endl;
                
                scalar ksPlus = Ks_*uTau/nuLam;
                ksPlus = max(1.0,ksPlus);
                
                forAll (U,I)
                {
                    vector& URef=U[I];
                    scalar& kRef=k[I];
                    scalar& omegaRef=omega[I];
                    scalar& nutwb=nut[I];

                    pointIndexHit pih = tss.nearest
                        (this->cellCentres()[I], span);
                    scalar y = 0.0;
                    if(pih.hit())
                    {
                        point nearestPoint = pih.hitPoint();
                        y = mag(this->cellCentres()[I]-nearestPoint);
                    }

                   
                    URef*=0;

                    scalar yPlus = y*uTau/nuLam;
                    yPlus = max(1.0,yPlus);
                    scalar uPlus = wFM.uPlus(yPlus,ksPlus);
                    URef.x() = uPlus*uTau;

                    
                    nutwb = max(0.0, (1.0/wFM.dupdyp(yPlus,ksPlus)-1.0)*nuLam);

                    wFM.kOmegaUpdate
                        (
                            yPlus,
                            uTau,
                            ksPlus,
                            kRef,
                            omegaRef,
                            nutwb
                        );   
                    if(Gamma[I]<0)
                    {
                        U[I]*=0;
                    }
                }
                evaluateCoupled(k);

                evaluateCoupled(omega);

                evaluateCoupled(nut);

                
            }
          
        
    }

}

void Foam::immersedBoundaryFvMesh::setInlet() const
{
    forAll(objectsList(),objectID)
    {
        
            const triSurfaceSearch& tss = triSurfaceSearchList()[objectID];

            word Inlet_ = objectDictList()[objectID].lookupOrDefault<word>("InletType","notfixedValue");

            word wallFunc_ = ibProperties().lookupOrDefault<word>("wallFunction","vanDriest");
            wallFunc_ = objectDictList()[objectID].lookupOrDefault<word>("wallFunction",wallFunc_);
            if(Inlet_=="fixedValue")
            {
                wallFuncModules wFM(wallFunc_,"Menter2001");

                volVectorField& U = const_cast<volVectorField&>
                (this->lookupObject<volVectorField>("U"));

                volVectorField::Boundary& Upatches = U.boundaryFieldRef();

         
                label patchID = this->boundaryMesh().findPatchID("inlet");

                volScalarField& Gamma = const_cast<volScalarField&>
                    (this->lookupObject<volScalarField>("Gamma"));

                scalar Ks_ = ibProperties().lookupOrDefault<scalar>("roughnessHeight",0);
                Ks_ = objectDictList()[objectID].lookupOrDefault<scalar>("roughnessHeight",Ks_);
                Ks_ =  ibProperties().lookupOrDefault<scalar>("roughnessHeight1",Ks_);

                scalar Intensity_ = ibProperties().lookupOrDefault<scalar>
                    ("turbulentIntensity",0.05);
                Intensity_ = objectDictList()[objectID].lookupOrDefault<scalar>
                    ("turbulentIntensity",Intensity_);

                const dictionary& transportProperties =
                      this->lookupObject<dictionary>
                      (
                         "transportProperties"
                      );
                dimensionedScalar nu(transportProperties.lookup("nu"));

                
                vector gravity= objectDictList()[objectID].lookup("gravity");
                vector gDir = gravity/mag(gravity);

                boundBox* meshBB = new boundBox(this->points(), false);

                scalar waterDepth = meshBB->max()[2];
                waterDepth = ibProperties().lookupOrDefault<scalar>
                    ("waterDepth",waterDepth);
                waterDepth = objectDictList()[objectID].lookupOrDefault<scalar>
                    ("waterDepth",waterDepth);

                scalar boundaryLayerThickness = ibProperties().lookupOrDefault<scalar>
                    ("boundaryLayerThickness",waterDepth);
                boundaryLayerThickness = objectDictList()[objectID].lookupOrDefault<scalar>
                    ("boundaryLayerThickness",boundaryLayerThickness);
                vector meanVelocity = objectDictList()[objectID].lookupOrDefault<vector>
                    ("meanVelocity",pTraits<vector>::zero);
                vector desiredFlowRate = objectDictList()[objectID].lookupOrDefault<vector>
                    ("desiredFlowRate",pTraits<vector>::zero);


                scalar sumArea = 0.0;
                scalarField IBWallDistance(Upatches[patchID].size(),-1.0);

                const pointField& fstart =  Upatches[patchID].patch().Cf();
                pointField fend = fstart + this->bounds().mag()*gDir;
                
                List<pointIndexHit> fhitInfo;
                tss.findLine(fstart, fend, fhitInfo);
                
                forAll (Upatches[patchID], I)
                {
                    label faceCellID = this->boundary()[patchID].faceCells()[I];
                    const pointIndexHit& pih = fhitInfo[I];
                    if (pih.hit())
                    {
                        IBWallDistance[I] = mag(fstart[I]-pih.hitPoint());

                        if(IBWallDistance[I]>cellSize(faceCellID)*0.5)
                        {
                            sumArea+=this->boundary()[patchID].magSf()[I];
                        }
                        else
                        {
                            sumArea+=this->boundary()[patchID].magSf()[I]
                                    *IBWallDistance[I]/cellSize(faceCellID);
                        }
                    }
                }
                reduce(sumArea, sumOp<scalar>());
                if(debug)
                {
                    Info<<"  Given"<<endl
                        <<tab<<"desired water depth: "<<waterDepth<<" m"<<endl
                        <<tab<<"desired boundary layer thickness: "<<boundaryLayerThickness<<" m"<<endl;
                    if(mag(desiredFlowRate)>0) // fixed flow rate
                    {
                        meanVelocity = desiredFlowRate/sumArea;
                        Info<<tab<<"desired fixed flow rate: "<<desiredFlowRate<<" m3/s"<<endl

                            <<tab<<"desired mean velocity: "<<meanVelocity<<" m/s"<<endl;
                    }
                    else // fixed mean velocity
                    {
                        desiredFlowRate = meanVelocity*sumArea;
                        Info<<tab<<"desired fixed mean velocity: "<<meanVelocity<<" m/s"<<endl

                            <<tab<<"desired flow rate: "<<desiredFlowRate<<" m3/s"<<endl;
                    }
                    if(Ks_>0)
                    {
                        Info<<tab<<"Roughness Height: "<<Ks_*100<<" cm" <<endl;
                    }
                    else
                    {
                        Info<<tab<<"Bed is smooth"<<endl;
                    }
                }
                // estimate utau for give water depth and mean velocity
                scalar U_edge = meanVelocity[0];

                forAll (Upatches[patchID], I)
                {
                    label faceCellID = this->boundary()[patchID].faceCells()[I];

                    vector& URef = Upatches[patchID][I];
  
                    URef.x() = U_edge;

                      
                    if(Gamma[faceCellID]<0.0)
                    {
                        URef.x() *=-0;
                    }
                     
                }


            }
        
    }

}



void Foam::immersedBoundaryFvMesh::evaluateU() const
{

    volVectorField& U = const_cast<volVectorField&>
    (this->lookupObject<volVectorField>("U"));
    vectorField& UI = U.primitiveFieldRef();

    word nutWORD = "nut";
    dictionary turbDict (this->lookupObject<dictionary>("turbulenceProperties"));
    word modelType = turbDict.lookup("simulationType");
    if(modelType == "laminar")
    {
        nutWORD = "nu";
    }

    volScalarField& nut = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>(nutWORD));

    surfaceScalarField& phi = const_cast<surfaceScalarField&>
        (this->lookupObject<surfaceScalarField>("phi"));

    volScalarField& Gamma = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("Gamma"));

    forAll(objectsList(),objectID)
    {
        BCtype_ = objectDictList()[objectID].lookupOrDefault<word>("BCtype","fixedWall");

        if(IBtypeList()[objectID]=="classic")
        {
            const labelList& ibCells=ibCellsList()[objectID];

            labelHashSet ibCellsSet(ibCells);

            // check every cell to see if it change from dry to wet
            //- Gamma live cells (1), IB cells (0.5), ghost cells (-0.5), and dead cells(-1)
            forAll(Gamma,I)
            {
                scalar oldGamma = Gamma.oldTime()[I];
                scalar newGamma = Gamma[I];
                if(oldGamma<0.5 and newGamma>0) // from dead cell to IB/live cell
                {
                    U[I] *=0;
                    nut[I] *=0;
                }
                
            }
            // set deadCell velocity to be zero
            const labelList& ibDeadCells(ibDeadCellsList()[objectID]);
            forAll(ibDeadCells,I)
            {
                label cellID=ibDeadCells[I];
                UI[cellID]=UI[cellID]*0;
                nut[cellID] *=0;
            }
            labelHashSet ibDeadCellsSet(ibDeadCells);
            // Evaluate the uncoupled patchField
            volVectorField::Boundary& Upatches = U.boundaryFieldRef();

            forAll (Upatches, patchI)
            {
                if (!Upatches[patchI].coupled())
                {
                    forAll(Upatches[patchI],I)
                    {
                        
                        label faceCellID= this->boundary()[patchI].faceCells()[I];
                        if(ibDeadCellsSet.found(faceCellID))
                        {
                            Upatches[patchI][I]=UI[faceCellID];
                        }
                    }
                }
            }

            phi=fvc::flux(U);
            ibCorrectPhi(phi);
        }
        else if(IBtypeList()[objectID]=="mix")
        {
            // set deadCell velocity to be zero
            const labelList& ibDeadCells(ibDeadCellsList()[objectID]);
            labelHashSet ghostCellsSet(ghostCellsList()[objectID]);
            forAll(ibDeadCells,I)
            {
                label cellID=ibDeadCells[I];
                if(!ghostCellsSet.found(cellID))
                {
                    UI[cellID]=UI[cellID]*0;
                    nut[cellID] *=0;
                }
            }
        }
    }
    setInlet();// calculate wall distance + velocity profile
    evaluateCoupled(U);

}

void Foam::immersedBoundaryFvMesh::manipulateUEqn
(
    fvVectorMatrix& UEqn
) const
{
    // boundary condition type
    forAll(objectsList(),objectID)
    {
        BCtype_ = objectDictList()[objectID].lookupOrDefault<word>("BCtype","fixedWall");

        if(IBtypeList()[objectID]=="classic")
        {

            ibCellReconstruction(UEqn,objectID);

        }
        else if(IBtypeList()[objectID]=="ghost-cell")
        {
            ghostCellReconstruction(UEqn,objectID);
        }
        else if(IBtypeList()[objectID]=="mix")
        {
            ghostCellReconstruction(UEqn,objectID);
        }
    }
}












void Foam::immersedBoundaryFvMesh::ibCellForcing
(
    fvVectorMatrix& UEqn,
    label objectID
) const
{
    volVectorField& U = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("U"));

    volVectorField& U_desired = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("U_desired"));


    dictionary turbDict (this->lookupObject<dictionary>("turbulenceProperties"));
    word modelType=turbDict.lookup("simulationType");
    if(modelType == "laminar")
    {
        fromSPointReconstruction(U,objectID);
    }
    else
    {
        surfaceScalarField& phi = const_cast<surfaceScalarField&>
            (this->lookupObject<surfaceScalarField>("phi"));

        volVectorField U_bak(U);
        U = U_desired;
        volVectorField bodyForce = fvc::ddt(U);
        U = U_bak;

        typedef incompressible::turbulenceModel incTub;
        const incTub& turbulence = const_cast<incTub&>
        (
            this->lookupObject<incTub>
            (
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    U.group()
                )
            )
        );
        surfaceScalarField nuEff
        (
            "nuEff",
            fvc::interpolate(turbulence.nuEff())
        );
        volVectorField& gradP = const_cast<volVectorField&>
            (this->lookupObject<volVectorField>("gradP"));
        volScalarField& p = const_cast<volScalarField&>
            (this->lookupObject<volScalarField>("p"));
        volScalarField p_bak=p;
        this->evaluateP();
        bodyForce +=
            (
                fvc::div(phi, U)
                - fvc::laplacian(nuEff, U)
                + gradP
            );
        p = p_bak;
        labelHashSet globalLiveCellsSet(globalLiveCells());
        forAll(bodyForce,cellI)
        {
            if(globalLiveCellsSet.found(cellI))
            {
                bodyForce[cellI] *=0;
            }
        }

        UEqn=UEqn-bodyForce;
        bodyForce.rename("bodyForce");
        if (time().outputTime())
        {
            bodyForce.write();
            turbulence.nuEff()->write();
        }
    }
}

void Foam::immersedBoundaryFvMesh::ibCellReconstruction
(
    fvVectorMatrix& UEqn,
    label objectID
) const
{
    dictionary turbDict (this->lookupObject<dictionary>("turbulenceProperties"));

    volVectorField& U = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("U"));

    evaluateCoupled(U);

    word modelType = turbDict.lookup("simulationType");



    manipulateMatrix(UEqn,objectID);

}











void Foam::immersedBoundaryFvMesh::ghostCellReconstruction
(
    fvVectorMatrix& UEqn,
    label objectID
) const
{
    dictionary turbDict (this->lookupObject<dictionary>("turbulenceProperties"));

    volVectorField& U = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("U"));

    volVectorField& U_desired = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("U_desired"));

    volScalarField& p = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("p"));

    surfaceScalarField& phi = const_cast<surfaceScalarField&>
        (this->lookupObject<surfaceScalarField>("phi"));

    const labelList& gCL = ghostCellsList()[objectID];
    word modelType = turbDict.lookup("simulationType");
    if(modelType == "laminar")
    {
        const pointField& imagePoints = imagePointsList()[objectID];
        const pointField& ghostHitPoints = ghostHitPointsList()[objectID];
        const pointField gc_centers(this->C(),ghostCellsList()[objectID]);


        vectorField UIPV(imagePointsValues(U,objectID));

        scalarField pIPV(imagePointsValues(p,objectID));

        pointField ghostCellCentres(this->C(),gCL);

        forAll(UIPV,I)
        {
            label cellID = gCL[I];
            vector Uw=pTraits<vector>::zero;
            const vector& Ui=UIPV[I];

            scalar Liw=mag(imagePoints[I]-ghostHitPoints[I]);
            scalar Ldw=mag(ghostCellCentres[I]-ghostHitPoints[I]);

            U_desired[cellID] = (Uw-Ui)/Liw*Ldw+Uw;
            p[cellID] = pIPV[I];
        }

        UIPV.clear();
        pIPV.clear();
    }

    volVectorField bodyForce
    (
        IOobject
        (
            "bodyForce",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedVector("zero", dimVelocity/dimTime, vector::zero)
    );

    const turbulenceModel& turbulence = this->lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            U.group()
        )
    );


    volVectorField U_bak(U);
    U = U_desired;
    bodyForce = fvc::ddt(U);
    U = U_bak;

    labelHashSet ghostCellsSet(gCL);

    surfaceScalarField nuEff
    (
        "nuEff",
        fvc::interpolate(turbulence.nuEff())
    );


    bodyForce +=
        (
            fvc::div(phi, U)
            - fvc::laplacian(nuEff, U)
            + fvc::grad(p)
        );

    forAll(bodyForce,cellI)
    {
        if(!ghostCellsSet.found(cellI))
        {
            bodyForce[cellI] *=0;
        }
    }

    UEqn=UEqn-bodyForce;
    if (time().outputTime())bodyForce.write();

    bodyForce.clear();

    // fix deadCells excluding ghost cells
    scalarField& Diag = UEqn.diag();

    const labelList& dce_org = ibDeadCellsList()[objectID]; // dead cell include ghost cells
    labelHashSet dceSet;
    forAll(dce_org,I)
    {
        label cellID=dce_org[I];
        if(!ghostCellsSet.found(cellID) and !dceSet.found(cellID))
        {
            dceSet.insert(cellID);
        }
    }
    labelList dce=dceSet.toc();
    // Estimate diagonal in live cells
    scalar liveDiag = 1;

    scalar gSumMagDiag = gSumMag(Diag);
    scalar gMaxDiag = gMax(Diag);

    if ( (dce.size()) < Diag.size())
    {
    // Problem: if dce.size() == Diag.size() in one processor, gSumMag(Diag) here may not be solved.
    // added by Y.C. Xu 2017/6
        liveDiag = gSumMagDiag/(Diag.size() - dce.size());

        // Correct for sign
        liveDiag *= sign(gMaxDiag);

    }

    forAll (dce, cellI)
    {
        if (mag(Diag[dce[cellI]]) < SMALL)
        {
            Diag[dce[cellI]] = liveDiag;
        }
    }

    vectorField deadCellsPsi
    (
        UEqn.psi(),
        dce
    );

}

template<class Type>
void Foam::immersedBoundaryFvMesh::fromSPointReconstruction
(
    GeometricField<Type,fvPatchField,volMesh>& psi,
    label objectID
) const
{
    GeometricField<Type,fvPatchField,volMesh>& U=psi;
    Field<Type>& UI = U.primitiveFieldRef();
    const pointField& samplingPoints = samplingPointsList()[objectID];
    const pointField& ibHitPoints = ibHitPointsList()[objectID];
    const pointField ibc_centers(this->C(),ibCellsList()[objectID]);
    const labelList& ibCells = ibCellsList()[objectID];

    label counter = 0;
    scalarField error(samplingPoints.size(), 0);
    scalar maxError = 0;

    vectorField ibn(ibc_centers-ibHitPoints);
    scalarField yIB(mag(ibn));
    ibn=ibn/(SMALL+mag(ibn));
    scalarField ySample(mag(samplingPoints-ibHitPoints));

    do
    {
        counter++;

        scalar maxMagUI = 0;

        Field<Type> samplingPointsValues
                (
                    stencilInterpolation
                    (
                        UI,
                        objectID,
                        ibCells,
                        samplingStencils(objectID).cellCells(),
                        samplingStencils(objectID).cellProcCells(),
                        samplingStencils(objectID).procCells(),
                        samplingStencils(objectID).procCentres(),
                        samplingStencils(objectID).weights(),
                        samplingStencils(objectID).procWeights()
                    )
                );

        forAll(ibCells, I)
        {
            label cellID = ibCells[I];
            Type oldUI = UI[cellID];
            Type ibValue = oldUI-oldUI;
            Type samplingValue = samplingPointsValues[I];
            UI[cellID] = (samplingValue-ibValue)/ySample[I]*yIB[I]+ibValue;
            Type tempUI = UI[cellID]-(UI[cellID]&ibn[I])*ibn[I];
            UI[cellID]  = tempUI/(mag(tempUI)+SMALL)*mag(UI[cellID]);
            error[I] = mag(UI[cellID] - oldUI);
            if(maxMagUI<mag(UI[cellID])) maxMagUI = mag(UI[cellID]);
        }

        error /= maxMagUI + SMALL;

        maxError = gMax(error);

        if(debug)
        {
            Info<<"Reconstruct "<<psi.name()<<", maxError: "<<maxError<<" counter: "<<counter<<endl;
        }
    }
    while (maxError > 0.0001 && counter < 10);
}

void Foam::immersedBoundaryFvMesh::fromIbFaceSPointReconstruction
(
    surfaceScalarField& phi,
    const volVectorField& U,
    label objectID
) const
{
    scalarField& phiI = phi.primitiveFieldRef();
    const pointField& ibFaceSamplingPoints = ibFaceSamplingPointsList()[objectID];
    const pointField& ibFaceHitPoints = ibFaceHitPointsList()[objectID];
    const pointField ibf_centers(this->faceCentres(),ibFacesList()[objectID]);
    const labelList& ibFaces = ibFacesList()[objectID];

    label counter = 0;
    scalarField error(ibFaceSamplingPoints.size(), 0);
    scalar maxError = 0;

    vectorField ibn(ibf_centers-ibFaceHitPoints);
    scalarField yIB(mag(ibn));
    ibn=ibn/(SMALL+mag(ibn));
    scalarField ySample(mag(ibFaceSamplingPoints-ibFaceHitPoints));

    do
    {
        counter++;

        scalar maxMagphiI = 0;

        vectorField ibFaceSamplingPointsValues
                (
                    stencilInterpolation
                    (
                        U,
                        objectID,
                        ibFaces,
                        ibFaceSamplingStencils(objectID).cellCells(),
                        ibFaceSamplingStencils(objectID).cellProcCells(),
                        ibFaceSamplingStencils(objectID).procCells(),
                        ibFaceSamplingStencils(objectID).procCentres(),
                        ibFaceSamplingStencils(objectID).weights(),
                        ibFaceSamplingStencils(objectID).procWeights()
                    )
                );

        forAll(ibFaces, I)
        {
            label faceID = ibFaces[I];
            scalar oldphiI = phiI[faceID];
            vector ibValue = pTraits<vector>::zero;
            vector samplingValue = ibFaceSamplingPointsValues[I];
            vector tempUI = (samplingValue-ibValue)/ySample[I]*yIB[I]+ibValue;
            phiI[faceID]  = tempUI&this->faceAreas()[faceID];

            error[I] = mag(phiI[faceID] - oldphiI);
            if(maxMagphiI<mag(phiI[faceID])) maxMagphiI = mag(phiI[faceID]);
        }

        error /= maxMagphiI + SMALL;

        maxError = gMax(error);

        if(debug)
        {
            Info<<"Reconstruct "<<phi.name()<<", maxError: "<<maxError<<" counter: "<<counter<<endl;
        }
    }
    while (maxError > 0.0001 && counter < 10);
}

template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::stencilInterpolation
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
) const
{
    if (psi.size() != this->nCells())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::Field<Type>\n"
            "immersedBoundaryFvMesh::stencilInterpolation\n"
            "(\n"
            "    const Field<Type>& psi\n"
            "    label objectID,\n"
            ") const"
        )   << "Field size does not correspond to cell centres "
            << "for patch " << ibProperties().subDict("objects").toc()[objectID] << nl
            << "Field size = " << psi.size()
            << " nCells = " << this->nCells()
            << abort(FatalError);
    }

    Field<Type> ibPsi(ibc.size(), pTraits<Type>::zero);

    labelHashSet ibcSet(ibc);
    // Do interpolation, local cell data
    forAll (ibc, cellI)
    {
        const labelList& curAddr = ibcc[cellI];
        const scalarList& curWeights = cellWeights[cellI];

        forAll (curAddr, ccI)
        {
            ibPsi[cellI] += curWeights[ccI]*psi[curAddr[ccI]];

        }
    }

    List<labelHashSet> ibcSetList(Pstream::nProcs());
    labelListList ibcList(Pstream::nProcs());


    // Parallel communication for psi
    FieldField<Field, Type> procCellValues = sendAndReceive(psi,procCells,CProc);

    // Do interpolation, cell data from other processors
    forAll (ibc, cellI)
    {
        const List<labelPair>& curProcCells = ibcProcC[cellI];
        const scalarList& curProcWeights = cellProcWeights[cellI];

        forAll (curProcCells, cpcI)
        {
            ibPsi[cellI] +=
                curProcWeights[cpcI]*
                procCellValues
                [
                    curProcCells[cpcI].first()
                ]
                [
                    curProcCells[cpcI].second()
                ];

        }
        if((ibcc[cellI].size()+ibcProcC[cellI].size()) < SMALL)
        {
            ibPsi[cellI] = SMALL* pTraits<Type>::one + pTraits<Type>::zero;
        }
    }

    return ibPsi;
}

template<class Type>
Foam::FieldField<Foam::Field, Type>
Foam::immersedBoundaryFvMesh::sendAndReceive
(
    const Field<Type>& psi,
    const labelListList& procCells,
    const vectorListList& CProc
) const
{


    FieldField<Field, Type> procPsi(Pstream::nProcs());

    forAll (procPsi, procI)
    {
        procPsi.set
        (
            procI,
            new Field<Type>
            (
                CProc[procI].size(),
                pTraits<Type>::zero
            )
        );
    }

    if (Pstream::parRun())
    {
        // Send
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Do not send empty lists
                if (!procCells[procI].empty())
                {
                    Field<Type> curPsi(psi, procCells[procI]);

                    // Parallel data exchange
                    OPstream toProc
                    (
                        Pstream::commsTypes::blocking,
                        procI,
                        curPsi.size()*sizeof(Type)

                    );
                    toProc << curPsi;
                }
            }
        }
        // Receive
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {

            if (procI != Pstream::myProcNo())
            {

                // Do not receive empty lists
                if (!procPsi[procI].empty())
                {

                    // Parallel data exchange
                    IPstream fromProc
                    (
                        Pstream::commsTypes::blocking,
                        procI,
                        procPsi[procI].size()*sizeof(Type)
                    );

                    fromProc >> procPsi[procI];

                }
            }
        }
    }
    return procPsi;
};





void immersedBoundaryFvMesh::manipulateMatrix
(
    fvVectorMatrix& eqn,
    label objectID
) const
{
    const labelList& dce = ibDeadCellsList()[objectID];
    const labelList& ibc = ibCellsList()[objectID];

    correctDiag(eqn,objectID);

    vectorField polyPsi(eqn.psi(), ibCellsList()[objectID]);
    eqn.setValues(ibc, polyPsi);
    
 
    // Correct equation for dead cells
    vectorField deadCellsPsi
    (
        eqn.psi(),
        dce
    );
    deadCellsPsi = deadCellsPsi*0;

    eqn.setValues(dce, deadCellsPsi);
    
}


void immersedBoundaryFvMesh::manipulateMatrix
(
    fvScalarMatrix& eqn,
    label objectID
) const
{
    const labelList& dce = ibDeadCellsList()[objectID];
    const labelList& ibc = ibCellsList()[objectID];

    correctDiag(eqn,objectID);


    scalarField polyPsi(eqn.psi(), ibCellsList()[objectID]);
    
    eqn.setValues(ibc, polyPsi);
    

   
    // Correct equation for dead cells
    scalarField deadCellsPsi
    (
        eqn.psi(),
        dce
    );
    
    deadCellsPsi = deadCellsPsi*0;

    eqn.setValues(dce, deadCellsPsi);
    
}


template<class Type>
void immersedBoundaryFvMesh::correctDiag
(
    fvMatrix<Type>& eqn,
    label objectID
) const
{
    scalarField& Diag = eqn.diag();

    const labelList& dce = ibDeadCellsList()[objectID];
    const labelList& ibc = ibCellsList()[objectID];

    // Estimate diagonal in live cells
    scalar liveDiag = 1;

    scalar gSumMagDiag = gSumMag(Diag);
    scalar gMaxDiag = gMax(Diag);

    if ( (dce.size()+ibc.size()) < Diag.size())
    {
    // Problem: if dce.size() == Diag.size() in one processor, gSumMag(Diag) here may not be solved.
    // added by Y.C. Xu 2017/6

        liveDiag = gSumMagDiag/(Diag.size() - dce.size());

        // Correct for sign
        liveDiag *= sign(gMaxDiag);

    }

    forAll (dce, cellI)
    {
        if (mag(Diag[dce[cellI]]) < SMALL)
        {
            Diag[dce[cellI]] = liveDiag;
        }
    }

 
}


template<class Type>
void immersedBoundaryFvMesh::ghostCorrectDiag
(
    fvMatrix<Type>& eqn,
    label objectID
) const
{
    scalarField& Diag = eqn.diag();

    const labelList& dce = ghostCellsList()[objectID]; // dead cell include ghost cells

    // Estimate diagonal in live cells
    scalar liveDiag = 1;

    scalar gSumMagDiag = gSumMag(Diag);
    scalar gMaxDiag = gMax(Diag);

    if ( (dce.size()) < Diag.size())
    {
    // Problem: if dce.size() == Diag.size() in one processor, gSumMag(Diag) here may not be solved.
    // added by Y.C. Xu 2017/6
        liveDiag = gSumMagDiag/(Diag.size() - dce.size());

        // Correct for sign
        liveDiag *= sign(gMaxDiag);

    }

    forAll (dce, cellI)
    {
        if (mag(Diag[dce[cellI]]) < SMALL)
        {
            Diag[dce[cellI]] = liveDiag;
        }
    }

    // Correct equation for dead cells

    Field<Type> deadCellsPsi
    (
        eqn.psi(),
        dce
    );

//    eqn.setValues(dce, deadCellsPsi);
}

void immersedBoundaryFvMesh::correctOffDiag
(
    fvScalarMatrix& eqn,
    label objectID
) const
{
    const labelList& ibFaces = ibFacesList()[objectID];
    const labelList& ibFaceCells = ibFaceCellsList()[objectID];
    const scalarField& ibGamma = ibGammaList()[objectID].primitiveField();

    const unallocLabelList& own = this->owner();
    const unallocLabelList& nei = this->neighbour();

    // Get delta coefficients
    const surfaceScalarField& dc = this->deltaCoeffs();
    const scalarField& dcI = dc.internalField();

    const surfaceVectorField& Sf = this->Sf();
    const surfaceScalarField& magSf = this->magSf();

    // assign new values of ibFaceGrad
    vectorField ibFaceGrad(ibFaces.size(), pTraits<vector>::zero);

    const volScalarField& psiI = eqn.psi();

    word name = psiI.internalField().name();

    volScalarField& psi = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>(name));

     // calculate the gradient on the internal faces
    const vectorField gradFieldSI = linearInterpolate(fvc::grad(psi));

     // calculate the gradient on the proc faces
    forAll(psi.boundaryFieldRef(),patchi)
    {
        if(psi.boundaryFieldRef()[patchi].coupled())
        {
            psi.boundaryFieldRef()[patchi].initEvaluate(Pstream::commsTypes::blocking);// not quite sure
        }
    }
    forAll(psi.boundaryFieldRef(),patchi)
    {
        if(psi.boundaryFieldRef()[patchi].coupled())
        {
            psi.boundaryFieldRef()[patchi].evaluate(Pstream::commsTypes::blocking);// not quite sure
        }
    }
    if(name!="p")
    {
        forAll (ibFaces, faceI)
        {
            const label curFace = ibFaces[faceI];

            if (curFace < nei.size())
            {
                ibFaceGrad[faceI] = gradFieldSI[curFace];
            }
        }
    }

    if (eqn.symmetric())
    {
        scalarField& diag = eqn.diag();
        scalarField& upper = eqn.upper();
        scalarField& source = eqn.source();

        forAll (ibFaces, faceI)
        {
            const label curFace = ibFaces[faceI];

            if (curFace < nei.size())
            {

                // calculate surface vector on ibFaces    added by Xu 7/2017
                const vector SfN = Sf[curFace]/magSf[curFace];

                // calculate value gradient on ibFaces o->n     added by Xu 7/2017
                const scalar faceGrad = (SfN&ibFaceGrad[faceI]) * pTraits<scalar>::one;

                // Internal face.  One side is an ibCell and another is a
                // live cell. Add gradient to the source of the live cell
                // and kill the off-diagonal coefficient
                if (ibGamma[own[curFace]] > SMALL)
                {
                    diag[own[curFace]] += upper[curFace];


                    source[own[curFace]] -=
                        upper[curFace]*faceGrad/dcI[curFace];

                }
                else
                {
                    diag[nei[curFace]] += upper[curFace];


                    source[nei[curFace]] +=
                        upper[curFace]*faceGrad/dcI[curFace];

                }

                upper[curFace] = 0.0;
            }
            else
            {

                label patchi = this->boundaryMesh().whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchi].empty())
                {
                    label patchFacei =
                        this->boundaryMesh()[patchi].whichFace(curFace);

                    const fvPatchField<scalar>& ptf = psi.boundaryField()[patchi];
                    const tmp<Field<scalar> > tpnf = ptf.patchNeighbourField();
                    const scalar& pnff = tpnf()[patchFacei];

                    tmp<Field<scalar> > tpic = eqn.internalCoeffs()[patchi].component(0);
                    const scalar picf = tpic()[patchFacei];

                    tmp<Field<scalar> > tpbc = eqn.boundaryCoeffs()[patchi];
                    const scalar pbcf = tpbc()[patchFacei];

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<scalar>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<scalar>::zero;

                    // Check if the live cell is on local or neighbour side
                    if (ibFaceCells[faceI] > -1)
                    {
                        if (ibGamma[ibFaceCells[faceI]] > SMALL)
                        {
                            diag[ibFaceCells[faceI]] += picf;

                            source[ibFaceCells[faceI]] +=
                                 cmptMultiply(pbcf, pnff);
                        }
                    }
                }
            }
        }
    }
    else if (eqn.asymmetric())
    {
        scalarField& diag = eqn.diag();
        scalarField& upper = eqn.upper();
        scalarField& lower = eqn.lower();
        Field<scalar>& source = eqn.source();

        forAll (ibFaces, faceI)
        {
            const label curFace = ibFaces[faceI];

            if (curFace < nei.size())
            {
                // calculate surface vector on ibFaces    added by Xu 7/2017
                const vector SfN = Sf[curFace]/magSf[curFace];

                // calculate value gradient on ibFaces o->n     added by Xu 7/2017
                const scalar faceGrad = (SfN&ibFaceGrad[faceI]) * pTraits<scalar>::one;

                // Internal face.  One side is an ibCell and another is a
                // live cell. Add gradient to the source of the live cell
                // and kill the off-diagonal coefficient
                if (ibGamma[own[curFace]] > SMALL)
                {
                    diag[own[curFace]] += upper[curFace];

                    source[own[curFace]] -=
                        upper[curFace]*faceGrad/dcI[curFace];
                }
                else
                {
                    diag[nei[curFace]] += lower[curFace];

                    source[nei[curFace]] +=
                        lower[curFace]*faceGrad/dcI[curFace];
                }

                upper[curFace] = 0.0;
                lower[curFace] = 0.0;
            }
            else
            {

                label patchi = this->boundaryMesh().whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchi].empty())
                {
                    label patchFacei =
                        this->boundaryMesh()[patchi].whichFace(curFace);

                    const fvPatchField<scalar>& ptf = psi.boundaryField()[patchi];
                    const tmp<Field<scalar> > tpnf = ptf.patchNeighbourField();
                    const scalar& pnff = tpnf()[patchFacei];

                    tmp<Field<scalar> > tpic = eqn.internalCoeffs()[patchi].component(0);
                    const scalar picf = tpic()[patchFacei];

                    tmp<Field<scalar> > tpbc = eqn.boundaryCoeffs()[patchi];
                    const scalar pbcf = tpbc()[patchFacei];

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<scalar>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<scalar>::zero;

                    // Check if the live cell is on local or neighbour side
                    if (ibFaceCells[faceI] > -1)
                    {
                        if (ibGamma[ibFaceCells[faceI]] > SMALL)
                        {
                            diag[ibFaceCells[faceI]] += picf;

                            source[ibFaceCells[faceI]] +=
                                 cmptMultiply(pbcf, pnff);
                        }
                    }
                }
            }
        }
    }
}


void immersedBoundaryFvMesh::manipulateMatrix
(
    fvVectorMatrix& eqn
) const
{
    forAll(objectsList(),objectID)
    {
       
       manipulateMatrix(eqn,objectID);
        
    }
}


void immersedBoundaryFvMesh::manipulateMatrix
(
    fvScalarMatrix& eqn
) const
{
    forAll(objectsList(),objectID)
    {
       
       manipulateMatrix(eqn,objectID);
        
    }
}



void immersedBoundaryFvMesh::correctDiag
(
    fvVectorMatrix& eqn
) const
{
    forAll(objectsList(),objectID)
    {
        BCtype_ = objectDictList()[objectID].lookupOrDefault<word>("BCtype","fixedWall");

        if(IBtypeList()[objectID]=="classic")
        {
            correctDiag(eqn,objectID);
        }
        else if(IBtypeList()[objectID]=="ghost-cell")
        {
        }
        else if(IBtypeList()[objectID]=="mix")
        {
            ghostCorrectDiag(eqn,objectID);
        }
    }
}

void immersedBoundaryFvMesh::correctDiag
(
    fvScalarMatrix& eqn
) const
{
    forAll(objectsList(),objectID)
    {
        BCtype_ = objectDictList()[objectID].lookupOrDefault<word>("BCtype","fixedWall");

        if(IBtypeList()[objectID]=="classic")
        {
            correctDiag(eqn,objectID);
        }
        else if(IBtypeList()[objectID]=="ghost-cell")
        {
        }
        else if(IBtypeList()[objectID]=="mix")
        {
            ghostCorrectDiag(eqn,objectID);
        }
    }
}

void immersedBoundaryFvMesh::correctOffDiag
(
    fvScalarMatrix& eqn
) const
{
    forAll(objectsList(),objectID)
    {
        BCtype_ = objectDictList()[objectID].lookupOrDefault<word>("BCtype","fixedWall");

        if(IBtypeList()[objectID]=="classic")
        {
            correctOffDiag(eqn,objectID);
        }
        else if(IBtypeList()[objectID]=="ghost-cell")
        {
        }
        else if(IBtypeList()[objectID]=="mix")
        {
            // not implemented
        }
    }
}


template<class Type>
void Foam::immersedBoundaryFvMesh::fieldExportToMesh
(
    const word& psiName,
    const Field<Type>& psi, //field value
    label objectID
)const
{
    if(IBtypeList()[objectID]=="classic")
    {
        const Field<Type> ibValues(psi,ibCellsList()[objectID]);
        // set only IB cell values to extrueded values
        hitPointExportToMesh(psiName,ibValues,objectID);
    }
    else if(IBtypeList()[objectID]=="mix")
    {
        const Field<Type> ghostValues(psi,ghostCellsList()[objectID]);
        // set only ghost cell values to extrueded values
        hitPointExportToMesh(psiName,ghostValues,objectID);
    }
}


// ************************************************************************* //



