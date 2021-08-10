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
    immersedBoundaryFvMeshTurbulence.C
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvMesh.H"
#include "omegaWallFunctionFvPatchScalarField.H"

#include "turbulentTransportModel.H"
#include "wallFuncModules.H"
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// k-epsilon model correction
void Foam::immersedBoundaryFvMesh::kEpsilonCorrection(const turbulenceModel& turbulence)const
{
    forAll(objectsList(),objectID)
    {
        kEpsilonCorrection(turbulence,objectID);
    }
}
void Foam::immersedBoundaryFvMesh::kEpsilonCorrection
(
    const turbulenceModel& turbulence,
    label objectID
)const
{
    volScalarField& epsilon = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("epsilon"));
    volScalarField& k = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("k"));
    volScalarField::Internal& G = const_cast< volScalarField::Internal&>
        (this->lookupObject< volScalarField::Internal>(turbulence.GName()));

    volVectorField& U = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("U"));
    volScalarField& nut = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("nut"));
    volScalarField& nu = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("nu"));

    scalar Cmu_(0.09);
    scalar kappa_(0.41);
    scalar E_(9.8);
    scalar Cs_(0.5);
    scalar Ks_ = ibProperties().lookupOrDefault<scalar>("roughnessHeight",0);
    Ks_ = objectDictList()[objectID].lookupOrDefault<scalar>("roughnessHeight",Ks_);

    scalar yPlusLam_ = 11.0;

    for (int i=0; i<10; i++)
    {
        yPlusLam_ = log(max(E_*yPlusLam_, 1))/kappa_;
    }

    const scalar Cmu25 = pow(Cmu_, 0.25);
    const scalar Cmu50 = sqrt(Cmu_);
    const scalar Cmu75 = pow(Cmu_, 0.75);

    const pointField& ibHitPoints = ibHitPointsList()[objectID];
    const pointField& samplingPoints = samplingPointsList()[objectID];
    const labelList& ibCells = ibCellsList()[objectID];

    const pointField ibc_centers(this->C(),ibCells);
    vectorField ibNormals = ibHitPoints-samplingPoints;
    ibNormals=ibNormals/mag(ibNormals);

    vectorField USample(samplingPointsValues(U,objectID));

    USample=USample-(USample&ibNormals)*ibNormals;

    scalarField USampleMag=mag(USample)+SMALL;

    // k at sampling point
    scalarField kSample(samplingPointsValues(k,objectID));

    // epsilon at sampling point
    scalarField epsilonSample(SMALL+samplingPointsValues(epsilon,objectID));

    // nu at sampling point
    scalarField nuSample(nu,ibCells);

    // IB distance
    scalarField yIB(mag(ibc_centers-ibHitPoints));

    // sampling distance
    scalarField ySample(mag(samplingPoints-ibHitPoints));

    scalarField gradUSampleMag=mag(USampleMag/(ySample+SMALL));
    Info<<average(USampleMag)<<tab<<average(ySample)<<tab<<average(kSample)<<endl;
    Info<<max(USampleMag)<<tab<<max(ySample)<<tab<<max(kSample)<<endl;

    scalarField kNew(k,ibCells);
    scalarField epsilonNew(epsilon,ibCells);
    scalarField nutNew(nut,ibCells);
    scalarField GNew(G,ibCells);
    scalarField UTanNew(ibCells.size(),pTraits<scalar>::zero);

    // Calculate yPlus for sample points
    scalarField ypd = Cmu25*ySample*sqrt(kSample)/nuSample;
    forAll(ypd, facei)
    {
        scalar kappaRe = kappa_*USampleMag[facei]*ySample[facei]/nuSample[facei];

        scalar yp = yPlusLam_;
        scalar ryPlusLam = 1.0/yp;

        int iter = 0;
        scalar yPlusLast = 0.0;

        do
        {
            yPlusLast = yp;
            yp = (kappaRe + yp)/(1.0 + log(E_*yp));

        } while (mag(ryPlusLam*(yp - yPlusLast)) > 0.01 && ++iter < 10 );

        ypd[facei] = max(0.0, yp);
    }


    scalarField ypIB=ypd;

    vectorField tauWall(ibCells.size(),pTraits<vector>::zero);

    forAll(ibCells, ibCellI)
    {
        const scalar nuLam = nuSample[ibCellI];

        // Calculate yPlus from k and laminar viscosity for the IB point

        const scalar yPlusSample = ypd[ibCellI];
        scalar uTau;

        if (yPlusSample > yPlusLam_)
        {
            // Calculate uTau from log-law, knowing sampled k and U
            uTau = Cmu25*sqrt(kSample[ibCellI]);

        }
        else
        {
            // Sampling point is in laminar sublayer
            // Bug fix: Xu, 21/Aug/2017
            uTau = sqrt(nuLam*USampleMag[ibCellI]/ySample[ibCellI]);
            //uTau = Cmu25*sqrt(kSample[ibCellI]);
        }

        // Calculate KsPlus
        scalar KsPlus = uTau*Ks_/nuLam;

        // Calculate Edash
        scalar Edash = E_;

        if (KsPlus > 2.25)
        {
            scalar fnRough=1.0 + Cs_*KsPlus;
            if(KsPlus < 90.0)
            {
                fnRough = pow
                (
                  (KsPlus - 2.25)/87.75 + Cs_*KsPlus,
                  sin(0.4258*(log(KsPlus) - 0.811))
                );
            }
            Edash /= fnRough;
        }

        // Set wall shear stress
        tauWall[ibCellI] = sqr(uTau)*USample[ibCellI]/(USampleMag[ibCellI] + SMALL);

        // Calculate yPlus for IB point
 
        scalar yPlusIB = yPlusSample*yIB[ibCellI]/(SMALL+ySample[ibCellI]);
        ypIB[ibCellI] = yPlusIB;

        // Calculate wall function data in the immersed boundary point
        if (yPlusIB > yPlusLam_)
        {

            scalar nutwb = max(0,nuLam*(yPlusIB*kappa_/log(Edash*yPlusIB) - 1.0));

            // Fix generation even though it if is not used
            GNew[ibCellI] =
                (nutwb + nuLam)*gradUSampleMag[ibCellI]
                *Cmu25*sqrt(k[ibCellI])
                /(kappa_*yIB[ibCellI]+SMALL);

            // Log-Law for tangential velocity
            UTanNew[ibCellI] =
                min
                (
                    USampleMag[ibCellI],
                    uTau*(1.0/kappa_*log(Edash*(SMALL+yPlusIB)))
                );

            // Calculate turbulent viscosity
            nutNew[ibCellI] = nutwb;


            // Calculate k in the IB cell from G = epsilon
            kNew[ibCellI] = (nutwb + nuLam)*gradUSampleMag[ibCellI]/Cmu50;

            
            epsilonNew[ibCellI] =
                //Cmu75*pow(kNew[ibCellI], 1.5)/(kappa_*yIB_tilda+SMALL);
                Cmu75*pow(kNew[ibCellI], 1.5)/(kappa_*yIB[ibCellI]+SMALL);
        }
           else
        {
            // G is zero
            G[ibCells[ibCellI]] = 0;

            // Laminar sub-layer for tangential velocity: uPlus = yPlus
            UTanNew[ibCellI] = min(USampleMag[ibCellI], uTau*yPlusIB);

          
            // Turbulent viscosity is zero
            nutNew[ibCellI] = SMALL;

            // k is zero gradient: use the sampled value
            kNew[ibCellI] = SMALL+kSample[ibCellI];
            
            epsilonNew[ibCellI] = Cmu_*sqr(kNew[ibCellI])/nuLam;
    
        }
    }
    scalarList UPout(5);
    scalarList nutPout(5);
    scalarList kPout(5);
    scalarList epsilonPout(5);
    scalarList GPout(5);
    scalarList yPout(5);
    scalarList yIBPout(5);
    vectorList tauWallPout(4);

    if (Pstream::parRun() and debug)//Start of mpi run
    {
        UPout[0] = gMin(UTanNew);
        nutPout[0] = gMin(nutNew);
        kPout[0] = gMin(kNew);
        epsilonPout[0] = gMin(epsilonNew);
        GPout[0] = gMin(GNew);
        yPout[0] = gMin(ypd);
        yIBPout[0] = gMin(ypIB);
        tauWallPout[0] = gMin(tauWall);

        UPout[1] = gMax(UTanNew);
        nutPout[1] = gMax(nutNew);
        kPout[1] = gMax(kNew);
        epsilonPout[1] = gMax(epsilonNew);
        GPout[1] = gMax(GNew);
        yPout[1] = gMax(ypd);
        yIBPout[1] = gMax(ypIB);
        tauWallPout[1] = gMax(tauWall);

        UPout[2] = gAverage(UTanNew);
        nutPout[2] = gAverage(nutNew);
        kPout[2] = gAverage(kNew);
        epsilonPout[2] = gAverage(epsilonNew);
        GPout[2] = gAverage(GNew);
        yPout[2] = gAverage(ypd);
        yIBPout[2] = gAverage(ypIB);
        tauWallPout[2] = gAverage(tauWall);
        Info<< "UTangentialNew nutNew kNew epsilonNew GNew yPlus yPlusIB tauWall" << endl;

        for (label I = 0; I < 3; I++)
        {
            if(I==0){Info<<"min ";}
            if(I==1){Info<<"max ";}
            if(I==2){Info<<"ave ";}
            Info<<UPout[I]<<" "
                <<nutPout[I]<<" "
                <<kPout[I]<<" "
                <<epsilonPout[I]<<" "
                <<GPout[I]<<" "
                <<yPout[I]<<" "
                <<yIBPout[I]<<" "
                <<tauWallPout[I]<<" "<<endl;
        }

    }//End of mpi run
    else  if (!Pstream::parRun())//Start of serial run
    {

        UPout[0] = gMin(UTanNew);
        nutPout[0] = gMin(nutNew);
        kPout[0] = gMin(kNew);
        epsilonPout[0] = gMin(epsilonNew);
        GPout[0] = gMin(GNew);
        yPout[0] = gMin(ypd);
        yIBPout[0] = gMin(ypIB);
        tauWallPout[0] = gMin(tauWall);

        UPout[1] = gMax(UTanNew);
        nutPout[1] = gMax(nutNew);
        kPout[1] = gMax(kNew);
        epsilonPout[1] = gMax(epsilonNew);
        GPout[1] = gMax(GNew);
        yPout[1] = gMax(ypd);
        yIBPout[1] = gMax(ypIB);
        tauWallPout[1] = gMax(tauWall);

        UPout[2] = gAverage(UTanNew);
        nutPout[2] = gAverage(nutNew);
        kPout[2] = gAverage(kNew);
        epsilonPout[2] = gAverage(epsilonNew);
        GPout[2] = gAverage(GNew);
        yPout[2] = gAverage(ypd);
        yIBPout[2] = gAverage(ypIB);
        tauWallPout[2] = gAverage(tauWall);
        Info<< "UTangentialNew nutNew kNew epsilonNew GNew yPlus yPlusIB tauWall" << endl;
        for (label I = 0; I < 3; I++)
        {
            if(I==0){Info<<"min ";}
            if(I==1){Info<<"max ";}
            if(I==2){Info<<"ave ";}
            Info<<UPout[I]<<" "
                <<nutPout[I]<<" "
                <<kPout[I]<<" "
                <<epsilonPout[I]<<" "
                <<GPout[I]<<" "
                <<yPout[I]<<" "
                <<yIBPout[I]<<" "
                <<tauWallPout[I]<<" "<<endl;
        }
    }//End of serial ru


    forAll(ibCells,I)
    {
        label cellID=ibCells[I];
        k[cellID]=kNew[I];
        epsilon[cellID]=epsilonNew[I];
        nut[cellID]=nutNew[I];
        G[cellID]=GNew[I];
        vector Utan =(tensor::I-sqr(ibNormals[I])) & U[cellID];
        Utan /=mag(Utan)+VSMALL;
        if(IBtypeList()[objectID]!="mix")
        {
            U[cellID]=Utan*UTanNew[I]+(U[cellID]&ibNormals[I])*ibNormals[I];
        }
    }


    wallShearStressListPtr_->set
    (
        objectID,
        tauWall
    );


    hitPointExportToMesh("ypIB",ypIB,objectID);
}

// k-omega model correction
void Foam::immersedBoundaryFvMesh::kOmegaCorrection(const turbulenceModel& turbulence)const
{
    forAll(objectsList(),objectID)
    {
        if(IBtypeList()[objectID]=="classic")
        {
            kOmegaIbCorrection(turbulence,objectID);
            Info << "using kOmegaIbCorrection  "<<endl;

        }
        
    }
}




// k-omega model correction
void Foam::immersedBoundaryFvMesh::nutCorrection()const
{
    forAll(objectsList(),objectID)
    {
        if(IBtypeList()[objectID]=="classic")
        {
        
            nutIbCorrection(objectID);
            Info << "using nutIbCorrection  "<<endl;
             
        }
        
    }
}




void Foam::immersedBoundaryFvMesh::nutIbCorrection( label objectID )const
{
    const labelList& ibCells = ibCellsList()[objectID];
    const scalarField& ibNut = nutNew(objectID);
    volScalarField& nut = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("nut"));
    
    forAll(ibCells, ibCellI)
    {
        label cellI = ibCells[ibCellI];  
        nut[cellI]=ibNut[ibCellI];
        
    }


}





//  wall function based on Rey Deleon, Inanc Senocak 2018: Simulation of turbulence flow over
// complex terrain using an immersed boundary method


void Foam::immersedBoundaryFvMesh::kOmegaIbCorrection
(
    const turbulenceModel& turbulence,
    label objectID
)const
{
    volScalarField& omega = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("omega"));
    volScalarField& k = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("k"));
    volScalarField::Internal& G = const_cast< volScalarField::Internal&>
        (this->lookupObject< volScalarField::Internal>(turbulence.GName()));
    volVectorField& U = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("U"));
    volScalarField& nut = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("nut"));
    volScalarField& nu = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("nu"));
    volVectorField& shearVelocity = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("shearVelocity"));
    volScalarField& yplusIB = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("yplusIB"));
    volScalarField& yplus = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("yplus"));
    scalar Cmu_(0.09);
    scalar kappa_(0.41);

    
    
    scalar Ks_ = ibProperties().lookupOrDefault<scalar>("roughnessHeight",0);
    Ks_ = objectDictList()[objectID].lookupOrDefault<scalar>("roughnessHeight",Ks_);
    
    scalar roughnessFactor_ = ibProperties().lookupOrDefault<scalar>("roughnessFactor",1);
    roughnessFactor_ = objectDictList()[objectID].lookupOrDefault<scalar>("roughnessFactor",roughnessFactor_);


    scalar roughnessConstant_ = ibProperties().lookupOrDefault<scalar>("roughnessConstant",0.5);
    roughnessConstant_ = objectDictList()[objectID].lookupOrDefault<scalar>("roughnessConstant_",roughnessConstant_);

  

    const scalar Cmu25 = pow(Cmu_, 0.25);

    const pointField& ibHitPoints = ibHitPointsList()[objectID];
    const pointField& samplingPoints = samplingPointsList()[objectID];

    const labelList& ibCells = ibCellsList()[objectID];

    const pointField ibc_centers(this->C(),ibCells);
    vectorField ibNormals = ibHitPoints-samplingPoints;
    ibNormals=ibNormals/mag(ibNormals);
    

    
    vectorField USample(samplingPointsValues(U,objectID));


    // remove normal velocity
    vectorField UTanSample=USample-(USample&ibNormals)*ibNormals;

    scalarField UTanSampleMag=mag(UTanSample)+SMALL;
      
    // k at sampling point
    scalarField kSample(samplingPointsValues(k,objectID)+SMALL);

    // epsilon at sampling point
    scalarField omegaSample(SMALL+samplingPointsValues(omega,objectID)+SMALL);

    // nu at sampling point
    scalarField nuSample(nu,ibCells);

    // nut at sampling point
    scalarField nutSample(samplingPointsValues(nut,objectID)+SMALL);
    
    // IB distance
    scalarField yIB(mag(ibc_centers-ibHitPoints)+SMALL);

    // sampling distance
    scalarField ySample(mag(samplingPoints-ibHitPoints)+SMALL);

    scalarField gradUSampleMag=mag(UTanSampleMag/(ySample+SMALL));
    if(debug)
    {
        Info<<"Interpolted values at sampled points"<<endl;
        Info<<"    "<<"U"<<tab<<"k"<<tab<<"ySample"<<tab<<"yIB"<<endl;
        Info<<"MIN "<<gMin(UTanSampleMag)<<tab<<gMin(kSample)<<tab<<gMin(ySample)<<tab<<gMin(yIB)<<endl;
        Info<<"AVE "<<gAverage(UTanSampleMag)<<tab<<tab<<gAverage(kSample)<<tab<<gAverage(ySample)<<tab<<gAverage(yIB)<<endl;
        Info<<"MAX "<<gMax(UTanSampleMag)<<tab<<tab<<gMax(kSample)<<tab<<gMax(ySample)<<tab<<gMax(yIB)<<endl;
    }
    scalarField kNew(k,ibCells);
    scalarField omegaNew(omega,ibCells);
    scalarField nutNew(nut,ibCells);
    scalarField utauNew(ibCells.size(),0);
    scalarField GNew(G,ibCells);
    scalarField UTanNew(ibCells.size(),pTraits<scalar>::zero);
    scalarField UTanOld(ibCells.size(),pTraits<scalar>::zero);
    // Calculate yPlus for sample points
    scalarField ypd = Cmu25*ySample*sqrt(kSample)/nuSample;
    scalarField ypIB=ypd;



    vectorField tauWall(ibCells.size(),pTraits<vector>::zero);
    
    word wallFunc_ = ibProperties().lookupOrDefault<word>("wallFunction","log law");

    wallFunc_ = objectDictList()[objectID].lookupOrDefault<word>("wallFunction",wallFunc_);
    
    // wallFuncModules is to perform wall functions related calculations.
    // 3 different velocity profiles are implemented: 
    // two-layer, vanDriest, Spalding, (with roughness)
    // It has multiple functions:
    // 1. uTau -- calcualte shear velocity (with roughness) based on U, y
    // 2. uPlus
    // 3. kOmegaUpdate -- include Tamaki2017
    // 4. dudy, dupdyp -- derivatives
    wallFuncModules wFM(wallFunc_,"Tamaki2017");
    
    
    forAll(ibCells, ibCellI)
    { 
        const scalar nuLam = nuSample[ibCellI];
 
        scalar& UT =  UTanNew[ibCellI];
     
        scalar& yPlusSample = ypd[ibCellI];

     	yPlusSample = wFM.yPlus
         (
	        UTanSampleMag[ibCellI],
	        ySample[ibCellI],
	        Ks_,
	        roughnessFactor_,
	        roughnessConstant_
         );
         
        scalar& yPlusIB = ypIB[ibCellI];
        yPlusIB = yPlusSample*yIB[ibCellI]/(SMALL+ySample[ibCellI]);
        scalar uTau = yPlusSample*nuLam/(ySample[ibCellI]+SMALL);
     
        scalar kPlus = 1.0/0.3;
        kNew[ibCellI] = kPlus*sqr(uTau); 
         
               
        scalar  omega_vis = 6.0*nuLam/0.075/sqr(ySample[ibCellI]+SMALL);
        scalar  omega_log = sqrt(kNew[ibCellI])/(Cmu25*kappa_*ySample[ibCellI]+SMALL);
        omegaNew[ibCellI] = sqrt(sqr(omega_vis)+sqr(omega_log));  

        scalar nutwb = kNew[ibCellI]/(omegaNew[ibCellI] +SMALL);
        
        scalar dupdyp = nuLam/(nuLam+nutwb);
        
        scalar dudy = sqr(uTau)/(nuLam+nutwb);
        
        UT = max(0.0,(UTanSampleMag[ibCellI] - dupdyp*uTau*(yPlusSample-yPlusIB)));

             
        GNew[ibCellI] =                 
            (nutwb + nuLam)*(dudy)
            *Cmu25*sqrt(kNew[ibCellI])
            /(kappa_*ySample[ibCellI]+SMALL);

       
        nutNew[ibCellI] = nutwb;
        
        utauNew[ibCellI] = uTau;
       
    
        tauWall[ibCellI] = sqr(uTau)*UTanSample[ibCellI]/(UTanSampleMag[ibCellI] + SMALL);
       
            
    }
    // calculate weights for boundary patch
    // in case that one cell is shared by IB cell and Boundary cell
    
        
    const volScalarField::Boundary& bf = omega.boundaryField();

    scalarField weights(omega.primitiveField().size(),0.0);
    
    forAll(bf, patchi)
    {
        if (isA<omegaWallFunctionFvPatchScalarField>(bf[patchi]))// or bf[patchi].patch().name()=="inlet")
        {
            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i)
            {
                label celli = faceCells[i];
                weights[celli]++;
            }
        }
    }

    forAll(ibCells,I)
    {
        scalar cellID=ibCells[I];
        weights[cellID]++;
        scalar w = 1.0/(SMALL+weights[cellID]);
        omegaNew[I] = omegaNew[I]*w+omega[cellID]*(1.0-w);
        GNew[I] = GNew[I]*w+G[cellID]*(1.0-w);
        nutNew[I] = nutNew[I]*w+nut[cellID]*(1.0-w);
        kNew[I] = kNew[I]*w+k[cellID]*(1.0-w);
        vector tmpUtan = U[cellID]-(U[cellID]&ibNormals[I])*ibNormals[I];

        if(IBtypeList()[objectID]=="classic")
        {
            UTanNew[I] = UTanNew[I]*w+mag(tmpUtan)*(1.0-w);
        }
    }
    
    
    
    

    scalarList UPout(5);
    scalarList nutPout(5);
    scalarList utauPout(5);
    scalarList kPout(5);
    scalarList omegaPout(5);
    scalarList GPout(5);
    scalarList yPout(5);
    scalarList yIBPout(5);
    vectorList tauWallPout(4);
    
    UPout[0] = gMin(UTanNew);
    
    nutPout[0] = gMin(nutNew);
    utauPout[0] = gMin(utauNew);
    kPout[0] = gMin(kNew);
    omegaPout[0] = gMin(omegaNew);
    GPout[0] = gMin(GNew);
    yPout[0] = gMin(ypd);
    yIBPout[0] = gMin(ypIB);
    tauWallPout[0] = gMin(tauWall);

    UPout[1] = gMax(UTanNew);
    nutPout[1] = gMax(nutNew);
    utauPout[1] = gMax(utauNew);
    kPout[1] = gMax(kNew);
    omegaPout[1] = gMax(omegaNew);
    GPout[1] = gMax(GNew);
    yPout[1] = gMax(ypd);
    yIBPout[1] = gMax(ypIB);
    tauWallPout[1] = gMax(tauWall);


    UPout[2] = gAverage(UTanNew);
    nutPout[2] = gAverage(nutNew);
    utauPout[2] = gAverage(utauNew);
    kPout[2] = gAverage(kNew);
    omegaPout[2] = gAverage(omegaNew);
    GPout[2] = gAverage(GNew);
    yPout[2] = gAverage(ypd);
    yIBPout[2] = gAverage(ypIB);
    tauWallPout[2] = gAverage(tauWall);

    Info<< "UTangentialNew nutNew utauNew kNew omegaNew GNew yPlus yPlusIB tauWall" << endl;

    for (label I = 0; I < 3; I++)
    {
        if(I==0){Info<<"min ";}
        if(I==1){Info<<"max ";}
        if(I==2){Info<<"ave ";}
        Info<<UPout[I]<<" "
            <<nutPout[I]<<" "
            <<utauPout[I]<<" "
            <<kPout[I]<<" "
            <<omegaPout[I]<<" "
            <<GPout[I]<<" "
            <<yPout[I]<<" "
            <<yIBPout[I]<<" "
            <<tauWallPout[I]<<" "<<endl;
    }
    volVectorField& U_desired = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("U_desired"));

    U_desired = U;
     // insert IB values
    forAll(ibCells,I)
    {
        label cellID = ibCells[I];
        omega[cellID] = omegaNew[I];
        nut[cellID] = nutNew[I];
       G[cellID] = GNew[I];
       shearVelocity[cellID] = utauNew[I]*UTanSample[I]/(UTanSampleMag[I] + SMALL);
       yplusIB[cellID] = ypIB[I];
       yplus[cellID] =ypd[I];
        vector Unormal = (U[cellID]&ibNormals[I])*ibNormals[I];

        vector Utan = UTanSample[I]/(UTanSampleMag[I]+SMALL);
        k[cellID] = kNew[I];
        Unormal = Unormal*(sqr(yIB[I])/sqr(ySample[I]));// Unormal should be removed or scaled?
        if(IBtypeList()[objectID]!="mix")
        {
            U[cellID] = Utan*UTanNew[I]+Unormal;
        }
    }

    // transfer wall shear stress to pointer
    wallShearStressListPtr_->set
    (
        objectID,
        tauWall
    );
    // transfer nut to pointer
    nutListPtr_->set
    (
        objectID,
        nutNew
    );
    
    if(debug)
    {
        hitPointExportToMesh("yPlus",ypd,objectID);
        hitPointExportToMesh("kNew",kNew,objectID);
        hitPointExportToMesh("omegaNew",omegaNew,objectID);
        hitPointExportToMesh("nutNew",nutNew,objectID);
   
	    hitPointExportToMesh("ypIB",ypIB,objectID);
        hitPointExportToMesh("yIB",yIB,objectID);
        hitPointExportToMesh("ySample",ySample,objectID);
        hitPointExportToMesh("UTanNew",UTanNew,objectID);
        hitPointExportToMesh("USampleMag",UTanSampleMag,objectID);
    }
    
    const dictionary& dict = objectDictList()[objectID].subDict("sediment");
    bool changeSTL(dict.lookupOrDefault<bool>("changeSTL",true));
    if(!changeSTL)
    {
        hitPointExportToMesh("wallShearStress",tauWall,objectID);
    }

    
}




template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::samplingPointsValues
(
    const Field<Type>& psi,
    label objectID
)const
{
    Field<Type> ibU
    (
        stencilInterpolation
        (
            psi,
            objectID,
            ibCellsList()[objectID],
            samplingStencils(objectID).cellCells(),
            samplingStencils(objectID).cellProcCells(),
            samplingStencils(objectID).procCells(),
            samplingStencils(objectID).procCentres(),
            samplingStencils(objectID).weights(),
            samplingStencils(objectID).procWeights()
        )
    );
    return ibU;
}

template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::farSamplingPointsValues
(   
    const Field<Type>& psi,
    label objectID
)const
{
    Field<Type> ibU
    (
        stencilInterpolation
        (
            psi,
            objectID,
            ibCellsList()[objectID],
            farSamplingStencils(objectID).cellCells(),
            farSamplingStencils(objectID).cellProcCells(),
            farSamplingStencils(objectID).procCells(),
            farSamplingStencils(objectID).procCentres(),
            farSamplingStencils(objectID).weights(),
            farSamplingStencils(objectID).procWeights()
        )
    );
    return ibU;
}





template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::ibFaceSamplingPointsValues
(
    const Field<Type>& psi,
    label objectID
)const
{
    Field<Type> ibU
    (
        stencilInterpolation
        (
            psi,
            objectID,
            ibFacesList()[objectID],
            ibFaceSamplingStencils(objectID).cellCells(),
            ibFaceSamplingStencils(objectID).cellProcCells(),
            ibFaceSamplingStencils(objectID).procCells(),
            ibFaceSamplingStencils(objectID).procCentres(),
            ibFaceSamplingStencils(objectID).weights(),
            ibFaceSamplingStencils(objectID).procWeights()
        )
    );
    return ibU;
}

template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::ibCellCentersValues
(
    const Field<Type>& psi,
    label objectID
)const
{
    Field<Type> ibU
    (
        stencilInterpolation
        (
            psi,
            objectID,
            ibCellsList()[objectID],
            ibCellStencils(objectID).cellCells(),
            ibCellStencils(objectID).cellProcCells(),
            ibCellStencils(objectID).procCells(),
            ibCellStencils(objectID).procCentres(),
            ibCellStencils(objectID).weights(),
            ibCellStencils(objectID).procWeights()
        )
    );
    return ibU;
}

template<class Type>
Foam::Field<Type> Foam::immersedBoundaryFvMesh::imagePointsValues
(
    const Field<Type>& psi,
    label objectID
)const
{
    Field<Type> ibU
    (
        stencilInterpolation
        (
            psi,
            objectID,
            ghostCellsList()[objectID],
            imageStencils(objectID).cellCells(),
            imageStencils(objectID).cellProcCells(),
            imageStencils(objectID).procCells(),
            imageStencils(objectID).procCentres(),
            imageStencils(objectID).weights(),
            imageStencils(objectID).procWeights()
        )
    );
    return ibU;
}

template<class Type>
void Foam::immersedBoundaryFvMesh::hitPointExportToMesh
(
    const word& psiName,
    const Field<Type>& psi, //hit pt values it can be ib hit or ghost hit
    label objectID,
    bool write //default=true
)const
{
    
    if(dual(objectID))
    {
        GeometricField< Type, fvPatchField, volMesh > PSI
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

        hitPointExportToMesh(psi,PSI,objectID,write,true);
    }
}

template<class Type>
void Foam::immersedBoundaryFvMesh::hitPointExportToMesh
(
    const Field<Type>& psi,
    GeometricField<Type, fvPatchField, volMesh>& PSI,
    label objectID,
    bool write, //default=true
    bool dual
)const
{
    if(dual)
    {
        
        if(IBtypeList()[objectID]=="classic")
        {
            setValueToDualMesh
            (
                objectID,
                PSI,// dualMesh
                mapFromIBHitToDualFace(psi,objectID)// dual patch value
            );
        }
        else if (IBtypeList()[objectID]=="mix")
        {
            setValueToDualMesh
            (
                objectID,
                PSI,// dualMesh
                mapFromGHOSTHitToDualFace(psi,objectID)// dual patch value
            );
        }
    }
    
    scalar timeValue=time().timeOutputValue()+1e-6;
    scalar timeValue0=time().timeOutputValue()-time().deltaT0Value()+1e-6;

    scalar extrudeMeshOutputTimeStep= ibProperties().lookupOrDefault<scalar>("extrudeMeshOutputTimeStep",-1);

    // write out
    if (time().outputTime() or floor(timeValue/extrudeMeshOutputTimeStep)>floor(timeValue0/extrudeMeshOutputTimeStep))
    {
         if(write and Pstream::myProcNo()==Pstream::masterNo()) parallelWrite(PSI);
    }
}






/*
//dynamicKEqnIBCorrection
void Foam::immersedBoundaryFvMesh::dynamicKEqnIBCorrection
(
    const turbulenceModel& turbulence
)const
{
    forAll(this->objectsList(),objectID)
    {
        if(IBtypeList()[objectID]=="classic")
        {
               dynamicKEqnIBCorrection(turbulence,objectID);
        }

    }

}






void Foam::immersedBoundaryFvMesh::dynamicKEqnIBCorrection
(
    const turbulenceModel& turbulence,
    label objectID
)const
{


    volScalarField& k = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("k"));
    volVectorField& U = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("U"));
    volScalarField& nu = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("nu"));
    volVectorField& shearVelocity = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("shearVelocity"));
    volScalarField& yplusIB = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("yplusIB"));
    volScalarField& yplus = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("yplus"));
    scalar Cmu_(0.09);
    scalar kappa_(0.41);

    
    
    scalar Ks_ = ibProperties().lookupOrDefault<scalar>("roughnessHeight",0);
    Ks_ = objectDictList()[objectID].lookupOrDefault<scalar>("roughnessHeight",Ks_);
    
    scalar roughnessFactor_ = ibProperties().lookupOrDefault<scalar>("roughnessFactor",1);
    roughnessFactor_ = objectDictList()[objectID].lookupOrDefault<scalar>("roughnessFactor",roughnessFactor_);


    scalar roughnessConstant_ = ibProperties().lookupOrDefault<scalar>("roughnessConstant",0.5);
    roughnessConstant_ = objectDictList()[objectID].lookupOrDefault<scalar>("roughnessConstant_",roughnessConstant_);

  

    const scalar Cmu25 = pow(Cmu_, 0.25);

    const pointField& ibHitPoints = ibHitPointsList()[objectID];
    const pointField& samplingPoints = samplingPointsList()[objectID];

    const labelList& ibCells = ibCellsList()[objectID];

    const pointField ibc_centers(this->C(),ibCells);
    vectorField ibNormals = ibHitPoints-samplingPoints;
    ibNormals=ibNormals/mag(ibNormals);
    

    
    vectorField USample(samplingPointsValues(U,objectID));


    // remove normal velocity
    vectorField UTanSample=USample-(USample&ibNormals)*ibNormals;

    scalarField UTanSampleMag=mag(UTanSample)+SMALL;
    // nu at sampling point
    scalarField nuSample(nu,ibCells);
      
    // k at sampling point
    scalarField kSample(samplingPointsValues(k,objectID)+SMALL);

    // IB distance
    scalarField yIB(mag(ibc_centers-ibHitPoints)+SMALL);

    // sampling distance
    scalarField ySample(mag(samplingPoints-ibHitPoints)+SMALL);


    if(debug)
    {
        Info<<"Interpolted values at sampled points"<<endl;
        Info<<"    "<<"U"<<tab<<"k"<<tab<<"ySample"<<tab<<"yIB"<<endl;
        Info<<"MIN "<<gMin(UTanSampleMag)<<tab<<gMin(kSample)<<tab<<gMin(ySample)<<tab<<gMin(yIB)<<endl;
        Info<<"AVE "<<gAverage(UTanSampleMag)<<tab<<tab<<gAverage(kSample)<<tab<<gAverage(ySample)<<tab<<gAverage(yIB)<<endl;
        Info<<"MAX "<<gMax(UTanSampleMag)<<tab<<tab<<gMax(kSample)<<tab<<gMax(ySample)<<tab<<gMax(yIB)<<endl;
    }
    scalarField kNew(k,ibCells);
    scalarField utauNew(ibCells.size(),0);
    scalarField UTanNew(ibCells.size(),pTraits<scalar>::zero);
    // Calculate yPlus for sample points
    scalarField ypd = Cmu25*ySample*sqrt(kSample)/nuSample;
    scalarField ypIB=ypd;



    vectorField tauWall(ibCells.size(),pTraits<vector>::zero);
    
    word wallFunc_ = ibProperties().lookupOrDefault<word>("wallFunction","vanDriest");

    wallFunc_ = objectDictList()[objectID].lookupOrDefault<word>("wallFunction",wallFunc_);
    
    // wallFuncModules is to perform wall functions related calculations.
    wallFuncModules wFM(wallFunc_,"Baykal2015");
    
    
    forAll(ibCells, ibCellI)
    { 
        const scalar nuLam = nuSample[ibCellI];
 
        scalar& UT =  UTanNew[ibCellI];
     
        scalar& yPlusSample = ypd[ibCellI];

     	yPlusSample = wFM.yPlus
         (
	        UTanSampleMag[ibCellI],
	        ySample[ibCellI],
	        Ks_,
	        roughnessFactor_,
	        roughnessConstant_
         );
         
        scalar& yPlusIB = ypIB[ibCellI];
        yPlusIB = yPlusSample*yIB[ibCellI]/(SMALL+ySample[ibCellI]);
        scalar uTau = yPlusSample*nuLam/(ySample[ibCellI]+SMALL);
     
        scalar kPlus = 1.0/0.3;
        kNew[ibCellI] = kPlus*sqr(uTau);         
        
        scalar dupdyp = 1/(kappa_* yPlusSample + +SMALL);
        
 
        
        UT = max(0.0,(UTanSampleMag[ibCellI] - dupdyp*uTau*(yPlusSample-yPlusIB)));
        
        utauNew[ibCellI] = uTau;
       
    
        tauWall[ibCellI] = sqr(uTau)*UTanSample[ibCellI]/(UTanSampleMag[ibCellI] + SMALL);
       
            
    }

    forAll(ibCells,I)
    {
        scalar cellID=ibCells[I];


        kNew[I] = k[cellID];
        vector tmpUtan = U[cellID]-(U[cellID]&ibNormals[I])*ibNormals[I];

        if(IBtypeList()[objectID]=="classic")
        {
            UTanNew[I] = mag(tmpUtan);
        }
    }
    
    
    
    

    scalarList UPout(5);
    scalarList utauPout(5);
    scalarList kPout(5);
    scalarList yPout(5);
    scalarList yIBPout(5);
    vectorList tauWallPout(4);
    
    UPout[0] = gMin(UTanNew);
    utauPout[0] = gMin(utauNew);
    kPout[0] = gMin(kNew);
    yPout[0] = gMin(ypd);
    yIBPout[0] = gMin(ypIB);
    tauWallPout[0] = gMin(tauWall);

    UPout[1] = gMax(UTanNew);
    utauPout[1] = gMax(utauNew);
    kPout[1] = gMax(kNew);
    yPout[1] = gMax(ypd);
    yIBPout[1] = gMax(ypIB);
    tauWallPout[1] = gMax(tauWall);


    UPout[2] = gAverage(UTanNew);
    utauPout[2] = gAverage(utauNew);
    kPout[2] = gAverage(kNew);
    yPout[2] = gAverage(ypd);
    yIBPout[2] = gAverage(ypIB);
    tauWallPout[2] = gAverage(tauWall);

    Info<< "UTangentialNew utauNew kNew yPlus yPlusIB tauWall" << endl;

    for (label I = 0; I < 3; I++)
    {
        if(I==0){Info<<"min ";}
        if(I==1){Info<<"max ";}
        if(I==2){Info<<"ave ";}
        Info<<UPout[I]<<" "
            <<utauPout[I]<<" "
            <<kPout[I]<<" "
            <<yPout[I]<<" "
            <<yIBPout[I]<<" "
            <<tauWallPout[I]<<" "<<endl;
    }
    volVectorField& U_desired = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("U_desired"));

    U_desired = U;
     // insert IB values
    forAll(ibCells,I)
    {
        label cellID = ibCells[I];
       shearVelocity[cellID] = utauNew[I]*UTanSample[I]/(UTanSampleMag[I] + SMALL);
       yplusIB[cellID] = ypIB[I];
       yplus[cellID] =ypd[I];
        vector Unormal = (U[cellID]&ibNormals[I])*ibNormals[I];
        vector Utan = UTanSample[I]/(UTanSampleMag[I]+SMALL);
        k[cellID] = kNew[I];
        Unormal = Unormal*(sqr(yIB[I])/sqr(ySample[I]));// Unormal should be removed or scaled?
        if(IBtypeList()[objectID]!="mix")
        {
            U[cellID] = Utan*UTanNew[I]+Unormal;
        }
    }

    // transfer wall shear stress to pointer
    wallShearStressListPtr_->set
    (
        objectID,
        tauWall
    );

    
    if(debug)
    {
        hitPointExportToMesh("yPlus",ypd,objectID);
        hitPointExportToMesh("kNew",kNew,objectID);   
	    hitPointExportToMesh("ypIB",ypIB,objectID);
        hitPointExportToMesh("yIB",yIB,objectID);
        hitPointExportToMesh("ySample",ySample,objectID);
        hitPointExportToMesh("UTanNew",UTanNew,objectID);
        hitPointExportToMesh("USampleMag",UTanSampleMag,objectID);
    }
    
    const dictionary& dict = objectDictList()[objectID].subDict("sediment");
    bool changeSTL(dict.lookupOrDefault<bool>("changeSTL",true));
    if(!changeSTL)
    {
        hitPointExportToMesh("wallShearStress",tauWall,objectID);
    }


}






















*/








//SmagorinskyCorrection
void Foam::immersedBoundaryFvMesh::SmagorinskyCorrection
(
    const turbulenceModel& turbulence,
    const volScalarField& k
)const
{
    forAll(this->objectsList(),objectID)
    {
        if(IBtypeList()[objectID]=="classic")
        {
               SmagorinskyCorrection(turbulence,k,objectID);
        }
        else if(IBtypeList()[objectID]=="mix")
        {
            SmagorinskyGhostCorrection(turbulence,k,objectID);
        }
    }

}






void Foam::immersedBoundaryFvMesh::SmagorinskyCorrection
(
    const turbulenceModel& turbulence,
    const volScalarField& k,
    label objectID
)const
{

    volVectorField& U = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("U"));
    volScalarField& nut = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("nut"));
    volScalarField& nu = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("nu"));

    volVectorField Umean (U);
    if(this->foundObject<volVectorField>("UMean"))
    {
        Umean = this->lookupObject<volVectorField>("UMean");
    }
    volVectorField UPrime2Mean (U);
    volScalarField kFull(k);

    if(this->foundObject<volVectorField>("UPrime2Mean"))
    {
        const volVectorField& UPrime2Mean = this->lookupObject<volVectorField>("UPrime2Mean");
        forAll(kFull,I)
        {
            scalar k2Part = std::abs(UPrime2Mean[I].x()+UPrime2Mean[I].y()+UPrime2Mean[I].z());
            kFull[I]=kFull[I]+Foam::sqrt(k2Part);
        }
    }
    scalar Cmu_(0.09);
    scalar kappa_(0.41);
    scalar E_(9.8);
    //scalar Cs_(0.5);
    scalar Ks_ = ibProperties().lookupOrDefault<scalar>("roughnessHeight",0);
    Ks_ = objectDictList()[objectID].lookupOrDefault<scalar>("roughnessHeight",Ks_);

    scalar yPlusLam_ = 11.0;

    for (int i=0; i<10; i++)
    {
        yPlusLam_ = log(max(E_*yPlusLam_, 1))/kappa_;
    }

    const scalar Cmu25 = pow(Cmu_, 0.25);
//    const scalar Cmu50 = sqrt(Cmu_);
//    const scalar Cmu75 = pow(Cmu_, 0.75);

    const pointField& ibHitPoints = ibHitPointsList()[objectID];
    const pointField& samplingPoints = samplingPointsList()[objectID];
    const labelList& ibCells = ibCellsList()[objectID];

    const pointField ibc_centers(this->C(),ibCells);
    vectorField ibNormals = ibHitPoints-samplingPoints;
    ibNormals=ibNormals/mag(ibNormals);

    //full scale U sample
    vectorField USample(samplingPointsValues(Umean,objectID));

    // remove normal velocity
    vectorField tmpUsample=USample-(USample&ibNormals)*ibNormals;
    USample = tmpUsample/(SMALL+mag(tmpUsample))*(mag(USample)+SMALL);
    tmpUsample.clear();
    scalarField USampleMag=mag(USample)+SMALL;

    // nu at sampling point
    scalarField nuSample(nu,ibCells);

    // k at sampling point
    scalarField kSample(samplingPointsValues(k,objectID)+SMALL);

    // nut at sampling point
    scalarField nutSample(samplingPointsValues(nut,objectID)+SMALL);

    // IB distance
    scalarField yIB(mag(ibc_centers-ibHitPoints));

    // sampling distance
    scalarField ySample(mag(samplingPoints-ibHitPoints));

    scalarField gradUSampleMag=mag(USampleMag/(ySample+SMALL));

    scalarField nutNew(nutSample);

    scalarField kNew(kSample);

    scalarField UTanNew(ibCells.size(),pTraits<scalar>::zero);

    // Calculate yPlus for sample points
    scalarField ypd(yIB.size(),0.0);
    scalarField ypIB(yIB.size(),0.0);

    vectorField tauWall(ibCells.size(),pTraits<vector>::zero);

    forAll(ibCells, ibCellI)
    {
        const scalar nuLam = nuSample[ibCellI];

        // eq2.10 from Baykal 2015
        scalar Uf=Cmu25*sqrt(kSample[ibCellI]);
        scalar oldUf=Uf;
        int iter = 0;
        do
        {
            scalar ycp=ySample[ibCellI]*Uf/nuLam;
            scalar RHS=0.0;
            scalar NN=60;
            label ii=0;
            do
            {
                scalar dyp=ycp/NN;
                scalar yp=dyp*(ii+0.5);
                scalar ksp=Ks_*Uf/nuLam;
                scalar dyccp=0.9*sqrt(ksp)-ksp*exp(-ksp/6.0);
                scalar CC=1-exp(-(dyccp+yp)/25.0);
                CC *=CC;
                scalar ff=1+pow((1+4*pow(kappa_*(yp+dyccp),2)*CC),0.5);
                RHS +=dyp*2/ff;
            } while(++ii <NN);
            Uf=sqrt(USampleMag[ibCellI]/RHS*Uf);
        } while (mag(oldUf - Uf) > 1e-7 && ++iter < 20);

        scalar uTau=Uf;

        // Calculate yPlus from k and laminar viscosity for the IB point

        scalar& yPlusSample = ypd[ibCellI];

        yPlusSample=ySample[ibCellI]*uTau/nuLam;

        // Calculate yPlus for IB point
        scalar yPlusIB = yPlusSample*yIB[ibCellI]/(SMALL+ySample[ibCellI]);
        ypIB[ibCellI] = yPlusIB;

        scalar& UT =  UTanNew[ibCellI];
        scalar ycp=yIB[ibCellI]*Uf/nuLam;
        scalar RHS=0.0;
        scalar NN=60;
        label ii=0;
        do
        {
            scalar dyp=ycp/NN;
            scalar yp=dyp*(ii+0.5);
            scalar ksp=Ks_*Uf/nuLam;
            scalar dyccp=0.9*sqrt(ksp)-ksp*exp(-ksp/6.0);
            scalar CC=1-exp(-(dyccp+yp)/25.0);
            CC *=CC;
            scalar ff=1+pow((1+4*pow(kappa_*(yp+dyccp),2)*CC),0.5);
            RHS +=dyp*2/ff;
        } while(++ii <NN);

        UT=RHS*Uf;

        // Set wall shear stress
        tauWall[ibCellI] = sqr(uTau)*USample[ibCellI]/(USampleMag[ibCellI] + SMALL);

        nutNew[ibCellI] =max(0.0, sqr(uTau)/(UT/yIB[ibCellI]));

       }
    scalarList UPout(5);
    scalarList nutPout(5);
    scalarList kPout(5);
    scalarList epsilonPout(5);
    scalarList yPout(5);
    scalarList yIBPout(5);
    vectorList tauWallPout(4);

    if (Pstream::parRun() and debug)//Start of mpi run
    {
        UPout[0] = gMin(UTanNew);
        nutPout[0] = gMin(nutNew);
        kPout[0] = gMin(kNew);
        yPout[0] = gMin(ypd);
        yIBPout[0] = gMin(ypIB);
        tauWallPout[0] = gMin(tauWall);

        UPout[1] = gMax(UTanNew);
        nutPout[1] = gMax(nutNew);
        kPout[1] = gMax(kNew);
        yPout[1] = gMax(ypd);
        yIBPout[1] = gMax(ypIB);
        tauWallPout[1] = gMax(tauWall);

        UPout[2] = gAverage(UTanNew);
        nutPout[2] = gAverage(nutNew);

        yPout[2] = gAverage(ypd);
        yIBPout[2] = gAverage(ypIB);
        tauWallPout[2] = gAverage(tauWall);
        Info<< "subGrid- UTangentialNew nutNew kNew yPlus yPlusIB tauWall" << endl;

        for (label I = 0; I < 3; I++)
        {
            if(I==0){Info<<"min ";}
            if(I==1){Info<<"max ";}
            if(I==2){Info<<"ave ";}
            Info<<UPout[I]<<" "
                <<nutPout[I]<<" "
                <<kPout[I]<<" "
                <<yPout[I]<<" "
                <<yIBPout[I]<<" "
                <<tauWallPout[I]<<" "<<endl;
        }

    }//End of mpi run
    else  if (!Pstream::parRun())//Start of serial run
    {

        UPout[0] = gMin(UTanNew);
        nutPout[0] = gMin(nutNew);
        kPout[0] = gMin(kNew);
        yPout[0] = gMin(ypd);
        yIBPout[0] = gMin(ypIB);
        tauWallPout[0] = gMin(tauWall);

        UPout[1] = gMax(UTanNew);
        nutPout[1] = gMax(nutNew);
        kPout[1] = gMax(kNew);
        yPout[1] = gMax(ypd);
        yIBPout[1] = gMax(ypIB);
        tauWallPout[1] = gMax(tauWall);

        UPout[2] = gAverage(UTanNew);
        nutPout[2] = gAverage(nutNew);
        kPout[2] = gAverage(kNew);
        yPout[2] = gAverage(ypd);
        yIBPout[2] = gAverage(ypIB);
        tauWallPout[2] = gAverage(tauWall);
        Info<< "subGrid- UTangentialNew nutNew kNew yPlus yPlusIB tauWall" << endl;
        for (label I = 0; I < 3; I++)
        {
            if(I==0){Info<<"min ";}
            if(I==1){Info<<"max ";}
            if(I==2){Info<<"ave ";}
            Info<<UPout[I]<<" "
                <<nutPout[I]<<" "
                <<kPout[I]<<" "
                <<yPout[I]<<" "
                <<yIBPout[I]<<" "
                <<tauWallPout[I]<<" "<<endl;
        }
    }//End of serial ru


    forAll(ibCells,I)
    {
        label cellID=ibCells[I];
        //k[cellID]=kNew[I];// k is ok
        nut[cellID]=nutNew[I];// nut is also ok
        vector Utan =(tensor::I-sqr(ibNormals[I])) & U[cellID];
        Utan /=mag(Utan)+VSMALL;
        U[cellID]=Utan*UTanNew[I];// U needs to be corrected by full scale Utau
    }

    // here tauWall is only
    wallShearStressListPtr_->set
    (
        objectID,
        tauWall
    );

    hitPointExportToMesh("yPlus",ypd,objectID);
    hitPointExportToMesh("kSample",kSample,objectID);
    hitPointExportToMesh("yIB",yIB,objectID);
    hitPointExportToMesh("ySample",ySample,objectID);
    hitPointExportToMesh("UTanNew",UTanNew,objectID);
    hitPointExportToMesh("nutNew",nutNew,objectID);
}

void Foam::immersedBoundaryFvMesh::SmagorinskyGhostCorrection
(
    const turbulenceModel& turbulence,
    const volScalarField& k,
    label objectID
)const
{
    // Uprime and Umean needs to be calculated, so as to correct U and k in LES, but k and nut only needs to modified by subgrid uTau, U needs to modified by full scale uTau, which is very essential here.

    volVectorField& U = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("U"));
    volScalarField& nut = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("nut"));
    const volScalarField& nu = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("nu"));

    const volScalarField& Gamma = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("Gamma"));


    volVectorField Umean(U);

    if(this->foundObject<volVectorField>("UMean"))
    {
        Umean = this->lookupObject<volVectorField>("UMean");
    }
    else
    {
        Info<<"Immersed boundary for LES requires UMean."<<endl;
    }

    volScalarField kFull(k);

    if(this->foundObject<volSymmTensorField>("UPrime2Mean"))
    {
        const volSymmTensorField& UPrime2Mean = this->lookupObject<volSymmTensorField>("UPrime2Mean");
        forAll(kFull,I)
        {
            scalar k2Part = std::abs(UPrime2Mean[I].xx()+UPrime2Mean[I].yy()+UPrime2Mean[I].zz());
            kFull[I]=0.5*k2Part;
        }
    }
    else
    {
        Info<<"Immersed boundary for LES requires UPrime2Mean."<<endl;
    }

    // check every cell to see if it change from dry to wet
    //- Gamma live cells (1), IB cells (0.5), ghost cells (-0.5), and dead cells(-1)
    forAll(Gamma,I)
    {
        scalar oldGamma=Gamma.oldTime()[I];
        scalar newGamma=Gamma[I];
        if(oldGamma<0 and newGamma>0) // from dead or ghost to fluid
        {
            U[I] *=0;
        }
    }

    scalar Cmu_(0.09);
    scalar kappa_(0.41);
    scalar E_(9.8);
    scalar Cs_(0.5);
    scalar Ks_ = ibProperties().lookupOrDefault<scalar>("roughnessHeight",0);
    Ks_ = objectDictList()[objectID].lookupOrDefault<scalar>("roughnessHeight",Ks_);

    scalar yPlusLam_ = 11.0;

    for (int i=0; i<10; i++)
    {
        yPlusLam_ = log(max(E_*yPlusLam_, 1))/kappa_;
    }

    const scalar Cmu25 = pow(Cmu_, 0.25);

    const pointField& ghostHitPoints = ghostHitPointsList()[objectID];
    const pointField& imagePoints = imagePointsList()[objectID];
    const labelList& ghostCells = ghostCellsList()[objectID];

    const pointField gc_centers(this->C(),ghostCells);
    vectorField ibNormals = ghostHitPoints-imagePoints;
    ibNormals=ibNormals/mag(ibNormals);

    // interpolate velocity at image points
    vectorField UIPV(imagePointsValues(U,objectID));
    vectorField UIPVn= (UIPV&ibNormals)*ibNormals;
    vectorField UIPVt=UIPV-UIPVn;
    scalarField magUIPVt = mag(UIPVt);

    // full scale k at image point
    scalarField kIPV(imagePointsValues(kFull,objectID)+SMALL);
    scalarField ksgsIPV(imagePointsValues(k,objectID)+SMALL);

    // nu at image point
    scalarField nuIPV(nu,ghostCells);

    // nut at image point
    scalarField nutIPV(imagePointsValues(nut,objectID)+SMALL);

   // image point distance
    scalarField y(mag(imagePoints-ghostHitPoints)+SMALL);
    scalarField yGhost(mag(gc_centers-ghostHitPoints)+SMALL);

    scalarField magGradUIPVt = magUIPVt/y;

    scalarField GNew(ghostCells.size(),0.0);

    scalarField kNew(k,ghostCells);

    scalarField nutNew(nut,ghostCells);
    scalarField UTanNew(magUIPVt/y*yGhost);

    vectorField tauWall(ghostCells.size(),pTraits<vector>::zero);

    Info<<"Interpolted values at image points"<<endl;
    Info<<"    "<<"U"<<tab<<"k"<<tab<<"ksgs"<<tab<<"magGradUIPVt"<<endl;
    Info<<"MIN "<<gMin(magUIPVt)<<tab<<gMin(kIPV)<<tab<<gMin(ksgsIPV)<<tab<<gMin(magGradUIPVt)<<endl;
    Info<<"AVE "<<gAverage(magUIPVt)<<tab<<tab<<gAverage(kIPV)<<tab<<gAverage(ksgsIPV)<<tab<<gAverage(magGradUIPVt)<<endl;
    Info<<"MAX "<<gMax(magUIPVt)<<tab<<tab<<gMax(kIPV)<<tab<<gMax(ksgsIPV)<<tab<<gMax(magGradUIPVt)<<endl;

    scalarField yPlusList(ghostCells.size(),0.0);
    // calculate y plus at image point using velocity
    if(Ks_>0)
    {
        // Rough Walls
        const scalar c_1 = 1/(90 - 2.25) + Ks_;
        static const scalar c_2 = 2.25/(90 - 2.25);
        static const scalar c_3 = 2.0*atan(1.0)/log(90/2.25);
        static const scalar c_4 = c_3*log(2.25);

        // If KsPlus is based on YPlus the extra term added to the law
        // of the wall will depend on yPlus
        forAll(yPlusList, facei)
        {
            const scalar Re = magUIPVt[facei]*y[facei]/nuIPV[facei];
            const scalar kappaRe = kappa_*Re;

            scalar yp = yPlusLam_;
            const scalar ryPlusLam = 1.0/yp;

            int iter = 0;
            scalar yPlusLast = 0.0;
            scalar dKsPlusdYPlus = Ks_/y[facei];

            // Additional tuning parameter - nominally = 1
            dKsPlusdYPlus *= Ks_;

            do
            {
                yPlusLast = yp;

                // The non-dimensional roughness height
                scalar KsPlus = yp*dKsPlusdYPlus;

                // The extra term in the law-of-the-wall
                scalar G = 0.0;

                scalar yPlusGPrime = 0.0;

                if (KsPlus >= 90)
                {
                    const scalar t_1 = 1 + Ks_*KsPlus;
                    G = log(t_1);
                    yPlusGPrime = Ks_*KsPlus/t_1;
                }
                else if (KsPlus > 2.25)
                {
                    const scalar t_1 = c_1*KsPlus - c_2;
                    const scalar t_2 = c_3*log(KsPlus) - c_4;
                    const scalar sint_2 = sin(t_2);
                    const scalar logt_1 = log(t_1);
                    G = logt_1*sint_2;
                    yPlusGPrime =
                        (c_1*sint_2*KsPlus/t_1) + (c_3*logt_1*cos(t_2));
                }

                scalar denom = 1.0 + log(E_*yp) - G - yPlusGPrime;
                if (mag(denom) > VSMALL)
                {
                    yp = (kappaRe + yp*(1 - yPlusGPrime))/denom;
                }
            } while
            (
                mag(ryPlusLam*(yp - yPlusLast)) > 0.0001
             && ++iter < 10
             && yp > VSMALL
            );

            yPlusList[facei] = max(0.0, yp);
        }
    }
    else
    {
        // Smooth Walls
        forAll(yPlusList, facei)
        {
            const scalar Re = magUIPVt[facei]*y[facei]/nuIPV[facei];
            const scalar kappaRe = kappa_*Re;

            scalar yp = yPlusLam_;
            const scalar ryPlusLam = 1.0/yp;

            int iter = 0;
            scalar yPlusLast = 0.0;

            do
            {
                yPlusLast = yp;
                yp = (kappaRe + yp)/(1.0 + log(E_*yp));

            } while (mag(ryPlusLam*(yp - yPlusLast)) > 0.0001 && ++iter < 10);

            yPlusList[facei] = max(0.0, yp);


        }
    }

//    // wall shear stress

    scalarField yPlusGhostList =  yPlusList*yGhost/y;

    forAll(tauWall,I)
    {
        scalar yPlusGhost=yPlusGhostList[I];
        // use full scale k to calcualte shear velocity
        scalar uTau=Cmu25*sqrt(kIPV[I]);
        uTau = magUIPVt[I]*kappa_/log(E_*yPlusGhost);

        // using instant velocity direction
        tauWall[I] = sqr(uTau)*UIPVt[I]/(magUIPVt[I] + SMALL);

        // Calculate Edash, for roughness wall function
        scalar nuLam=nuIPV[I];
        scalar KsPlus = uTau*Ks_/nuLam;
        scalar Edash=E_;
        if (KsPlus > 2.25)
        {
            scalar fnRough=1.0 + Cs_*KsPlus;
            if(KsPlus < 90.0)
            {
                fnRough = pow
                (
                  (KsPlus - 2.25)/87.75 + Cs_*KsPlus,
                  sin(0.4258*(log(KsPlus) - 0.811))
                );
            }
            Edash = E_/fnRough;
        }

        // calculate Utan
        if(yPlusGhost>yPlusLam_)
        {
            // Log-Law for tangential velocity
            UTanNew[I] = uTau*(1.0/kappa_*log(Edash*(SMALL+yPlusGhost)));

            const scalar Re = UTanNew[I]*yGhost[I]/nuIPV[I] + ROOTVSMALL;
            nutNew[I] = nuIPV[I]*(sqr(yPlusGhost)/Re - 1);
        }
        else
        {
            // Laminar sub-layer for tangential velocity: uPlus = yPlus
            UTanNew[I] = uTau*yPlusGhost;
        }
        UTanNew[I]=magUIPVt[I]*yGhost[I]/y[I];
    }


    scalarList UPout(5);
    scalarList nutPout(5);
    scalarList kPout(5);
    scalarList GPout(5);
    scalarList yPlusPout(5);
    scalarList yPlusGPout(5);
    vectorList tauWallPout(4);

    UPout[0] = gMin(UTanNew);
    nutPout[0] = gMin(nutNew);
    kPout[0] = gMin(kNew);
    GPout[0] = gMin(GNew);
    yPlusPout[0] = gMin(yPlusList);
    yPlusGPout[0] = gMin(yPlusGhostList);
    tauWallPout[0] = gMin(tauWall);

    UPout[1] = gMax(UTanNew);
    nutPout[1] = gMax(nutNew);
    kPout[1] = gMax(kNew);
    GPout[1] = gMax(GNew);
    yPlusPout[1] = gMax(yPlusList);
    yPlusGPout[1] = gMax(yPlusGhostList);
    tauWallPout[1] = gMax(tauWall);

    UPout[2] = gAverage(UTanNew);
    nutPout[2] = gAverage(nutNew);
    kPout[2] = gAverage(kNew);
    GPout[2] = gAverage(GNew);
    yPlusPout[2] = gAverage(yPlusList);
    yPlusGPout[2] = gAverage(yPlusGhostList);
    tauWallPout[2] = gAverage(tauWall);

    Info<< "UTangentialNew nutNew kNew GNew yPlus yPlusGhost tauWall" << endl;
    for (label I = 0; I < 3; I++)
    {
        if(I==0){Info<<"min ";}
        if(I==1){Info<<"max ";}
        if(I==2){Info<<"ave ";}
        Info<<UPout[I]<<" "
            <<nutPout[I]<<" "
            <<kPout[I]<<" "
            <<GPout[I]<<" "
            <<yPlusPout[I]<<" "
            <<yPlusGhostList[I]<<" "
            <<tauWallPout[I]<<" "<<endl;
    }


    // apply BC at image point to ghost cells
    // make sure grad Utan between this two point equals to tauWall/(nut+nu)

    vector Uw0=pTraits<vector>::zero;

    volVectorField& U_desired = const_cast<volVectorField&>
        (this->lookupObject<volVectorField>("U_desired"));

    U_desired==U;

    const labelList& idc=ibDeadCellsList()[objectID];

    forAll(idc,I)
    {
        U_desired[idc[I]] *=0;
    }
    // set value to ghost cells
    forAll(ghostCells,I)
    {
        label cellID=ghostCells[I];
        nut[cellID]=nutNew[I];


        const vector& Uin=UIPVn[I];


        U_desired[cellID] = (Uw0-Uin)+Uw0; // add normal part

    }

    // here tauWall is only
    wallShearStressListPtr_->set
    (
        objectID,
        tauWall
    );

    hitPointExportToMesh("kNew",kNew,objectID);
    hitPointExportToMesh("nutNew",nutNew,objectID);
    hitPointExportToMesh("yPlus",yPlusList,objectID);

}


void immersedBoundaryFvMesh::makeFunctionObjectList( const label& objectID) const
{
    dictionary FOdict(objectDictList()[objectID]);
    functionObjectList FOList(time(),FOdict);
}
// ************************************************************************* //

