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

#include "wallFuncModules.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "vectorTools.H"
#include "simpleObjectRegistry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
 defineTypeNameAndDebug(wallFuncModules, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallFuncModules::wallFuncModules
(
    const word& name,
    const word& turName
)
:
name_(name),
turbName_(turName),
Cs_(0.5),
Cmu_(0.09),
Cmu25(pow(Cmu_, 0.25)),
Cmu50(pow(Cmu_, 0.5)),
Cmu75(pow(Cmu_, 0.75)),
E_(9.8),
kappa_(0.41),
nu_(1.0e-6),
beta1_(0.075),
yPlusLam_(yPlusLam(kappa_, E_))

{
    Info<<"Use "<<name<<" wall function"<<endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //



Foam::wallFuncModules::~wallFuncModules()
{
    clearOut();
}


void Foam::wallFuncModules::clearOut()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::wallFuncModules::calcUtau
(
     const scalarField& USampleMag, // tangential velocity at sample point
     const scalarField& kSample, // tke at sample point
     const scalarField& ySample, // sample point wall distance
     const scalar& Ks_, // roughness height
     scalarField& uTau
)const
{
    forAll(uTau,I)
    {
        uTau[I]=this->uTau
            (
                USampleMag[I],
                kSample[I],
                ySample[I],
                Ks_
            );
    }
}




Foam::scalar Foam::wallFuncModules::yPlus
(
     const scalar& USampleMag, // tangential velocity at sample point
     const scalar& ySample, // sample point wall distance
     const scalar& Ks_, // roughness height
     const scalar& roughnessFactor_,
     const scalar& roughnessConstant_

 )const
{
    if(Ks_>0)
    {
   
	    const scalar c_1 = 1/(90 - 2.25) + roughnessConstant_;
        static const scalar c_2 = 2.25/(90 - 2.25);
        static const scalar c_3 = 2.0*atan(1.0)/log(90/2.25);
        static const scalar c_4 = c_3*log(2.25);

	    scalar kappaRe = kappa_*USampleMag*ySample/nu_;
	    scalar yp = yPlusLam_;
	    scalar ryPlusLam = 1.0/yp;
	
	    int iter = 0;
	    scalar yPlusLast = 0.0;
	    scalar dKsPlusdYPlus = Ks_/ySample;

	    dKsPlusdYPlus *= roughnessFactor_;
	
	    do
	    {
		    yPlusLast = yp;
		    scalar KsPlus = yp*dKsPlusdYPlus;
		    scalar G = 0.0;

                    scalar yPlusGPrime = 0.0;
      		if (KsPlus >= 90)
                    {
                            const scalar t_1 = 1 + roughnessConstant_*KsPlus;
                            G = log(t_1);
                            yPlusGPrime = roughnessConstant_*KsPlus/t_1;
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


	    }while	
	    (
                mag(ryPlusLam*(yp - yPlusLast)) > 0.0001
             && ++iter < 10
             && yp > VSMALL
	    );
            
            scalar yPlus = max (0.0, yp);
            return yPlus;
    }
    else
    {
        scalar Re = USampleMag*ySample/nu_;
        scalar kappaRe = kappa_*Re;
        scalar yp = yPlusLam_;
	    scalar ryPlusLam = 1.0/yp;
	    int iter = 0;
	    scalar yPlusLast = 0.0;
        do
	    {
		    yPlusLast = yp;
		    
            
            yp = (kappaRe + yp)/(1 + log(E_*yp));
                   

	    }while	
	    (
                mag(ryPlusLam*(yp - yPlusLast)) > 0.0001
             && ++iter < 10
             && yp > VSMALL
	    );
            

            scalar yPlus = max (0.0, yp);
            return yPlus;
    }
}    







Foam::scalar Foam::wallFuncModules::uTau
(
     const scalar& USampleMag, // tangential velocity at sample point
     const scalar& kSample, // tke at sample point
     const scalar& ySample, // sample point wall distance
     const scalar& Ks_ // roughness height
)const
{
    scalar error=1.0;
    scalar counter=0;

    if(USampleMag<1.0e-7)
    {
        return 0.0;
    }

    //estimate initial utau
    scalar utau = Cmu25*sqrt(kSample);
    scalar oldUtau=utau;
    scalar ksPlus=Ks_*utau/nu_;

    scalar yPlus=ySample*utau/nu_;
    if(name()=="two-layer" or name()=="vanDriest")
    {
        while(error>0.5e-2 and counter <40)
        {
            oldUtau=utau;
            ksPlus=Ks_*utau/nu_;
            yPlus = max(1.0,ySample*utau/nu_);
            scalar tmp_up=this->uPlus(yPlus,ksPlus);
    

            if(tmp_up==0)
            {
                utau=0.0;
                error=0.0;
            }
            else
            {
                utau=(USampleMag/tmp_up);
                error=mag((oldUtau-utau)/oldUtau);
            }
            counter++;
        }

    }
    else if(name()=="Spalding")
    {
        scalar uPlus=USampleMag/utau;
        while(error>0.5e-2 and counter <20)
        {
            oldUtau=utau;
            ksPlus=max(1.0,Ks_*utau/nu_);
            uPlus=USampleMag/utau;
            uPlus=min(uPlus,50.0/kappa_);
            yPlus = max(1.0,ySample*utau/nu_);
            scalar newF=yPlus-spalding(uPlus,ksPlus);
            scalar dNewF=ySample/nu_+spalding_dydu(uPlus,ksPlus)*USampleMag/sqr(utau);
            utau=utau-newF/dNewF;
            oldUtau+=SMALL;
            error=mag((oldUtau-utau)/oldUtau);
            counter++;
        }
    }
    else
    {
        this->name_="two-layer";
        utau=this->uTau
            (
                USampleMag,
                kSample,
                ySample,
                Ks_
      
            );
    }

    return utau;
}

  Foam::scalar Foam::wallFuncModules::estimateUTau
 (
     scalar& U_edge, // meanVelocity for incoming flow / edge velocity
     const scalar& waterDepth, // water depth
     const scalar& boundaryLayerThickness, // boundaryLayerThickness
     const scalar& nuLam, // nuLam
     const scalar& Ks_ // roughness height

 )const
 {
    scalar uTau=0;
    scalar counter=0;
    scalar Umean_esti=0;
    const scalar meanVelocity=U_edge;

    while (mag(Umean_esti-meanVelocity)>0.0001 and counter<100)
    {
        counter++;
        uTau=this->uTau(U_edge,0.001,boundaryLayerThickness,Ks_);

        // calculate mean velocity based on this estimated uTau
        scalarField U(100,0);
        for(label I = 1; I < 100; I++)
        {
            scalar H=I/100.0*waterDepth;
            H=min(H,boundaryLayerThickness);
            scalar yPlus = max(1.0,H*uTau/nuLam);
            scalar ksPlus = max(1.0,Ks_*uTau/nuLam);
            scalar uPlus = this->uPlus(yPlus,ksPlus);
            U[I]=uPlus*uTau;
        }
        Umean_esti=average(U);
        U_edge=meanVelocity/Umean_esti*U_edge;
    }
    return uTau;
 }
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "velocityProfiles.C"
#include "turbModels.C"
// ************************************************************************* //
