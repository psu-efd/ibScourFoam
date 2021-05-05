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
Foam::scalar Foam::wallFuncModules::uPlus
(
    const scalar& yPlus,
    const scalar& ksPlus
)const
{
    if(name()=="two-layer")
    {
        return two_layer(yPlus,max(1.0,ksPlus));
    }
    else if(name()=="Spalding")
    {
        return inverse_spalding(yPlus,max(1.0,ksPlus));
    }
    else if(name()=="vanDriest")
    {
        return vanDriest(yPlus,max(1.0,ksPlus));
    }
    else
    {
        return two_layer(yPlus,ksPlus);
    }
}

Foam::scalar Foam::wallFuncModules::dupdyp
(
    const scalar& yPlus,
    const scalar& ksPlus
)const
{
    if(name()=="two-layer")
    {
        return two_layer_dudy(max(1.0,yPlus),max(1.0,ksPlus));
    }
    else if(name()=="Spalding")
    {
        return 1.0/(SMALL+spalding_dydu(inverse_spalding(max(1.0,yPlus),ksPlus),max(1.0,ksPlus)));
    }
    else if(name()=="vanDriest")
    {
        return vanDriest_dudy(max(1.0,yPlus),max(1.0,ksPlus));
    }
    else
    {
        return two_layer_dudy(max(1.0,yPlus),ksPlus);
    }
}

Foam::scalar Foam::wallFuncModules::Edash
(
    const scalar& ksPlus
)const
{
    if(ksPlus<2.25)
    {
        return E_;
    }

    else
    {
        if(debug>2)
        {
            if(ksPlus>90.0)
            {
                WarningInFunction
                    << "ksPlus is "<<ksPlus<<" larger than "
                    << 90.0<< endl;
            }
        }
        return E_/pow((ksPlus - 2.25)/87.75 + Cs_*ksPlus,sin(0.4258*(log(ksPlus) - 0.811)));
    }

};

Foam::scalar Foam::wallFuncModules::yPlusLam
(
    const scalar kappa,
    const scalar E
)const
{
     scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(max(E*ypl, 1))/kappa;
    }

    return ypl;

}


Foam::scalar Foam::wallFuncModules::yPlusLam
(
    const scalar& ksPlus
)const
{
    scalar yPlusLam_=11.1;
    scalar error=1.0;
    scalar counter=0;
    scalar oldYPlusLam_=yPlusLam_;
    scalar Edash=this->Edash(ksPlus);
    while (error>1.0e-3 and counter < 10)
    {
        oldYPlusLam_=yPlusLam_;
        yPlusLam_=1.0/kappa_*log(max(Edash*yPlusLam_,1.0));
        error=mag((oldYPlusLam_-yPlusLam_)/(oldYPlusLam_+SMALL));
        counter++;
    }
    return yPlusLam_;
}

Foam::scalar Foam::wallFuncModules::two_layer
(
    const scalar& yPlus,
    const scalar& ksPlus
)const
{

    if(yPlus<ksPlus)
    {
        return 0.0;
    }
    else if(yPlus<yPlusLam(ksPlus))
    {
        return yPlus;
    }
    else
    {
        return 1.0/kappa_*log(Edash(ksPlus)*yPlus);
    }
}



Foam::scalar Foam::wallFuncModules::two_layer_dudy
(
    const scalar& yPlus,
    const scalar& ksPlus
)const
{
    scalar B = 5.2;
    scalar DeltaB = 0;

    if(ksPlus <= 2.5)
    {    
       DeltaB = 0;
    }
    else if(ksPlus>2.5 and ksPlus<90)
    {
       DeltaB = (B-8.5+1/0.41*log(ksPlus))*sin(0.4258*(log(ksPlus)-0.811));
    }
    else
    {
       DeltaB = B-8.5+1/0.41*log(ksPlus);
    }


    return 1/0.41*log(yPlus)+B-DeltaB; 
}








Foam::scalar Foam::wallFuncModules::spalding  // return yPlus
(
    const scalar& uPlus,
    const scalar& ksPlus
)const
{
    scalar D=log(Edash(ksPlus))/kappa_;
    scalar kappaUPlus=kappa_*uPlus;
    return uPlus+exp(-kappa_*D)*
        (
            exp(kappaUPlus)-1-kappaUPlus-sqr(kappaUPlus)/(1.0*2.0)-pow3(kappaUPlus)/(1.0*2.0*3.0)
        );
}

Foam::scalar Foam::wallFuncModules::spalding_dydu
(
    const scalar& uPlus,
    const scalar& ksPlus
)const
{
    scalar D=log(Edash(ksPlus))/kappa_;
    scalar partA=exp(-kappa_*D)*
        (
            kappa_*exp(kappa_*uPlus)-kappa_-sqr(kappa_)*uPlus-pow3(kappa_)*sqr(uPlus)/(1.0*2.0)
        );
    scalar partB=(spalding(uPlus,ksPlus)-uPlus)*kappa_*2.5/uPlus;
    return 1.0+partA+partB;
}

Foam::scalar Foam::wallFuncModules::inverse_spalding
(
    const scalar& yPlus,
    const scalar& ksPlus
)const
{
    // estimate uplus
    scalar up = two_layer(yPlus,ksPlus);

    scalar error=1.0;
    scalar counter=0;
    scalar oldUp=up;
    while(error>1e-4 and counter<100)
    {
        oldUp=up;
        scalar newF=yPlus-spalding(up,ksPlus);
        up+=SMALL;
        scalar dNewF=-yPlus/up-spalding_dydu(up,ksPlus);
        dNewF+=SMALL;
        up=up-newF/dNewF*0.75;
        oldUp=oldUp+SMALL;
        error=mag((oldUp-up)/oldUp);
        counter++;
    }

    return up;
}

Foam::scalar Foam::wallFuncModules::vanDriest_dudy
(
    const scalar& yPlus,
    const scalar& ksPlus
)const
{
    scalar yPlusR=yPlus;
    yPlusR=max(0.0,yPlusR);
    if(ksPlus>5.0)
    {
        yPlusR=yPlusR+0.9*(sqrt(ksPlus)-ksPlus*exp(-ksPlus/6.0));
    }
    scalar partA=1.0-exp(-yPlusR/25.0);
    scalar partB=sqr(2*kappa_*yPlusR);
    return 2.0*1.0/(1.0+sqrt(1.0+sqr(partA)*partB));
}

Foam::scalar Foam::wallFuncModules::vanDriest
(
    const scalar& yPlus,
    const scalar& ksPlus
)const
{


    if(ksPlus<SMALL)
    {
        if(yPlus<11.0)
        {
             return yPlus;
        }
        else
        {
            return 1.0/kappa_*log(9.0*yPlus);
        }
   }
   else
   {

        scalar NN=max(yPlus*0.5,30);
        scalar dyPlus=yPlus/(NN+1);
       scalar i=0;
        scalar SUM=0.0;
        while(i<NN)
        {
            SUM +=
                dyPlus*vanDriest_dudy((i+(i+1))*0.5*dyPlus,ksPlus);
            i++;
        }
        return SUM;
        
        if(yPlus<ksPlus)
       {
           return 0.0;
       }
   }
}


// ************************************************************************* //
