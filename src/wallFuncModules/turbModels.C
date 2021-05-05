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
void Foam::wallFuncModules::kOmegaUpdate
(
    const scalar& y,
    const scalar& uTau,
    const scalar& ksPlus,
    scalar& k,
    scalar& omega,
    scalar& nut
)const
{
//    word omega_name="Menter2001";
    word omega_name=turbName();
    // Wilcox2008
    if(omega_name=="Wilcox2008")
    {
        Omega_Wilcox2008(y,uTau,ksPlus,k,omega,nut);
    }
    else if (omega_name=="Menter2001")
    {
        Omega_Menter2001(y,uTau,ksPlus,k,omega,nut);
    }
    else if (omega_name=="Baykal2015")
    {
        Omega_Baykal2015(y,uTau,ksPlus,k,omega,nut);
    }
    else
    {
        turbName()="Baykal2015";
        kOmegaUpdate(y,uTau,ksPlus,k,omega,nut);
    }
}


void Foam::wallFuncModules::Omega_Wilcox2008
(
    const scalar& y,
    const scalar& uTau,
    const scalar& ksPlus,
    scalar& k,
    scalar& omega,
    scalar& nut
)const
{
    // update omega
    scalar ksPlus_tmp=100.0/ksPlus;
    omega = sqr(ksPlus_tmp*2);
    if(ksPlus>5)
    {
        omega=ksPlus_tmp+(sqr(2*ksPlus_tmp)-ksPlus_tmp)*exp(5.0-ksPlus);
    }

    // update k
    k=omega*nut;
}


void Foam::wallFuncModules::Omega_Menter2001
(
    const scalar& y,
    const scalar& uTau,
    const scalar& ksPlus,
    scalar& k,
    scalar& omega,
    scalar& nut
)const
{
    if(y<SMALL)
    {
        Omega_Baykal2015(y,uTau,ksPlus,k,omega,nut);
        return;
    }
    if(nut==0)
    {
        k=1e-10;
        omega=6.0*nu_/(beta1_*sqr(y));
        return;
    }
    // update k
    k=Cmu25*kappa_*y/nut;

    scalar yPlus=y*uTau/nu_;

    const scalar omegaVis = 6.0*nu_/(beta1_*sqr(y));
    const scalar omegaLog = sqrt(k)/(Cmu25*kappa_*y);

    if(yPlus>yPlusLam(ksPlus))
    {
        omega=omegaLog;
    }
    else
    {
        omega=omegaVis;
    }

    omega=sqrt(sqr(omegaVis) + sqr(omegaLog));
    // update k
    k=omega*nut;
}

void Foam::wallFuncModules::Omega_Baykal2015
(
    const scalar& y,
    const scalar& uTau,
    const scalar& ksPlus,
    scalar& k,
    scalar& omega,
    scalar& nut
)const
{
    scalar yPlus=y*uTau/nu_;

    yPlus=max(1.0,yPlus);

    omega=sqr(uTau)/nu_*max(96.885/sqr(yPlus),1.0/sqrt(0.09)/kappa_/yPlus);
    // update k
    k=omega*nut;
}
// ************************************************************************* //
