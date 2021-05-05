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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sedModules::roulund2005
(
    const scalar& D50, // grain size
    const scalar& ssG, // specific Gravity
    const scalar& mus, // static friction coefficient
    const scalar& mud, // dynamic friction coefficient
    const vectorField& triNormals, // triangle normals
    const vectorField& wss, // wall shear stress
    const scalar& sNC0, // critical shields number for flat bed
    const vector& gravity, // gravity
    scalarField& beta, // bed slope (to be updated)
    scalarField& phi, // angle between utau and steepest bed slope (to be updated)
    scalarField& sN, // magnitude of shields number
    scalarField& sNC, // critical shields number (to be updated)
    vectorField& q0 // bed load flux (to be updated)

)const
{
    vector gravityNormal=gravity/mag(gravity);
    sN=mag(wss)/((ssG-1)*mag(gravity)*D50);
    forAll(beta,I)
    {
        // compute beta
        vector fN=triNormals[I];

         roulund2005
         (
             D50, // grain size
             ssG, // specific Gravity
             mus, // static friction coefficient
             mud, // dynamic friction coefficient
             fN, // bed normal
             wss[I], // wall shear stress
             sNC0, // critical shields number for flat bed
             gravity, // gravity
             gravityNormal, // gravityNormal
             beta[I], // bed slope (to be updated)
             phi[I], // direction between sN_dir and steepest bed slope (to be updated)
             sN[I], // magnitude of shields number
             sNC[I], // critical shields number (to be updated)
             q0[I] // bed load flux (to be updated)

         );
    }

}




void Foam::sedModules::roulund2005
(
     const scalar& D50, // grain size
     const scalar& ssG, // specific Gravity
     const scalar& mus, // static friction coefficient
     const scalar& mud, // dynamic friction coefficient
     const vector& fN, // bed normal
     const vector& wss, // wall shear stress
     const scalar& sNC0, // critical shields number for flat bed
     const vector& gravity, // gravity
     const vector& gravityNormal, // gravityNormal
     scalar& beta, // bed slope (to be updated)
     scalar& phi, // direction between sN_dir and steepest bed slope (to be updated)
     scalar& sN, // magnitude of shields number
     scalar& sNC, // critical shields number (to be updated)
     vector& q0 // bed load flux (to be updated)
)const
{
    scalar pi_value=3.141592653;

    vector xyPlaneNormal(0,0,1);
    beta = Foam::acos(fN & xyPlaneNormal);

    if(beta<0.0)  beta += pi_value;
    if(beta>pi_value/2.0)  beta -=pi_value;

    if(beta < 0) beta = 0;

    if(beta > mus) beta = mus;


    //calculate maximum slope       s=-(g*n)/(n*n)*n+g
    vector slope = -(gravity&fN)*fN+gravity;


    if((slope&gravity)<0) slope*=-1.0;      

    if(mag(slope)<1e-6)
    { 
        slope = vector(1,0,0);
    }
    slope = slope/(mag(slope)+SMALL);
    

    if(mag(wss)<1e-6)
    {
	phi = 0;
    } 
    else
    {
        phi = Foam::acos((wss & slope)/(mag(wss)+SMALL)/(mag(slope)+SMALL));
    }

    if(mag(wss)<SMALL) phi=0;

    // Roulund et al 2005, Eq. 25
    scalar sqrt_base=1-pow(Foam::sin(phi)*Foam::tan(beta)/mus,2);

    sqrt_base=max(sqrt_base,0);

    sNC = sNC0*
        (
            Foam::cos(beta)
           *Foam::sqrt(sqrt_base)
           -Foam::cos(phi)*Foam::sin(beta)/mus
        );

    sNC=max(sNC,0);
    scalar qstar=0;
    if(sN>sNC)
    {
        qstar=18.74*(sN-sNC)*(sqrt(sN)-0.7*sqrt(sNC))*sqrt((ssG-1)*mag(gravity)*D50)*D50;      
    }
    else
    {
        qstar=0;
    }
    q0 = qstar*wss/(mag(wss)+SMALL);
   

}


void Foam::sedModules::roulund2005
(
     const scalar& D50, // grain size
     const scalar& ssG, // specific Gravity
     const scalar& mus, // static friction coefficient
     const scalar& mud, // dynamic friction coefficient
     const vector& fN, // bed normal
     const vector& wss, // wall shear stress
     const vector& Vs, // settling velocity
     const scalar& Cb, // near-bed concentration
     const scalar& sNC0, // critical shields number for flat bed
     const vector& gravity, // gravity
     const vector& gravityNormal, // gravityNormal
     scalar& Zb, // referenceHeight (to be updated)
     scalar& beta, // bed slope (to be updated)
     scalar& phi, // direction between sN_dir and steepest bed slope (to be updated)
     scalar& sN, // magnitude of shields number
     scalar& sNC, // critical shields number (to be updated)
     vector& q0, // bed load flux (to be updated)
     scalar& E, // entrainment rate (to be updated)
     scalar& D // deposition rate (to be updated)
)const
{
    scalar pi_value=3.141592653;


    vector slope = -(gravity&fN)/(fN&fN)*fN+gravity;
    if(mag(slope)<1e-6)
    { 
        slope = wss - (gravityNormal&wss)*gravityNormal;
    }
    slope = slope/(mag(slope)+SMALL);
    if((slope&gravity)<0) slope*=-1.0;
    beta = Foam::vectorTools::radAngleBetween(slope,gravity);
    beta = degToRad(90-radToDeg(beta));
    if(beta> 0.5*pi_value) beta=0.5*pi_value-beta;

    phi = Foam::vectorTools::radAngleBetween(slope,wss);

    if(mag(wss)<SMALL) phi=0;

    // Roulund et al 2005, Eq. 25
    scalar sqrt_base=1-pow(Foam::sin(phi)*Foam::tan(beta)/mus,2);
    sqrt_base=max(sqrt_base,0);
    sNC = sNC0*
        (
            Foam::cos(beta)
           *Foam::sqrt(sqrt_base)
           -Foam::cos(phi)*Foam::sin(beta)/mus
        );

    // only if phi>repose angle
    sNC=max(sNC,0);
    scalar pEF = 0.0;
    vector utau=pTraits<vector>::zero;
    if(mag(wss)>0)
    {
        utau=Foam::sqrt(mag(wss))*wss/mag(wss);
    }
    scalar Cbstar=0;
    if(sN>sNC)
    {
        // Roulund et al 2005, Eq. 24
        pEF=1.0/6.0*pi_value*mud/(SMALL+sN-sNC);
        pEF=pow(1+pow(pEF,4),-0.25);

        // Roulund et al 2005, Eq. 23
        vector ub=10.0*(1.0-0.7*Foam::sqrt(sNC/sN))*utau;
        q0=1.0/6.0*pi_value*D50*pEF*ub;
        Cbstar=1.0/6.0*pi_value*D50*pEF;
    }
    else
    {
        q0*=0;
    }
    vector q0Tan = q0-(q0&gravityNormal)*gravityNormal;
    q0=q0Tan/(SMALL+mag(q0Tan))*mag(q0);
    phi=radToDeg( phi);


    // deposition rate
    D=Cb*(Vs&gravityNormal);

    // entrainment rate
    E=Cbstar*(Vs&gravityNormal);
}












void Foam::sedModules::EDUpdateValues
(
    const scalar& D50, // grain size
     const scalar& ssG, // specific Gravity
     const scalar& mus, // static friction coefficient
     const scalar& mud, // dynamic friction coefficient
     const vectorField& fN, // bed normal
     const vectorField& wss, // wall shear stress
     const vectorField& Vs, // settling velocity
     const scalarField& Cb, // near-bed concentration
     const scalar& sNC0, // critical shields number for flat bed
     const vector& gravity, // gravity
     scalarField& Zb, // referenceHeight (to be updated)
    scalarField& E, // entrainment rate (to be updated)
    scalarField& D // deposition rate (to be updated)
)const
{
    vector gravityNormal=gravity/mag(gravity);
    scalarField sN=mag(wss)/((ssG-1)*mag(gravity)*D50);
    forAll(E,I)
    {
        scalar beta=0;
        scalar phi=0;
        scalar sNC=0;
        vector q0=pTraits<vector>::zero;
        roulund2005
        (
            D50, // grain size
            ssG, // specific Gravity
            mus, // static friction coefficient
            mud, // dynamic friction coefficient
            fN[I], // bed normal
            wss[I], // wall shear stress
            Vs[I], // settling velocity
            Cb[I], // near-bed concentration
            sNC0, // critical shields number for flat bed
            gravity, // gravity
            gravityNormal, // gravityNormal
            Zb[I], // referenceHeight (to be updated)
            beta, // bed slope (to be updated)
            phi, // direction between sN_dir and steepest bed slope (to be updated)
            sN[I], // magnitude of shields number
            sNC, // critical shields number (to be updated)
            q0, // bed load flux (to be updated)
            E[I], // entrainment rate (to be updated)
            D[I] // deposition rate (to be updated)
        );
    }
}
// ************************************************************************* //
