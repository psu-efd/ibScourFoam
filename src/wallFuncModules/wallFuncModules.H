 /*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
     \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
      \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
     Foam::wallFuncModules

 Description
     Provide wall functions options for immersed boundary method
     Different velocity profiles to calc utau
         two-layer wall model
         van Driest
         Spalding
     k-omega estimate omega_w
         Roulund2005
         Baykal2015

 SourceFiles
     wallFuncModules.C

 \*---------------------------------------------------------------------------*/

#ifndef wallFuncModules_H
#define wallFuncModules_H


#include "DynamicList.H"
#include "meshSearch.H"
#include "SubField.H"
#include "volFields.H"
#include "fvCFD.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class wallFuncModules Declaration
\*---------------------------------------------------------------------------*/

class wallFuncModules
{
public:

private:
 // Private member data

     //- velocity profile name
     mutable word name_;

     //- turbulenceModel name
     mutable word turbName_;

         // constants
         mutable scalar Cs_;
         mutable scalar Cmu_;
         mutable scalar Cmu25;
         mutable scalar Cmu50;
         mutable scalar Cmu75;
         mutable scalar E_;
         mutable scalar kappa_;
         mutable scalar nu_;
         mutable scalar beta1_;
	 mutable scalar yPlusLam_;


     //- Disallow default bitwise copy construct
     wallFuncModules(const wallFuncModules&);

     //- Disallow default bitwise assignment
     void operator=(const wallFuncModules&);

    void makeBoundaryFacesAndEdges
    (
    ) const;

public:

 //- Runtime type information
 TypeName("wallFuncModules");

 //- Return the mesh sub-directory name (usually "triSurface")
 static word meshSubDir;


 // Constructors

     //- Construct from velocity profile name  1-two-layer 2-vanDriest 3-Spalding
     wallFuncModules
     (
         const word& name,
         const word& turbName
     );

 //- Destructor
 virtual ~wallFuncModules();

     //- Clear storage
     void clearOut();


 // Member Functions
     //- name
     const word& name() const
     {
         return name_;
     }

     word& turbName() const
     {
         return turbName_;
     }

     //- calculate nut
     void calcUtau
     (
         const scalarField& USampleMag, // tangential velocity at sample point
         const scalarField& kSample, // tke at sample point
         const scalarField& ySample, // sample point wall distance
         const scalar& Ks_, // roughness height
         scalarField& uTau
     )const;

     scalar uTau
     (
         const scalar& USampleMag, // tangential velocity at sample point
         const scalar& kSample, // tke at sample point
         const scalar& ySample, // sample point wall distance
         const scalar& Ks_ // roughness height
     )const;
    scalar yPlus
    (
         const scalar& USampleMag, // tangential velocity at sample point
         const scalar& ySample, // sample point wall distance
         const scalar& Ks_, // roughness height
         const scalar& roughnessFactor_,
         const scalar& roughnessConstant

    )const;

    scalar yPlusLam
    (
        const scalar kappa,
        const scalar E
    )const;
    

     scalar uPlus
     (
         const scalar& yPlus,
         const scalar& ksPlus
     )const;

     scalar Edash
     (
         const scalar& ksPlus
     )const;

     scalar yPlusLam
     (
         const scalar& ksPlus
     )const;

     scalar two_layer
     (
         const scalar& yPlus,
         const scalar& ksPlus
     )const;

         scalar two_layer_dudy
         (
             const scalar& yPlus,
             const scalar& ksPlus
         )const;

     scalar spalding // uplus
     (
         const scalar& uPlus,
         const scalar& ksPlus
     )const;

         scalar spalding_dydu //derivative
         (
             const scalar& uPlus,
             const scalar& ksPlus
         )const;

         scalar inverse_spalding
         (
             const scalar& yPlus,
             const scalar& ksPlus
         )const;

     scalar vanDriest
     (
         const scalar& yPlus,
         const scalar& ksPlus
     )const;

         scalar vanDriest_dudy
         (
             const scalar& yPlus,
             const scalar& ksPlus
         )const;

    scalar dupdyp
    (
         const scalar& yPlus,
         const scalar& ksPlus
    )const;

     //- estimate utau based on incoming mean velocity and water depth
      scalar estimateUTau
     (
         scalar& U_edge, // meanVelocity for incoming flow / edge velocity
         const scalar& waterDepth, // water depth
         const scalar& boundaryLayerThickness, // boundaryLayerThickness
         const scalar& nuLam, // nuLam
         const scalar& Ks_ // roughness height
     )const;


    // turbulence models
    void kOmegaUpdate
    (
        const scalar& y,
        const scalar& uTau,
        const scalar& ksPlus,
        scalar& k,
        scalar& omega,
        scalar& nut
    )const;

        void Omega_Wilcox2008
        (
            const scalar& y,
            const scalar& uTau,
            const scalar& ksPlus,
            scalar& k,
            scalar& omega,
            scalar& nut
        )const;

        void Omega_Menter2001
        (
            const scalar& y,
            const scalar& uTau,
            const scalar& ksPlus,
            scalar& k,
            scalar& omega,
            scalar& nut
        )const;

        void Omega_Baykal2015
        (
            const scalar& y,
            const scalar& uTau,
            const scalar& ksPlus,
            scalar& k,
            scalar& omega,
            scalar& nut
        )const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
