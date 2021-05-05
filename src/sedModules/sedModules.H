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
     Foam::sedModules

 Description
     Provide sediment transport modules, processing functions

 SourceFiles
     sedModules.C

 \*---------------------------------------------------------------------------*/

#ifndef sedModules_H
#define sedModules_H


#include "DynamicList.H"
#include "meshSearch.H"
#include "SubField.H"
#include "volFields.H"
#include "fvCFD.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sedModules Declaration
\*---------------------------------------------------------------------------*/

class sedModules
{
public:

private:
 // Private member data

     //- include faces
     mutable word* namePtr_;

     //- Disallow default bitwise copy construct
     sedModules(const sedModules&);

     //- Disallow default bitwise assignment
     void operator=(const sedModules&);

     //- update values
     void roulund2005
     (
         const scalar& D50, // grain size
         const scalar& ssG, // specific Gravity
         const scalar& mus, // static friction coefficient
         const scalar& mud, // dynamic friction coefficient
         const vectorField& gradEta, // bed elevation gradient
         const vectorField& wss, // wall shear stress
         const scalar& sNc0, // critical shields number for flat bed
         const vector& gravity, // gravity
         scalarField& beta, // bed slope (to be updated)
         scalarField& phi, // direction between sN_dir and steepest bed slope (to be updated)
         scalarField& sN, // magnitude of shields number
         scalarField& sNc, // critical shields number (to be updated)
         vectorField& q0 // bed load flux (to be updated)
     )const;

     void roulund2005
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
     )const;
     
     void roulund2005
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
         scalar& Zb, //referenceHeight (to be updated)
         scalar& beta, // bed slope (to be updated)
         scalar& phi, // direction between sN_dir and steepest bed slope (to be updated)
         scalar& sN, // magnitude of shields number
         scalar& sNC, // critical shields number (to be updated)
         vector& q0, // bed load flux (to be updated)
         scalar& E, // entrainment rate (to be updated)
         scalar& D // deposition rate (to be updated)
     )const;


     //- update values
     void EDUpdateValues
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
         scalarField& Zb, //referenceHeight (to be updated)
        scalarField& E, // entrainment rate (to be updated)
        scalarField& D // deposition rate (to be updated)
     )const;

    void makeBoundaryFacesAndEdges
    (
    ) const;

public:

 //- Runtime type information
 TypeName("sedModules");

 //- Return the mesh sub-directory name (usually "triSurface")
 static word meshSubDir;


 // Constructors

     //- Construct from triSurface
     sedModules(const word& name);

 //- Destructor
 virtual ~sedModules();

     //- Clear storage
     void clearOut();


 // Member Functions
     //- name
     const word& name() const
     {
         return *namePtr_;
     }

     //- update values
     void bedLoadUpdateValues
     (
         const scalar& D50, // grain size
         const scalar& ssG, // specific Gravity
         const scalar& mus, // static friction coefficient
         const scalar& mud, // dynamic friction coefficient
         const vectorField& gradEta, // bed elevation gradient
         const vectorField& wss, // wall shear stress
         const scalar& sNC0, // critical shields number for flat bed
         const vector& gravity, // gravity
         scalarField& beta, // bed slope (to be updated)
         scalarField& phi, // direction between sN_dir and steepest bed slope (to be updated)
         scalarField& sN, // magnitude of shields number
         scalarField& sNC, // critical shields number (to be updated)
         vectorField& q0 // bed load flux (to be updated)
     )const;
     
    void bedLoadUpdateValues
    (
         const scalar& D50, // grain size
         const scalar& ssG, // specific Gravity
         const scalar& mus, // static friction coefficient
         const scalar& mud, // dynamic friction coefficient
         const surfaceVectorField& triNormals, // triangle normals
         const surfaceVectorField& wss, // wall shear stress
         const scalar& sNC0, // critical shields number for flat bed
         const vector& gravity, // gravity
         surfaceScalarField& beta, // bed slope (to be updated)
         surfaceScalarField& phi, // direction between sN_dir and steepest bed slope (to be updated)
         surfaceScalarField& sN, // magnitude of shields number
         surfaceScalarField& sNC, // critical shields number (to be updated)
         surfaceVectorField& q0 // bed load flux (to be updated)
    )const;
    
    void bedLoadUpdateValues
    (
         const scalar& D50, // grain size
         const scalar& ssG, // specific Gravity
         const scalar& mus, // static friction coefficient
         const scalar& mud, // dynamic friction coefficient
         const volVectorField& triNormals, // triangle normals
         const volVectorField& wss, // wall shear stress
         const scalar& sNC0, // critical shields number for flat bed
         const vector& gravity, // gravity
         volScalarField& beta, // bed slope (to be updated)
         volScalarField& phi, // direction between sN_dir and steepest bed slope (to be updated)
         volScalarField& sN, // magnitude of shields number
         volScalarField& sNC, // critical shields number (to be updated)
         volVectorField& q0 // bed load flux (to be updated)
  
    )const;    

     //- update values entrainment and deposition
     void ibEDUpdateValues
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
     )const;



     // sand-slide algorithm
     // diffusion-based
     bool sandSlide
     (
           volScalarField& eta, //
        const scalar& mus, // static friction coefficient
        const scalar& K0, // diffusivity coefficient
        const scalar& error1,
        const scalar& error2,
        const scalar& timeStep,
        const scalar& maxIter
     ) const;

     

     // sand-slide algorithm
     // geometric-based
     bool sandSlide
     (
           volScalarField& eta, //
        const scalar& mus, // static friction coefficient
        const scalar& alpha, // acceleration on mus
        const scalar& beta, // relaxation factor
        const scalar& maxIter
     ) const;

     //- fix bed load flux accross wall  sum==0;
     void fixWallFlux
     (
        surfaceScalarField& qFlux // bed load flux
     )const;

     //- change inlet flux to ensure inlet eta does not change
     void feedInlet
     (
        surfaceScalarField& qFlux // bed load flux
     )const;

     //- check mass balance of bed changes = bedLoadflux + deposition - entrainment
     void checkMassBlance
     (
        volScalarField& eta, // elevation
        const surfaceScalarField& q0, // bed load flux
        const volScalarField& E, // entrainment
        const volScalarField& D, // deposition
        const scalarField& triAreas,
        const scalar deltaT
     )const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
