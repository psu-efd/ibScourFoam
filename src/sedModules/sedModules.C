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

#include "sedModules.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "vectorTools.H"
#include "simpleObjectRegistry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(sedModules, 0);
	//addToRunTimeSelectionTable(searchableSurface, sedModules, dict);
	//word sedModules::meshSubDir = "triSurface";
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

	Foam::sedModules::sedModules
(
 const word& name
 )
	:
		namePtr_(new word(name))
{
	Info<<"Use "<<name<<" sediment transport models"<<endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //



Foam::sedModules::~sedModules()
{
	clearOut();
}


void Foam::sedModules::clearOut()
{
	deleteDemandDrivenData(namePtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::sedModules::bedLoadUpdateValues
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
 scalarField& phi, // direction between sN_dir and steepest bed slope (to be updated)
 scalarField& sN, // magnitude of shields number
 scalarField& sNC, // critical shields number (to be updated)
 vectorField& q0 // bed load flux (to be updated)
 )const
{
	if(name()=="Roulund2005")
	{
		//roulund2005(D50,ssG,mus,mud,triNormals,wss,Vs,Cb,sNC0,gravity,Zb,beta,phi,sN,sNC,q0,E,D);
		  roulund2005(D50,ssG,mus,mud,triNormals,wss,sNC0,gravity,beta,phi,sN,sNC,q0);
	}
	else
	{
		Info<<"Cannot find sedment modules of "<<name()<<endl;
	}
}

void Foam::sedModules::bedLoadUpdateValues
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
 )const
{
	
	bedLoadUpdateValues(D50,ssG,mus,mud,triNormals.primitiveField(),wss.primitiveField(),
			sNC0,gravity,beta.primitiveFieldRef(),phi.primitiveFieldRef(),
			sN.primitiveFieldRef(),sNC.primitiveFieldRef(),q0.primitiveFieldRef()
			);
	forAll(wss.boundaryField(),patchI)
	{
		bedLoadUpdateValues(D50,ssG,mus,mud,triNormals.boundaryField()[patchI],wss.boundaryField()[patchI],
				sNC0,gravity,beta.boundaryFieldRef()[patchI],phi.boundaryFieldRef()[patchI],
				sN.boundaryFieldRef()[patchI],sNC.boundaryFieldRef()[patchI],q0.boundaryFieldRef()[patchI]
			    );       
	}
}

void Foam::sedModules::bedLoadUpdateValues
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

 )const
{
	bedLoadUpdateValues(D50,ssG,mus,mud,triNormals.primitiveField(),wss.primitiveField(),
			sNC0,gravity,beta.primitiveFieldRef(),phi.primitiveFieldRef(),
			sN.primitiveFieldRef(),sNC.primitiveFieldRef(),q0.primitiveFieldRef()
			);
	forAll(wss.boundaryField(),patchI)
	{
		bedLoadUpdateValues(D50,ssG,mus,mud,triNormals.boundaryField()[patchI],wss.boundaryField()[patchI],
				sNC0,gravity,beta.boundaryFieldRef()[patchI],phi.boundaryFieldRef()[patchI],
				sN.boundaryFieldRef()[patchI],sNC.boundaryFieldRef()[patchI],q0.boundaryFieldRef()[patchI]
				);       
	}
}

void Foam::sedModules::ibEDUpdateValues
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
	if(name()=="Roulund2005")
	{
		EDUpdateValues(D50,ssG,mus,mud,fN,wss,Vs,Cb,sNC0,gravity,Zb,E,D);
	}
	else
	{
		Info<<"Cannot find sedment modules of "<<name()<<endl;
	}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "bedLoadModules.C"
#include "massBalanceModules.C"
#include "sandSlideModules.C"
// ************************************************************************* //
