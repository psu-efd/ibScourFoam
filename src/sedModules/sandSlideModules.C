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
#include "wallPolyPatch.H"

#include "SortableList.H"
// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * //



// diffusion-based
bool Foam::sedModules::sandSlide
(
    volScalarField& eta, // grain size
    const scalar& mus, // static friction coefficient
    const scalar& K0, // diffusivity coefficient
    const scalar& err1,
    const scalar& err2, 
    const scalar& timeStep,
    const scalar& maxIter
)const
{
    bool checkSandSlide=true;
    const fvMesh& mesh=eta.mesh();

    // change SolverPerformance to 0;
    simpleObjectRegistryEntry* objPtr = debug::debugObjects().lookupPtr("SolverPerformance");
    List<simpleRegIOobject*>& objects = *objPtr;
    forAll(objects,I)
    {
        IStringStream is("0");
        objects[I]->readData(is);
    }
    volScalarField::Boundary& etaPatches=eta.boundaryFieldRef();
    forAll(etaPatches,patchi)
    {
        etaPatches[patchi]=etaPatches[patchi].patchInternalField();
        etaPatches.set
        (
            patchi,
            fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchi], eta)
        );
    }
    volScalarField oldEta=eta*0;
    surfaceScalarField K
    (
        IOobject
        (
            "K",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "0",
            sqr(dimLength)/dimTime,
            K0
        )
    );


    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();

    scalar iterCounter = 0;
    scalar slopeCounter = 0;
    word  gradientScheme = "LGS+GNGS";

    while(checkSandSlide)
    {
     
        // ensure numerical stability
        forAll(eta.boundaryFieldRef(),patchi)
        {
            eta.boundaryFieldRef()[patchi]=eta.boundaryFieldRef()[patchi].patchInternalField();
        }
        
        slopeCounter = 0;
         
        if(gradientScheme.match("LGS"))
        {
            forAll(K.primitiveFieldRef(),I)
            {
		point owner(mesh.cellCentres()[faceOwner[I]]);
                point neigh(mesh.cellCentres()[faceNeighbour[I]]);
                vector length = owner -neigh;
                length.z() = 0;
                scalar dEta = eta.primitiveField()[faceOwner[I]]-eta.primitiveField()[faceNeighbour[I]];
                scalar tanPhi = mag(dEta/mag(length));
                
		if(tanPhi<mus )
                {
                    K[I]=0;
                }
                else
                {

                   K[I]=K0;
                   {
                      slopeCounter++;
                   }

                }
    
                
	    }
        }
        else if(gradientScheme.match("GNGS") or gradientScheme.match("LGS+GNGS"))
	{
             surfaceScalarField sGradEta = mag(linearInterpolate(fvc::grad(eta)));
             forAll(K.primitiveFieldRef(),I)
             {

                if(sGradEta[I]<mus )
                {
                    K[I]=0;
                }
                else
                {
                   K[I]=K0;
                   {
                      slopeCounter++;
                   }

                }

             }

        }
        else
        {
              FatalErrorIn("gradientScheme")
              << "The specified gradientScheme is not valid." << nl
              << abort(FatalError);

        }
        
        forAll(etaPatches,patchi)
        {
            etaPatches[patchi]=etaPatches[patchi].patchInternalField();
        } 
        volScalarField dEta = mag(eta-oldEta);
    
        if(gradientScheme.match("LGS+GNGS") && ((gMax(dEta)<err1) ))
        {
            gradientScheme = "LGS";
        }        
        

        if(returnReduce(slopeCounter, maxOp<label>())==0 || (gMax(dEta)<err2) )
        {
            break;
        }
     


        iterCounter++;
        
        oldEta = eta;
	
	    scalar rDeltaT = 1.0/timeStep;
        fvScalarMatrix DDT(eta,eta.dimensions()*dimVol/dimTime);
        DDT.diag()=  rDeltaT*mesh.V();
        DDT.source() = rDeltaT*oldEta.primitiveField()*mesh.V();

        fvScalarMatrix sandSlideEqn
        (
     
            DDT-fvm::laplacian(K,eta)
        );
        
        
       
        sandSlideEqn.solve();
        
        
    }
      


   Info<<"Iter. "<<iterCounter<<", total number of large slope cells: "<<returnReduce(slopeCounter, maxOp<label>())<<endl;
   Info<<"K0 = "<<K0<<endl;
               
    // change SolverPerformance to 1;
    forAll(objects,I)
    {
        IStringStream is("1");
        objects[I]->readData(is);
    }
    
    
    return checkSandSlide;
    
}

// geometric-based
bool Foam::sedModules::sandSlide
(
    volScalarField& eta, // grain size
    const scalar& mus, // static friction coefficient
    const scalar& alpha, // acceleration on mus
    const scalar& beta, // relaxation factor
    const scalar& maxIter
)const
{
    bool checkSandSlide=true;
    const fvMesh& mesh=eta.mesh();

    // change SolverPerformance to 0;
    simpleObjectRegistryEntry* objPtr = debug::debugObjects().lookupPtr("SolverPerformance");
    List<simpleRegIOobject*>& objects = *objPtr;
    forAll(objects,I)
    {
        IStringStream is("0");
        objects[I]->readData(is);
    }
    volScalarField::Boundary& etaPatches=eta.boundaryFieldRef();
    forAll(etaPatches,patchi)
    {
        etaPatches[patchi]=etaPatches[patchi].patchInternalField();
        etaPatches.set
        (
            patchi,
            fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchi], eta)
        );
    }
    volScalarField oldEta=eta;
    scalar counter=0;

    const labelListList& cellCells = mesh.cellCells();
    const pointField& C = mesh.C();
    const scalarField& Ah = mesh.V();

    while( counter < maxIter && checkSandSlide)
    {
        // ensure numerical stability
        forAll(eta.boundaryFieldRef(),patchi)
        {
            eta.boundaryFieldRef()[patchi]=eta.boundaryFieldRef()[patchi].patchInternalField();
        }

        checkSandSlide = false;
        label slopeCounter = 0;
        volVectorField gradEta = fvc::grad(eta);
        volScalarField magGradEta = mag(gradEta);
        SortableList<scalar> slopes (magGradEta.primitiveField());
        labelList newList (Ah.size());
        sortedOrder(slopes,newList);
       

        forAll(C,I0)
        {
            label I = newList[newList.size()-I0-1];
//           label I = newList[I0];
//           label I=I0;
    
            point C_p=C[I];
            C_p.z()=0.0;
            const labelList& myCells = cellCells[I];
            scalarField Spi(myCells.size(),0.0);
            scalarField Lpi(myCells.size(),0.0);
            scalar Ci=0.0;
            scalar Cp=Ah[I];
            forAll(myCells,II)
            {
                label neighID=myCells[II];
                point C_n=C[neighID];
                C_n.z()=0.0;
                Lpi[II]=mag(C_p-C_n);
                Spi[II]=(eta[I]-eta[neighID])/Lpi[II];
                if(mag(Spi[II])>mus)
                {
                    checkSandSlide=true;
                    slopeCounter++;
                }
                Spi[II]=max(-mus*alpha,Spi[II]);
                Spi[II]=min(mus*alpha,Spi[II]);
                Ci+=Ah[neighID]*(eta[I]-eta[neighID]-Lpi[II]*Spi[II]);
                Cp+=Ah[neighID];


            }
            scalar totalMass=0.0;
            scalar dEtap=0.0;
            scalar oldEdap=eta[I];
            if(mag(Cp)>SMALL)
            {
                dEtap=-Ci/Cp;
            }
            totalMass+=Ah[I]*dEtap;

            eta[I]+=dEtap*beta;

            forAll(myCells,II)
            {
                label neighID=myCells[II];
                scalar dEtan=oldEdap+dEtap-eta[neighID]-Lpi[II]*Spi[II];
                eta[neighID]+=dEtan*beta;
                totalMass+=Ah[neighID]*dEtan;
            }


            if(mag(totalMass/Ah[I])>0.000001)
            {
                Info<<totalMass/Ah[I]<<tab<<eta[I]<<tab<<dEtap;
                Info<<endl;
            }

        }

         if(returnReduce(slopeCounter, maxOp<label>())>0) checkSandSlide=true;
         counter++;
        if(!checkSandSlide)
        {
            if(debug)
            {
                Info<<"Iter. "<<counter<<", total number of large slope cells: "<<slopeCounter<<endl;
            }
            Info<<"Total iteration number of sandSlide: "<<counter<<endl;

            // change SolverPerformance to 1;
            forAll(objects,I)
            {
                IStringStream is("1");
                objects[I]->readData(is);
            }
            return checkSandSlide;
        }

        if(debug)
        {
            Info<<"Iter. "<<counter<<", total number of large slope cells: "<<slopeCounter<<endl;
        }
    }
    Info<<"Total iteration number of sandSlide: "<<counter<<endl;

    // change SolverPerformance to 1;
    forAll(objects,I)
    {
        IStringStream is("1");
        objects[I]->readData(is);
    }
    return checkSandSlide;

}
// ************************************************************************* //
