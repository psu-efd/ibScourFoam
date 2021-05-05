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


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * //
void Foam::sedModules::fixWallFlux
(
       surfaceScalarField& qFlux // bed load flux
)const
{
    const fvMesh& mesh = qFlux.mesh();
    surfaceScalarField::Boundary& qFluxPatch=qFlux.boundaryFieldRef();
    forAll(mesh.boundaryMesh(),patchI)
    {
        if (isA<wallPolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            scalar fluxIn = 0;
            scalar fluxOut = 0;
            forAll(qFluxPatch[patchI],I)
            {
                if(qFluxPatch[patchI][I]>0)
                {
                    fluxIn += qFluxPatch[patchI][I];
                }
                else
                {
                    fluxOut += qFluxPatch[patchI][I];
                }
            }
            scalar fluxCor = 1.0;
            scalar sumFlux=fluxIn+fluxOut;
            if(sumFlux>0)
            {
                fluxCor = 1.0-sumFlux/fluxIn;
                fluxIn = 0.0;
                fluxOut = 0.0;
                forAll(qFluxPatch[patchI],I)
                {
                    if(qFluxPatch[patchI][I]>0)
                    {
                        qFluxPatch[patchI][I] *= fluxCor;
                        fluxIn += qFluxPatch[patchI][I];
                    }
                    else
                    {
                        fluxOut += qFluxPatch[patchI][I];
                    }
                }
            }
            else if(sumFlux<0)
            {
                fluxCor = 1.0-sumFlux/fluxOut;
                fluxIn = 0.0;
                fluxOut = 0.0;

                forAll(qFluxPatch[patchI],I)
                {
                    if(qFluxPatch[patchI][I]<0)
                    {
                        qFluxPatch[patchI][I] *= fluxCor;
                        fluxIn += qFluxPatch[patchI][I];
                    }
                    else
                    {
                        fluxOut += qFluxPatch[patchI][I];
                    }
                }
            }

            //if(debug)
            {
                Info<<"Scaling bed-load flux on BC of "<<mesh.boundaryMesh()[patchI].name()
                    <<" with factor "<<fluxCor<<endl;
            }
        }
    }
}

void Foam::sedModules::feedInlet
(
       surfaceScalarField& qFlux // bed load flux
)const
{
    const fvMesh& mesh = qFlux.mesh();
    volScalarField divQ=fvc::div(qFlux);
    surfaceScalarField::Boundary& qFluxPatch=qFlux.boundaryFieldRef();

    label patchID_inlet = mesh.boundaryMesh().findPatchID("inlet");
    label patchID_outlet = mesh.boundaryMesh().findPatchID("outlet");
    scalar fluxOut=0.0;
    scalar fluxIn=0.0;
    scalar fluxCor = 1.0;
    if(patchID_inlet<0 or patchID_outlet<0) return;
    // outlet BC
    {
        const labelList& faceCells=mesh.boundaryMesh()[patchID_outlet].faceCells();

        forAll(faceCells,I)
        {
            scalar& BCflux = qFluxPatch[patchID_outlet][I];
            fluxOut +=BCflux;
        }
    }

    // inlet BC
    {
        const labelList& faceCells=mesh.boundaryMesh()[patchID_inlet].faceCells();
        forAll(faceCells,I)
        {
            scalar& BCflux = qFluxPatch[patchID_inlet][I];
            fluxIn+=BCflux;
        }
    }
    scalar sumFlux=fluxIn+fluxOut;

    // if (sumFlux>0) needs more feed in
    // if (sumFlux<0) needs less feed in
    if(fluxIn!=0)
    {
        fluxCor = fluxCor-sumFlux/fluxIn;
        qFluxPatch[patchID_inlet]*=fluxCor;
    }
    else
    {
        qFluxPatch[patchID_inlet] = fluxOut*qFluxPatch[patchID_outlet].patch().magSf()
                /sum(qFluxPatch[patchID_outlet].patch().magSf());
    }
}

void Foam::sedModules::checkMassBlance
(
    volScalarField& eta, // elevation
    const surfaceScalarField& qFlux, // bed load flux
    const volScalarField& E, // entrainment
    const volScalarField& D, // deposition
    const scalarField& triAreas,
    const scalar deltaT
)const
{
    if(Pstream::myProcNo()!=Pstream::masterNo()) return;
    scalarField H=eta.mesh().V()/triAreas;
    scalarField oldEta = eta.oldTime().primitiveField();
    scalar sumOldEta = sum(oldEta*triAreas);
    scalar sumEta = sum(eta*triAreas);
    scalar sumE = sum(E*triAreas);
    scalar sumD = sum(D*triAreas);
    const surfaceScalarField::Boundary& qFluxPatch=qFlux.boundaryField();
    scalar sumQFlux = 0.0; //outward

    forAll(qFluxPatch,patchI)
    {
        forAll(qFluxPatch[patchI],I)
        {
            label cellID=qFluxPatch[patchI].patch().faceCells()[I];
            sumQFlux += qFluxPatch[patchI][I]/H[cellID];
        }
    }
    scalar deltaEta = (sumEta-sumOldEta)/deltaT;
    scalar deltaFlux = -sumQFlux+sumD-sumE;
    Info<<"Checking bed-load flux balance: "<<((deltaEta-deltaFlux)/(deltaEta+SMALL)*100.0)<<"%"<<endl;
}
// ************************************************************************* //
