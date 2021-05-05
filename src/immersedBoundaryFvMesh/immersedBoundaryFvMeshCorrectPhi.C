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
    immersedBoundaryFvMeshCorrectPhi.C
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvMesh.H"
#include "MeshWave.H"
#include "SortableList.H"
#include "fvCFD.H"
#include "turbulenceModel.H"
#include "wallFuncModules.H"
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void immersedBoundaryFvMesh::adjustDivPhiHbyA
(
    volScalarField& divPhiHbyA
) const
{
    forAll(objectsList(),objectID)
    {
        adjustDivPhiHbyA(divPhiHbyA,objectID);
    }
}

void immersedBoundaryFvMesh::adjustDivPhiHbyA
(
    volScalarField& divPhiHbyA,
    const label& objectID
) const
{
    const labelList& ibc = ibCellsList()[objectID];
    const pointField& ibp = ibHitPointsList()[objectID];
    forAll(ibc,I)
    {
        label cellID = ibc[I];
        scalar yIB = mag(ibp[I]-this->C()[cellID]);
        scalar ratio = 0.5+yIB/cellSize(cellID);
        divPhiHbyA[cellID]/=min(pow(ratio,2.0),1.0);
    }
}


//- Adjust the fluxes on immersed boundary
void immersedBoundaryFvMesh::ibFacesAdjust
(
    surfaceScalarField& phi,
    const volVectorField& U
) const
{
    const double Oldtime1=time().elapsedCpuTime();
    forAll(objectsList(),objectID)
    {
        ibFacesAdjust(phi,U,objectID);
    }
    const double Oldtime2=time().elapsedCpuTime();
    if(debug>1)
    {
        Info<<"ibFacesAdjust Executation Time = "<<Oldtime2-Oldtime1<< " s"<<endl;
    }
}

//- Adjust the fluxes on immersed boundary
void immersedBoundaryFvMesh::ibFacesAdjust
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const label& objectID
) const
{
    dictionary turbDict (this->lookupObject<dictionary>("turbulenceProperties"));
    word modelType=turbDict.lookup("simulationType");
    if(modelType=="laminar")
    {
        fromIbFaceSPointReconstruction(phi,U,objectID);
        return;
    }
    volScalarField& k = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("k"));

    const pointField& ibFaceHitPoints = ibFaceHitPointsList()[objectID];
    const pointField& ibFaceSamplingPoints = ibFaceSamplingPointsList()[objectID];
    const labelList& ibFaces = ibFacesList()[objectID];

    const pointField ibf_centers(this->faceCentres(),ibFaces);
    vectorField ibNormals = ibFaceHitPoints-ibFaceSamplingPoints;
    ibNormals=ibNormals/mag(ibNormals);

    vectorField USample(ibFaceSamplingPointsValues(U,objectID));
    // remove normal velocity
    vectorField tmpUsample=USample-(USample&ibNormals)*ibNormals;
    scalarField USampleMag=mag(tmpUsample);
    // k at sampling point
    scalarField kSample(ibFaceSamplingPointsValues(k,objectID)+SMALL);
    // IB distance
    scalarField yIB(mag(ibf_centers-ibFaceHitPoints)+SMALL);
    // sampling distance
    scalarField ySample(mag(ibFaceSamplingPoints-ibFaceHitPoints)+SMALL);

    word wallFunc_ = ibProperties().lookupOrDefault<word>("wallFunction","vanDriest");
    wallFunc_ = objectDictList()[objectID].lookupOrDefault<word>("wallFunction",wallFunc_);
    wallFuncModules wFM(wallFunc_,"Baykal2015");
    scalar Ks_ = ibProperties().lookupOrDefault<scalar>("roughnessHeight",0);
    Ks_ = objectDictList()[objectID].lookupOrDefault<scalar>("roughnessHeight",Ks_);

    const dictionary& transportProperties =
          this->lookupObject<dictionary>
          (
             "transportProperties"
          );
    dimensionedScalar nu(transportProperties.lookup("nu"));

    scalar nuLam=nu.value();

    forAll(ibFaces,I)
    {
        label faceID = ibFaces[I];
        vector UT(pTraits<vector>::zero);

        if(USampleMag[I]<sqrt(SMALL) or ySample[I]<sqrt(SMALL) or yIB[I]<1e-5)
        {
            UT*=0;
        }
        else
        {
            scalar uTau=wFM.uTau
                (
                    USampleMag[I],
                    kSample[I],
                    ySample[I],
                    Ks_
                );

            scalar yPlus=yIB[I]*uTau/nuLam;

            scalar ksp=Ks_*uTau/nuLam;
            ksp=max(1.0,ksp);

            if(mag(tmpUsample[I])>0)
            {
                UT=uTau*wFM.uPlus(yPlus,ksp)*tmpUsample[I]/mag(tmpUsample[I]);
            }
            else
            {
                UT*=0;
            }

        }
        vector Unormal = (USample[I]&ibNormals[I])*ibNormals[I];
        Unormal = Unormal*(yIB[I]/ySample[I]);
        UT=UT+Unormal;
        if(faceID<ibFaces.size())
        {
            phi[faceID]=UT&this->faceAreas()[faceID];
        }
    }
}

//- Adjust the fluxes on immersed boundary to obey continuity.
void immersedBoundaryFvMesh::ibCorrectPhi
(
    surfaceScalarField& phi
) const
{
    forAll(objectsList(),objectID)
    {
        phi = phi*sGammaList()[objectID];
    }
}

//- Adjust the fluxes on immersed boundary to obey continuity.
void immersedBoundaryFvMesh::immersedBoundaryAdjustPhi
(
    surfaceScalarField& phi,
    const volVectorField& U
) const
{
    forAll(objectsList(),objectID)
    {

        if(IBtypeList()[objectID]=="classic")
        {
            immersedBoundaryAdjustPhi(phi,U,objectID);
        }
    }
}

//- Adjust the fluxes on immersed boundary to obey continuity.
void immersedBoundaryFvMesh::immersedBoundaryAdjustPhi
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const label& objectID
) const
{

    const fvMesh& mesh = phi.mesh();

    const triSurface& surf(objectsList()[objectID]);
    // Sum the flux through the immersed boundary
    vectorField UrefValue(surf.size(),pTraits<vector>::zero);
    const scalar triFlux = 0;
    scalarField& phiInternal = phi.primitiveFieldRef();
    scalar fluxIn = 0;
    scalar fluxOut = 0;
    scalar fixedFlux = 0;

    // Get all IB faces
    const labelList& ibFaces = ibFacesList()[objectID];
    const boolList& ibFaceFlips = ibFaceFlipsList()[objectID];
    forAll (ibFaces, faceI)
    {
        const label curFace = ibFaces[faceI];
        const bool curFlip = ibFaceFlips[faceI];

        if (mesh.isInternalFace(curFace))
        {
            const scalar curFlux = phiInternal[curFace];

            if (!curFlip)
            {
                // Face points out of the live cell
                if (curFlux >= 0)
                {
                    // Flux out of the live cell
                    fluxOut += curFlux;
                }
                else
                {
                    // Flux into the live cell
                    fluxIn -= curFlux;
                }
            }
            else
            {
                // Face points into the live cell: flip it
                if (curFlux >= 0)
                {
                    // Flux into of the live cell
                    fluxIn += curFlux;
                }
                else
                {
                    // Flux out the live cell
                    fluxOut -= curFlux;
                }
            }
        }
           else
        {
            const label patchID =
                mesh.boundaryMesh().whichPatch(curFace);
            const label faceID =
                mesh.boundaryMesh()[patchID].whichFace(curFace);

            const scalar curFlux =
                phi.boundaryFieldRef()[patchID][faceID];

            // Note: only coupled patches may carry flux
            // In order to avoid double summation and
            // inconsistencies in the correction,
            // coupled face fluxes will NOT be corrected,
            // but only accounted for in the summation.
            if (mesh.boundaryMesh()[patchID].coupled())
            {
                // Only do the master side; slave will
                // be handled on the other side of the couple
                if (!curFlip)
                {
                    fixedFlux += curFlux;
                }
            }
        }
    }
    reduce(fluxIn, sumOp<scalar>());
    reduce(fluxOut, sumOp<scalar>());
    reduce(fixedFlux, sumOp<scalar>());
    scalar imbalance = (fluxIn - fluxOut + fixedFlux) - triFlux;

    if (debug)
    {
        Info<< "triFlux = " << triFlux
            << " fluxIn = " << fluxIn << " fluxOut = " << fluxOut
            << " fixedFlux = " << fixedFlux
            << " imbalance = " << imbalance
            << endl;
    }
    scalar massCorr = 1.0;

    if (mag(imbalance) > SMALL)
    {
        // Scaling required: scale to match the smaller of two
        // fluxes
        if (fluxIn > fluxOut)
        {
            // Scale down incoming flux
            // Note change of sign: imbalance is negative
            massCorr = 1 - imbalance/(fluxIn + SMALL);

            if (debug)
            {
                Info<< "Scaling down incoming flux with factor = "
                    << massCorr << endl;
            }

            scalar newFluxIn = 0;

            // Visit all incoming flux faces and re-scale the flux
            forAll (ibFaces, faceI)
            {
                const label curFace = ibFaces[faceI];
                const bool curFlip = ibFaceFlips[faceI];

                if (mesh.isInternalFace(curFace))
                {
                    // Take reference to current flux
                    scalar& curFlux = phiInternal[curFace];

                    if (!curFlip)
                    {
                        // Face points out of the live cell
                        if (curFlux < 0)
                        {
                            // Flux out of the live cell
                            curFlux *= massCorr;
                            newFluxIn -= curFlux;
                        }
                    }
                    else
                    {
                        // Face points into the live cell: flip it
                        if (curFlux >= 0)
                        {
                            // Flux out the live cell
                            curFlux *= massCorr;
                            newFluxIn += curFlux;
                        }
                    }
                }
              
            }
        }
        else
        {
            // Scale down outgoing flux
            massCorr = 1 + imbalance/(fluxOut + SMALL);

            if (debug)
            {
                Info<< "Scaling down outgoing flux with factor = "
                    << massCorr << endl;
            }

            scalar newFluxOut = 0;

            // Visit all outgoing flux faces and re-scale the flux
            forAll (ibFaces, faceI)
            {
                const label curFace = ibFaces[faceI];
                const bool curFlip = ibFaceFlips[faceI];

                if (mesh.isInternalFace(curFace))
                {
                    // Take reference to current flux
                    scalar& curFlux = phiInternal[curFace];

                    if (!curFlip)
                    {
                        // Face points out of the live cell
                        if (curFlux >= 0)
                        {
                            // Flux out of the live cell
                            curFlux *= massCorr;
                            newFluxOut += curFlux;
                        }
                    }
                    else
                    {
                        // Face points into the live cell: flip it
                        if (curFlux < 0)
                        {
                            // Flux out the live cell
                            curFlux *= massCorr;
                            newFluxOut -= curFlux;
                        }
                    }
                }
            
            }
        }
    }
}
// ************************************************************************* //



