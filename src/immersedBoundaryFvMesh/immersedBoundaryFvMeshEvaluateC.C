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
    immersedBoundaryFvMeshEvaluateC.C
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvMesh.H"
#include "MeshWave.H"
#include "SortableList.H"
#include "fvCFD.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::immersedBoundaryFvMesh::evaluateC() const
{

    volScalarField& C = const_cast<volScalarField&>
    (this->lookupObject<volScalarField>("C"));

    volVectorField& Vs = const_cast<volVectorField&>
    (this->lookupObject<volVectorField>("Vs"));

    surfaceScalarField& phi = const_cast<surfaceScalarField&>
    (this->lookupObject<surfaceScalarField>("phi"));

    volVectorField& U = const_cast<volVectorField&>
    (this->lookupObject<volVectorField>("U"));

    //ibCorrectPhi(phi);
    phi=fvc::flux(U);
    forAll(objectsList(),objectID)
    {
        BCtype_ = objectDictList()[objectID].lookupOrDefault<word>("BCtype","fixedWall");

        if(IBtypeList()[objectID]=="classic" and objectDictList()[objectID].found("sediment"))
        {


            forAll(entrainment(objectID),I)
            {
                label cellID = ibCellsList()[objectID][I];
                C[cellID] = entrainment(objectID)[I]/mag(Vs[cellID]);
            }

            const labelList& ibDeadCells(ibDeadCellsList()[objectID]);
            forAll(ibDeadCells,I)
            {
                label cellID=ibDeadCells[I];
                C[cellID]=0;
            }
            labelHashSet ibDeadCellsSet(ibDeadCells);
            // Evaluate the coupled patchField, to replace evaluateCoupled()
            volScalarField::Boundary& Upatches = C.boundaryFieldRef();
            forAll (Upatches, patchI)
            {
                if (!Upatches[patchI].coupled())
                {
                    forAll(Upatches[patchI],I)
                    {
                        label faceCellID= this->boundary()[patchI].faceCells()[I];
                        if(ibDeadCellsSet.found(faceCellID))
                        {
                            C.boundaryFieldRef()[patchI][I]=C[faceCellID];
                            Vs.boundaryFieldRef()[patchI][I]=Vs[faceCellID];
                        }
                    }

                }
            }
        }
        else if(IBtypeList()[objectID]=="mix" and objectDictList()[objectID].found("sediment"))
        {

        }
    }
}


template<class Type>
void Foam::immersedBoundaryFvMesh::fromSPointReconstructionNeumann
(
    GeometricField<Type,fvPatchField,volMesh>& psi,
    label objectID
) const
{
    GeometricField<Type,fvPatchField,volMesh>& U=psi;
    Field<Type>& UI = U.primitiveFieldRef();
    const pointField& samplingPoints = samplingPointsList()[objectID];
    const labelList& ibCells = ibCellsList()[objectID];

    label counter = 0;
    scalarField error(samplingPoints.size(), 0);
    scalar maxError = 0;

    do
    {
        counter++;

        scalar maxMagUI = 0;

        Field<Type> samplingPointsValues
                (
                    stencilInterpolation
                    (
                        UI,
                        objectID,
                        ibCells,
                        samplingStencils(objectID).cellCells(),
                        samplingStencils(objectID).cellProcCells(),
                        samplingStencils(objectID).procCells(),
                        samplingStencils(objectID).procCentres(),
                        samplingStencils(objectID).weights(),
                        samplingStencils(objectID).procWeights()
                    )
                );


        forAll(ibCells, I)
        {
            label cellID = ibCells[I];
            Type oldUI = UI[cellID];
            Type samplingValue = samplingPointsValues[I];
            UI[cellID] =samplingValue;

            error[I] = mag(UI[cellID] - oldUI);
            if(maxMagUI<mag(UI[cellID])) maxMagUI = mag(UI[cellID]);
        }

        reduce(maxMagUI, maxOp<scalar>());

        error /= maxMagUI + SMALL;

        maxError = gMax(error);


        if(debug)
        {
            Info<<"Reconstruct "<<psi.name()<<", maxError: "<<maxError<<" counter: "<<counter<<endl;
        }
    }
    while (maxError > 0.001 && counter < 10);


}

void Foam::immersedBoundaryFvMesh::updateVs( const label& objectID) const
{
    // read information from dictionary
    const dictionary& dict = objectDictList()[objectID].subDict("sediment");
    vector gravity(dict.lookup("gravity"));


    scalar d50(readScalar(dict.lookup("grainDiameter")));
    scalar ssG(dict.lookupOrDefault<scalar>("specificGravity",2.69));

    word VsModel(dict.lookupOrDefault<word>("VsModel","constant"));
    vector V0(dict.lookup("V0"));


    volVectorField& Vs = const_cast<volVectorField&>
    (this->lookupObject<volVectorField>("Vs"));

    const dictionary& transportProperties =
          this->lookupObject<dictionary>
          (
             "transportProperties"
          );
    dimensionedScalar nu(transportProperties.lookup("nu"));

    if (VsModel == "constant")  //constant settling velocity
    {
        Vs.primitiveFieldRef() = V0;
    }
    else if (VsModel == "simple") //hindered
    {
        V0 = gravity*sqr(d50)/(18.0*nu.value())*(ssG);
        Vs.primitiveFieldRef() = V0;
    }
    else if (VsModel == "vanRijn1993") //vanRijn1993
    {
        if(d50>1e-6 and d50<1e-4)
        {
            V0 = gravity*sqr(d50)/(18.0*nu.value())*(ssG);
        }
        else if(d50<1e-3)
        {
            scalar temp1=1+(0.01*(ssG)*mag(gravity)*d50*d50*d50)/sqr(nu.value());
            V0 = gravity/mag(gravity)*10*nu.value()/d50*(Foam::sqrt(temp1)-1);
        }
        else
        {
            V0 = gravity/mag(gravity)*1.1*Foam::sqrt((ssG)*mag(gravity)*d50);
        }
        Vs.primitiveFieldRef() = V0;
    }
    else if (VsModel == "Fredse1992") //Fredse and Deigaard 1992 pp.200
    {
        scalar A1=1.4;
        scalar A2=36*nu.value()/d50;
        scalar A3=4*(ssG-1)*mag(gravity)*d50/3.0;

        V0 = gravity/mag(gravity)*Foam::sqrt(A3/A1+A2*A2/4.0/A1/A1);
        Vs.primitiveFieldRef() = V0;
    }
    else
    {
        FatalErrorIn("Foam::immersedBoundaryFvMesh::updateVs()")
            << "Unknown VsModel : " << VsModel
            << abort(FatalError);
    }

    Vs.correctBoundaryConditions();
}

// ************************************************************************* //

