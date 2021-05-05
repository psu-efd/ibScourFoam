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
    immersedBoundaryFvMeshPressureCorrector.C
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvMesh.H"
#include "MeshWave.H"
#include "SortableList.H"
#include "fvCFD.H"
#include "turbulenceModel.H"
#include "wallFuncModules.H"
#include "globalMeshData.H"

//#include "LUinvert.H"
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::immersedBoundaryFvMesh::evaluateGradP
(
    volVectorField& psi
) const
{
    forAll(objectsList(), objectID)
    {
        
        if(IBtypeList()[objectID]=="classic")
        {

            forAll(ibCellsList()[objectID],I)
            {

                vector& gradPRef=psi[ibCellsList()[objectID][I]];
                gradPRef *=0;

            }


            const labelList& ibDeadCells(ibDeadCellsList()[objectID]);
            forAll(ibDeadCells,I)
            {
                vector& gradPRef=psi[ibDeadCells[I]];
                gradPRef= gradPRef*0.0;
            }

        }
        else if(IBtypeList()[objectID]=="mix")
        {
            const labelList& dce_org = ibDeadCellsList()[objectID]; // dead cell include ghost cells
            labelHashSet dceSet;
            labelHashSet ghostCellsSet(ghostCellsList()[objectID]);
            forAll(dce_org,I)
            {
                label cellID=dce_org[I];
                if(!ghostCellsSet.found(cellID) and !dceSet.found(cellID))
                {
                    dceSet.insert(cellID);
                }
            }
            labelList ibDeadCells(dceSet.toc());
            forAll(ibDeadCells,I)
            {
                label cellID=ibDeadCells[I];
                psi[cellID]=psi[cellID]*0;
            }

        }
    }
}


void Foam::immersedBoundaryFvMesh::evaluateP() const
{

    volScalarField& p = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("p"));
        
        
        forAll(objectsList(), objectID)
    {
        if(IBtypeList()[objectID]=="classic")
        {

           
            const labelList& ibDeadCells(ibDeadCellsList()[objectID]);
            forAll(ibDeadCells,I)
            {
                 scalar& pRef=p[ibDeadCells[I]];
                pRef *=0;
            }

        }
       
    }

}

void Foam::immersedBoundaryFvMesh::evaluatek() const
{

    volScalarField& k = const_cast<volScalarField&>
        (this->lookupObject<volScalarField>("k"));
        
    
    
    forAll(objectsList(),objectID)
    {
            const vectorField& ibWallShearStress = this->wallShearStress(objectID); 
            forAll(ibCellsList()[objectID],I)
            {
                
                k[ibCellsList()[objectID][I]] = mag(ibWallShearStress[I]) /sqrt(0.09);
            }
            const labelList& ibDeadCells(ibDeadCellsList()[objectID]);
            forAll(ibDeadCells,I)
            {
               k[ibDeadCells[I]] = 0 ;
                
            };
    }
}




// ************************************************************************* //



