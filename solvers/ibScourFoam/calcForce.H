//calculate/ force magnitude


const fvPatchList& patches = mesh.boundary();



forAll(patches, patchi)
{
    const fvPatch& currPatch = patches[patchi];

    if (isA<wallFvPatch>(currPatch))
    {
        
        const volSymmTensorField Reff = turbulence->devReff();
        const vectorField& Sfp = mesh.Sf().boundaryField()[patchi];

        const symmTensorField& Reffp = Reff.boundaryField()[patchi];

        vectorField wallShearStress(Sfp & Reffp);
            
        vector  wallShear(pTraits<vector>::zero);
        vector  pressure(pTraits<vector>::zero);
        vector  totalWallShear(pTraits<vector>::zero);
        vector  totalPressure(pTraits<vector>::zero);
        
        forAll(currPatch, faceI)
        {
            label faceCellID = currPatch.faceCells()[faceI];
            if(gamma[faceCellID] > 0)
            {
                wallShear = wallShear + 1000*wallShearStress[faceI];
                pressure = pressure + 1000*(p.boundaryFieldRef()[patchi][faceI])*Sfp[faceI];
                
                
            }
        }
        
        totalWallShear = returnReduce(wallShear, maxOp<vector>());
        totalPressure = returnReduce(pressure, maxOp<vector>());
        
        Info<<"The friction on the object is    "<< totalWallShear <<endl;
        Info<<"The form drag on the object is    "<< totalPressure <<endl;
        
    }
}
