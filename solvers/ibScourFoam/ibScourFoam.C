/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    ibPisoFoam

Description
    Immersed boundary method solver for coupling between incompressible turbulent flow,
    and morphodynamics governed by Exner equation.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity
    - immersed boundary method with multiple objects
    - immersed boundary wall functions for high Reynold number
    - sand slide algorithms
    - local scour modelling

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "immersedBoundaryFvMesh.H"
#include "wallFvPatch.H"
#include "LduMatrix.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createImmersedBoundaryFvMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // Update IB information including:
    //      readDict(); markCells(); makeStencilsInfo();
    //      makeTriAddressing(); makeExtrudedMesh();
    // IB method is triggered here, and controlled by updateIB_
    // after updating IB information, updateIB_ will be false
    mesh.update();

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        const double Oldtime1 = runTime.elapsedCpuTime();
        mesh.update();

        const double Oldtime2 = runTime.elapsedCpuTime();
        
        // evaluateU has following functions:
        // 1. take care of dry/wet cells (only consider one object yet)
        // 2. set deadCell velocity to be zero
        // 3. setInlet (with velocity profile, only for sediment)
        // located around Line 584 in immersedBoundaryFvMeshUpdateCellValues.C 
        mesh.evaluateU();
        #include "immersedBoundaryCourantNo.H"

        // Pressure-velocity PISO corrector
        {
            Info<<"*** Starting Predictor Step"<<endl;
            #include "UEqn.H"

            Info<<"*** Starting Correct Step Loop"<<endl;
            // --- PISO loop
            while (piso.correct())
            {
                #include "pEqn.H"
            }

        }
        

        laminarTransport.correct();
        
      
        
        // Specialized RASModels are available, such as
        // kEpsilonIB, kOmegaIB, kOmegaSSTIB, kOmegaSSTLMIB, kOmegaSSTSASIB
        // IB wall function is implemented, but major part can be found in 
        // immersedBoundaryFvMesh::kOmegaCorrection() and 
        // immersedBoundaryFvMesh::kEpsilonCorrection()
        // which are located in immersedBoundaryFvMeshTurbulence.C
        // kOmegaCorrection() has detailed comments.
        
        
        turbulence->correct();
        
        #include "calcForce.H"
        
        Info<<"*** Starting Sediment Step"<<endl;
        const double Oldtime3 = runTime.elapsedCpuTime();
       

        
        // For scour modeling, this is the main program of solving Exner Eqn.
        // It is located in immersedBoundaryFvMeshPostEvaluation.C
        mesh.postEvaulation();
        
        // It is pre treatment for suspended sediment transport.
        // Not well implemented yet.
   //     mesh.evaluateC();

        // Compute suspended sediment transport
   //     #include "CEqn.H"

        // Write out integral pressure and viscous
   //     mesh.pressureOutput();
        runTime.write();

        const double Oldtime4 = runTime.elapsedCpuTime();

        if(mesh.debug)
        {
            Info<< "Total_Time = " << Oldtime4-Oldtime1 << " s" <<endl;
            Info<< "IB_Update_Time = " << Oldtime2-Oldtime1 << " s" <<endl;
            Info<< "IB_CFD_Time = " << Oldtime3-Oldtime2 << " s" <<endl;
            Info<< "Post_Time = " << Oldtime4-Oldtime3 << " s" <<endl;
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    // clear out all pointers
    mesh.finalClearOut();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
