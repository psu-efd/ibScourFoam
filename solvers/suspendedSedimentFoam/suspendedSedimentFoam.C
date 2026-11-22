/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Xiaofeng Liu, Penn State University
-------------------------------------------------------------------------------
    Derived from pimpleFoam of OpenFOAM:
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of ibScourFoam, built on OpenFOAM.

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
    suspendedSedimentFoam

Description
    Transient solver for incompressible, turbulent flow (pimpleFoam, static
    or moving mesh) with the transport of suspended sediment of one grain
    size (volume concentration C), one-way coupled: the sediment does not
    change the density or the flow.

    The mesh may move (dynamicMotionSolverFvMesh, as in pimpleFoam): with
    the deformingBed function object of libibScour the bed patch follows
    the Exner equation fed by the bed load and by E - D of the
    sedimentBedFlux condition, i.e. the body-fitted scour of Liu & Garcia
    (2008) with bed load and suspended load. The C
    equation uses the flux relative to the mesh motion.

    The sediment moves with the drift flux phi + ws.Sf (settling velocity
    ws along g from a settlingVelocityModel) and diffuses with
    DC = nu/ScC + nut/sigmaC (CEqn.H). Erosion and deposition at the bed
    come from the patch condition sedimentBedFlux of libsuspendedSediment
    (modes referenceLevel, HRN, LRN, Robin); the function object
    sedimentBedExchange of the same library writes E - D per bed face and
    the sediment budget.

    The C equation is solved once per time step after the PIMPLE loop.
    With solveFlow no the flow fields are frozen (only C is solved), with
    solveSediment no only the flow is solved (precursor run).

    Derived from pimpleFoam; the settling enters as the drift flux
    phiC = phi + ws.Sf of the C equation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CorrectPhi.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "uniformDimensionedFields.H"
#include "settlingVelocityModel.H"
#include "sedimentBedFluxFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow with suspended"
        " sediment (settling, erosion and deposition at the bed)."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "createSedimentFields.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (solveFlow && pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Absolute flux from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        // --- Suspended sediment (one-way coupled)
        #include "CEqn.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
