/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 ibScourFoam developers (Penn State)
-------------------------------------------------------------------------------
License
    This file is part of ibScourFoam, built on OpenFOAM.

    ibScourFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ibScourFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ibScourFoam.  If not, see <http://www.gnu.org/licenses/>.

Application
    ibScourInterFoam

Group
    grpMultiphaseSolvers

Description
    Transient solver for two incompressible, immiscible fluids (water and
    air, VoF interface capturing) around immersed solids and over an
    immersed erodible bed: interFoam with the immersed boundary of
    libibScour.

    The volume fraction is advected by MULES (default) or by isoAdvector
    (entry 'alphaScheme isoAdvector;' in the alpha solver dictionary of
    fvSolution). No flux crosses the blocked faces between fluid and solid
    cells, so no water enters the solids; the gravity and pressure fluxes
    are reconstructed from the open faces of every cell, which keeps a
    fluid at rest over an immersed bed at rest. The velocity and turbulence
    quantities of the IB cells follow from the wall model
    (constant/immersedBoundaryDict), as in ibScourPimpleFoam.

    Stage 1a: fixed solids and bed. The sediment transport with
    a wet/dry bed (stage 1b), the soil under the free surface (stage 2) and
    structures in two phases (stage 3) follow.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "correctedSnGrad.H"
#include "profiling.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "isoAdvection.H"
#include "immersedBoundary.H"
#include "ibWallModel.H"
#include "ibTurbulenceConstraint.H"
#include "ibSuspendedLoad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " (VoF) around immersed solids and over an immersed bed."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createImmersedBoundary.H"
    #include "createSuspendedLoad.H"
    #include "CourantNo.H"
    #include "alphaCourantNo.H"
    #include "setInitialDeltaT.H"

    const bool solveFlow = ib.dict().getOrDefault<bool>("solveFlow", true);
    if (!solveFlow)
    {
        Info<< "solveFlow false: the flow is not solved (bed-only test)"
            << nl << endl;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Re-classify the cells if a solid or the bed moved (no-op for
        // static solids), and initialise fresh cells
        ib.update();
        ib.initialiseFreshCells(U, alpha1, p_rgh);

        // --- Pressure-velocity PIMPLE corrector loop
        while (solveFlow && pimple.loop())
        {
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            // Wall-model values of the IB cells from the current velocity
            ib.wallModel().correct(U);

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                addProfiling(peqnProf, "pEqn");
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                addProfiling(turbProf, "turbulence::correct");
                turbulence->correct();
                if (nutPtr)
                {
                    ib.wallModel().correctNut(*nutPtr);
                }
            }
        }

        ib.checkFlux(phi, U);

        // Suspended load in the water only (alpha C with the water flux);
        // the net deposition of the wet IB cells goes to the Exner step
        if (suspendedPtr)
        {
            if (!solveFlow)
            {
                ib.wallModel().correct(U);
            }
            if (nutPtr)
            {
                nutField = *nutPtr;
            }
            if (ib.hasBed())
            {
                tmp<scalarField> tErodible(ib.bed().erodibleIbCells(ib));
                if (tErodible)
                {
                    suspendedPtr->setErodible(tErodible());
                }
            }
            suspendedPtr->correctDiffusivity(mixture.nu()(), nutField);
            suspendedPtr->solve
            (
                useIsoAdvector ? advectorPtr->alphaPhi() : alphaPhi10,
                nuWater,
                alpha1
            );
            if (ib.hasBed())
            {
                ib.bed().setSuspendedExchange
                (
                    ib,
                    suspendedPtr->netDeposition()
                );
            }
        }

        // Morphodynamics: wet IB cells only drive the bed (wetThreshold)
        ib.evolveBed(runTime.deltaTValue());

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
