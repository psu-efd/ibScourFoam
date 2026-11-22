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
    ibScourPimpleFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    around immersed solids (immersed boundary method), based on pimpleFoam.

    The immersed boundary (libibScour) blocks the faces between fluid and
    solid cells, so the flux is conservative and no mass crosses the solid
    surfaces. The velocity and turbulence quantities of the IB cells follow
    from a wall model (constant/immersedBoundaryDict).

    With a bed model (constant/immersedBoundaryDict, entry 'bed') the bed
    evolves by the Exner equation from the wall shear stress of the IB cells
    and the immersed boundary is re-classified when it has moved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "correctedSnGrad.H"
#include "profiling.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "immersedBoundary.H"
#include "ibWallModel.H"
#include "ibTurbulenceConstraint.H"
#include "ibSuspendedLoad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " around immersed solids (conservative sharp immersed boundary)."
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
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // Testing: with solveFlow false the flow fields are frozen and
    // only the classification, the fresh-cell initialisation and the bed
    // evolution run
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
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Re-classify the cells if a solid or the bed moved (no-op for
        // static solids), and initialise fresh cells
        ib.update();
        ib.initialiseFreshCells(U);

        // --- Pressure-velocity PIMPLE corrector loop
        while (solveFlow && pimple.loop())
        {
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
                laminarTransport.correct();
                turbulence->correct();
                if (nutPtr)
                {
                    ib.wallModel().correctNut(*nutPtr);
                }
            }
        }

        ib.checkFlux(phi, U);

        // Suspended load: concentration with the bed exchange of the IB
        // cells; the net deposition goes to the Exner equation below
        if (suspendedPtr)
        {
            if (!solveFlow)
            {
                // Frozen flow: wall-model values for the bed exchange
                ib.wallModel().correct(U);
            }
            if (ib.hasBed())
            {
                tmp<scalarField> tErodible(ib.bed().erodibleIbCells(ib));
                if (tErodible)
                {
                    suspendedPtr->setErodible(tErodible());
                }
            }
            suspendedPtr->correctDiffusivity
            (
                laminarTransport.nu()(),
                turbulence->nut()()
            );
            suspendedPtr->solve
            (
                phi,
                gAverage(laminarTransport.nu()().primitiveField())
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

        // Morphodynamics: Exner equation of the bed from the wall shear
        // stress of this step (re-classification at the next step)
        ib.evolveBed(runTime.deltaTValue());

        if
        (
            runTime.writeTime()
         && ib.dict().getOrDefault<bool>("writeMomentumTerms", false)
        )
        {
            #include "writeMomentumTerms.H"
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
