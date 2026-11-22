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

Description
    Pressure corrector of interFoam with blocked immersed-boundary faces:
    no flux (volumetric, gravity or pressure) crosses a face between fluid
    and solid cells, the solid cells are pinned and the velocity is
    reconstructed from the fluxes of the open faces.

\*---------------------------------------------------------------------------*/

{
    volScalarField rAU("rAU", 1.0/UEqn.A());

    // rAU on the faces with the blocked faces removed
    const surfaceScalarField rAUf(ib.maskedLaplacianCoeff(rAU));

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi))
    );
    MRF.makeRelative(phiHbyA);
    ib.maskFlux(phiHbyA);

    if (p_rgh.needReference())
    {
        adjustPhi(phiHbyA, U, p_rgh);
    }

    // Gravity and surface-tension flux on the open faces only
    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // phig is built with the laplacian coefficient, which stays finite on
    // the faces between two solid cells, so a density jump inside a solid
    // (the free surface through a cylinder) puts a buoyancy flux on faces
    // that carry no flow; it would survive the subset pressure solve and
    // throttle the time step through the Courant number
    ib.maskFlux(phiHbyA);

    // Update the pressure BCs to ensure flux consistency; the unmasked
    // rAU is used because the boundary faces of solid cells have
    // rAUf = 0 (their flux and velocity are zero as well)
    constrainPressure(p_rgh, U, phiHbyA, rAU, MRF);

    if (pimple.nCorrPISO() <= 1)
    {
        tUEqn.clear();
    }

    if (ib.pressureSubset())
    {
        ib.solvePressure
        (
            p_rgh,
            rAUf,
            phiHbyA,
            phi,
            pimple,
            pRefCell,
            getRefCellValue(p_rgh, pRefCell)
        );
    }
    else
    {
        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix p_rghEqn
            (
                fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
            );

            p_rghEqn.setReference
            (
                pRefCell,
                getRefCellValue(p_rgh, pRefCell)
            );
            ib.maskNonOrthCorrection(p_rghEqn, ib.ibTouchMask());
            ib.fixSolidCells(p_rghEqn);

            p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA - p_rghEqn.flux();
                ib.maskFlux(phi);
            }
        }
    }

    #include "continuityErrs.H"

    p_rgh.relax();

    // Velocity from the open-face fluxes: the pressure flux is
    // phiHbyA - phi; the blocked faces (rAUf = 0) are masked before the
    // division and carry nothing in the reconstruction
    U = HbyA
      + rAU*ib.reconstruct
        (
            ib.faceMask()*(phig - (phiHbyA - phi))
           /(rAUf + dimensionedScalar(rAUf.dimensions(), VSMALL))
        );
    U.correctBoundaryConditions();
    fvOptions.correct(U);

    // Diagnostic (PIMPLE ibMomentumBudget): balance of the IB cell with
    // the largest velocity, to find what feeds a runaway pit cell
    if (pimple.dict().getOrDefault<bool>("ibMomentumBudget", false))
    {
        const labelList& ibc = ib.ibCells();
        const vectorField& Ui = U.primitiveField();
        label iMax = -1;
        scalar uMax = -1;
        forAll(ibc, i)
        {
            if (mag(Ui[ibc[i]]) > uMax)
            {
                uMax = mag(Ui[ibc[i]]);
                iMax = i;
            }
        }
        const scalar gUMax = returnReduce(uMax, maxOp<scalar>());
        // Field construction evaluates the processor patches: all ranks
        const volScalarField AU(UEqn.A());
        if (iMax >= 0 && uMax == gUMax)
        {
            const label c = ibc[iMax];
            const scalar A = AU[c];
            const scalar Vc = mesh.V()[c];
            const scalar sink =
                ib.wallModel().wallCoeff()[iMax]*ib.ibWallArea()[iMax]
               *rho[c]/Vc;
            scalar sumPhi = 0;
            scalar maxPhi = 0;
            label nOpen = 0;
            for (const label facei : mesh.cells()[c])
            {
                if (facei < mesh.nInternalFaces())
                {
                    const scalar f = phi[facei];
                    sumPhi += (mesh.owner()[facei] == c ? f : -f);
                    maxPhi = max(maxPhi, mag(f));
                    nOpen += (ib.faceMask()[facei] > 0.5);
                }
            }
            Pout<< "IB budget: cell " << c << " at " << mesh.C()[c]
                << " dist " << ib.distance()[c]
                << " |U| " << uMax << " |U_old| " << mag(U.oldTime()[c])
                << " |HbyA| " << mag(HbyA[c])
                << " |rAU recon| " << mag(Ui[c] - HbyA[c])
                << " A " << A << " rho/dt " << rho[c]/runTime.deltaTValue()
                << " sink " << sink << " rAU " << rAU[c]
                << " openFaces " << nOpen << " maxPhi " << maxPhi
                << " sumPhi " << sumPhi << " V " << Vc << endl;
        }
    }
    ib.correctU(U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}
