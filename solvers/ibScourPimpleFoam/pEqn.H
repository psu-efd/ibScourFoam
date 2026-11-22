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
    Pressure corrector with blocked immersed-boundary faces: no flux crosses
    a face between fluid and solid cells, and solid cells are pinned.

\*---------------------------------------------------------------------------*/

volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

// Cut faces: flux of HbyA at the open-face centroid of the faces of cut
// cells (cutCentroids faceInterpolation; no-op when off)
ib.correctCutFlux(HbyA, phiHbyA);

if (pimple.ddtCorr())
{
    if (ib.faceInterpolation())
    {
        // Cut faces: the Rhie-Chow transient correction compares the old
        // flux with the flux of the old velocity; on cut faces both must
        // be the open-face (masked, open-face centroid) fluxes, otherwise
        // the closed part of the face enters the correction every step
        surfaceScalarField phiU0(fvc::flux(U.oldTime()));
        ib.correctCutFlux(U.oldTime(), phiU0);
        ib.maskFlux(phiU0);

        const surfaceScalarField phiCorr(phi.oldTime() - phiU0);

        tmp<fv::ddtScheme<vector>> ddtU
        (
            fv::ddtScheme<vector>::New
            (
                mesh,
                mesh.ddtScheme("ddt(" + U.name() + ')')
            )
        );

        phiHbyA +=
            MRF.zeroFilter
            (
                fvc::interpolate(rAU)
               *ddtU.ref().fvcDdtPhiCoeff(U.oldTime(), phi.oldTime(), phiCorr)
               *phiCorr/mesh.time().deltaT()
            );
    }
    else
    {
        phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi));
    }
}

MRF.makeRelative(phiHbyA);

// No flux through faces of solid cells
ib.maskFlux(phiHbyA);

if (p.needReference())
{
    adjustPhi(phiHbyA, U, p);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAU, MRF);

// Laplacian coefficient with the blocked faces removed
const surfaceScalarField rAUf(ib.maskedLaplacianCoeff(rAU));

if (ib.pressureSubset())
{
    // Pressure solved on the fluid-only mesh subset (no solid cells in the
    // linear system); sets p and the conservative flux phi
    ib.solvePressure(p, rAUf, phiHbyA, phi, pimple, pRefCell, pRefValue);
}
else
{
    // Non-orthogonal pressure corrector loop on the full mesh with
    // penalised solid cells
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        // Solid cells are decoupled from the fluid
        ib.maskNonOrthCorrection(pEqn, ib.ibTouchMask());
        ib.fixSolidCells(pEqn);

        pEqn.solve(p.select(pimple.finalInnerIter()));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();

            // Faces between solid cells carry a correction flux of the
            // decoupled solid block; it is not part of the fluid
            ib.maskFlux(phi);
        }
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

ibGradP = ib.grad(p);
U = HbyA - rAU*ibGradP;
U.correctBoundaryConditions();
fvOptions.correct(U);

// IB and solid cells keep their constrained values; the flux is not
// modified (it remains the conservative, projected flux)
ib.correctU(U);

// ************************************************************************* //
