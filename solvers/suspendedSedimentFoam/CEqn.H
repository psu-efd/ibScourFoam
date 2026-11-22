/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Xiaofeng Liu, Penn State University
-------------------------------------------------------------------------------
    Part of suspendedSedimentFoam (no code from OpenFOAM in this file).
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

Description
    Suspended-sediment transport: the volume concentration C is carried by
    the drift flux

        phiC = phi + ws & Sf

    (flow plus settling)
    and diffuses with DC = nu/ScC + nut/sigmaC:

        ddt(C) + div(phiC, C) - laplacian(DC, C) = fvOptions(C)

    With the bed condition option "diffusivity firstCell" (or logLaw) the
    bed-face value of DC is replaced by the value the condition provides.

    The settling part of phiC is kept only on patches where C has the
    condition sedimentBedFlux: there the sediment deposits (D = beta C_f)
    and is entrained (E). On all other non-coupled patches (walls, surface,
    inlets, outlets) phiC = phi, so sediment does not settle through them.
    With the condition in mode nearestCell the face carries no flux either
    and the exchange is a source of the wall cell, E |Sf|/V explicit and
    the deposition w_n C_P |Sf|/V implicit (Liu and Garcia 2008 Eq. 15; the
    form of ibSuspendedLoad).

    The equation is conservative: the change of int C dV equals the net
    boundary flux, which the function object sedimentBedExchange records.

\*---------------------------------------------------------------------------*/

if (solveSediment)
{
    if (!fixedDiffusivity)
    {
        DC = turbulence->nu()/ScC + turbulence->nut()/sigmaC;
    }

    // Bed faces with "diffusivity firstCell" (eps_s of the wall cell across
    // the near-wall layer, Liu 2014) or "logLaw" instead of the patch value
    forAll(C.boundaryField(), patchi)
    {
        const auto* bedPtr = dynamic_cast<const sedimentBedFluxFvPatchScalarField*>
        (
            &C.boundaryField()[patchi]
        );

        if (bedPtr)
        {
            tmp<scalarField> tGamma(bedPtr->faceDiffusivity());

            if (tGamma)
            {
                DC.boundaryFieldRef()[patchi] = tGamma;
            }
        }
    }

    surfaceScalarField phiC("phiC", phi + (ws & mesh.Sf()));

    surfaceScalarField::Boundary& phiCbf = phiC.boundaryFieldRef();

    // Bed exchange as a source of the wall cells (mode nearestCell)
    volScalarField::Internal CSu
    (
        IOobject("CSu", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(dimless/dimTime, Zero)
    );
    volScalarField::Internal CSp
    (
        IOobject("CSp", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(dimless/dimTime, Zero)
    );

    forAll(phiCbf, patchi)
    {
        if (phiCbf[patchi].coupled())
        {
            continue;
        }

        auto* bedPtr = dynamic_cast<sedimentBedFluxFvPatchScalarField*>
        (
            &C.boundaryFieldRef()[patchi]
        );

        if (!bedPtr)
        {
            phiCbf[patchi] = phi.boundaryField()[patchi];
        }
        else if (bedPtr->cellExchange())
        {
            phiCbf[patchi] = phi.boundaryField()[patchi];

            // E and w_n of this step (the matrix constructor would update
            // the condition only after the sources are assembled)
            bedPtr->updateCoeffs();

            const fvPatch& p = bedPtr->patch();
            const labelUList& cells = p.faceCells();
            const scalarField& magSf = p.magSf();
            const scalarField& E = bedPtr->E();
            const scalarField& wn = bedPtr->betaF();

            forAll(cells, facei)
            {
                const scalar AoV = magSf[facei]/mesh.V()[cells[facei]];
                CSu[cells[facei]] += E[facei]*AoV;
                CSp[cells[facei]] -= wn[facei]*AoV;
            }
        }
    }

    fvScalarMatrix CEqn
    (
        fvm::ddt(C)
      + fvm::div(phiC, C)
      - fvm::laplacian(DC, C)
     ==
        fvOptions(C)
      + CSu
      + fvm::Sp(CSp, C)
    );

    CEqn.relax();

    fvOptions.constrain(CEqn);

    CEqn.solve();

    fvOptions.correct(C);

    Info<< "C: min = " << gMin(C.primitiveField())
        << ", max = " << gMax(C.primitiveField())
        << ", volume integral = " << fvc::domainIntegrate(C).value() << endl;
}

// ************************************************************************* //
