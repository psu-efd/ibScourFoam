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
    Diagnostics (immersedBoundaryDict writeMomentumTerms true): write the
    terms of the momentum equation per unit volume, (A U - b)/V, with the
    final U, p and phi and the discretisation of UEqn.H (cut-face
    diffusion, solid decoupling, wall sink, IB pressure gradient):

        momDdt + momConv + momDiff + momSink + momGradP = momResidual

    with momDiff = momDiffLap + momDiffExpl: the implicit Laplacian
    -laplacian(nuEff, U) (cut-face scaled, decoupled from the solids) and
    the explicit part -div(nuEff dev2(T(grad(U)))) of divDevReff (not
    masked).

\*---------------------------------------------------------------------------*/

{
    fvVectorMatrix ddtEqn(fvm::ddt(U));
    fvVectorMatrix convEqn(fvm::div(phi, U));
    tmp<fvVectorMatrix> tDiffEqn(turbulence->divDevReff(U));
    fvVectorMatrix& diffEqn = tDiffEqn.ref();

    ib.scaleCutDiffusion(diffEqn, turbulence->nuEff()());

    ib.decoupleSolid(ddtEqn);
    ib.decoupleSolid(convEqn);
    ib.decoupleSolid(diffEqn);

    // Implicit Laplacian part alone, with the same cut-face and solid
    // treatment
    const volScalarField nuEffDiag(turbulence->nuEff());
    fvVectorMatrix lapEqn(-fvm::laplacian(nuEffDiag, U));
    ib.scaleCutDiffusion(lapEqn, nuEffDiag);
    ib.decoupleSolid(lapEqn);

    volVectorField momDdt("momDdt", ddtEqn & U);
    volVectorField momConv("momConv", convEqn & U);
    volVectorField momDiff("momDiff", diffEqn & U);
    volVectorField momDiffLap("momDiffLap", lapEqn & U);
    volVectorField momDiffExpl("momDiffExpl", momDiff - momDiffLap);
    volVectorField momSink("momSink", ib.wallSinkTerm(U));
    volVectorField momGradP("momGradP", ib.grad(p));
    volVectorField momResidual
    (
        "momResidual",
        momDdt + momConv + momDiff + momSink + momGradP
    );

    momDdt.write();
    momConv.write();
    momDiff.write();
    momDiffLap.write();
    momDiffExpl.write();
    momSink.write();
    momGradP.write();
    momResidual.write();
}

// ************************************************************************* //
