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
    Momentum predictor with the velocity of IB and solid cells fixed.

\*---------------------------------------------------------------------------*/

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
    // Cut faces: explicit non-orthogonal correction of the diffusion flux
    // (cutCentroids skewCorrection; zero when off)
  + ib.cutSkewCorrection(U, turbulence->nuEff()())
    // Cut faces: convection flux at the open-face centroid
    // (cutCentroids faceInterpolation; zero when off)
  + ib.cutConvectionCorrection(phi, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

// Cut faces: diffusion through the open part of the faces only
ib.scaleCutDiffusion(UEqn, turbulence->nuEff()());

// Faces of IB and solid cells: the explicit non-orthogonal correction of
// the momentum Laplacian (no flux correction is stored for U) is rebuilt
// and removed there; see immersedBoundary::maskNonOrthCorrection
if (ib.maskNonOrthCorrection())
{
    bool corrected = false;
    bool limited = false;
    for
    (
        const token& tok
      : static_cast<const tokenList&>(mesh.laplacianScheme("laplacian(nuEff,U)"))
    )
    {
        if (tok.isWord())
        {
            corrected = corrected || (tok.wordToken() == "corrected");
            limited = limited || (tok.wordToken() == "limited");
        }
    }
    if (corrected && !limited)
    {
        UEqn +=
            fvc::div
            (
                ib.ibTouchMask()*fvc::interpolate(turbulence->nuEff()())
               *mesh.magSf()*fv::correctedSnGrad<vector>(mesh).correction(U)
            );
    }
    else if (corrected)
    {
        WarningInFunction
            << "maskNonOrthCorrection: the limited correction of the"
            << " momentum Laplacian is not removed on the IB faces" << endl;
    }
}

UEqn.relax();

fvOptions.constrain(UEqn);

// IB cells: wall-model velocity; solid cells: solid velocity. Applied to
// UEqn so that H/A of the pressure equation also carries these values.
ib.constrainU(UEqn);

if (pimple.momentumPredictor())
{
    // Pressure gradient with zero normal gradient at blocked faces
    fvVectorMatrix UEqnP(UEqn == -ib.grad(p));

    // Fixed values again (the source changed); the wall sink is already
    // in UEqn
    ib.constrainU(UEqnP, false);

    UEqnP.solve();

    fvOptions.correct(U);
}

// ************************************************************************* //
