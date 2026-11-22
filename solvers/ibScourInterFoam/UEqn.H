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
    Momentum equation of interFoam with the immersed boundary: the wall
    sink of the IB cells is density weighted and the pressure and gravity
    fluxes of the predictor are reconstructed from the open faces.

\*---------------------------------------------------------------------------*/

MRF.correctBoundaryVelocity(U);

// Dynamic viscosity of the mixture for the cut-face corrections (the
// diffusion term of divDevRhoReff is laplacian(rho nuEff, U)); nut from
// the registry when a turbulence model holds one
const volScalarField muEff
(
    "muEff",
    rho
   *(
        mixture.nu()
      + (
            nutPtr
          ? *nutPtr
          : volScalarField::New("nutZero", mesh, dimensionedScalar(dimViscosity, Zero))()
        )
    )
);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(rho, U)
 ==
    fvOptions(rho, U)
    // Cut faces (height-field bed with apertures): explicit corrections of
    // the diffusion and convection fluxes on the faces of cut cells
    // (cutCentroids options; zero when off)
    // (the library labels them with the single-phase dimensions; the
    // values are the mass-weighted ones with muEff and rhoPhi)
  + dimensionedScalar(dimDensity, 1)*ib.cutSkewCorrection(U, muEff)
  + dimensionedScalar(dimDensity, 1)*ib.cutConvectionCorrection(rhoPhi, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

// Cut faces: diffusion through the open part of the faces only
ib.scaleCutDiffusion(UEqn, muEff);

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
      : static_cast<const tokenList&>
        (
            mesh.laplacianScheme("laplacian((rho*nuEff),U)")
        )
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
                ib.ibTouchMask()*fvc::interpolate(muEff)*mesh.magSf()
               *fv::correctedSnGrad<vector>(mesh).correction(U)
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
ib.constrainU(UEqn, rho);

if (pimple.momentumPredictor())
{
    fvVectorMatrix UEqnP
    (
        UEqn
     ==
        ib.reconstruct
        (
            (
                mixture.surfaceTensionForce()
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );
    ib.constrainU(UEqnP, rho, false);
    UEqnP.solve();

    fvOptions.correct(U);
}
