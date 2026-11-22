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

\*---------------------------------------------------------------------------*/

#include "wallLawFaces.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::wallLawFaces::evaluate
(
    const ibWallLaw& law,
    const fvPatchVectorField& Up,
    const scalarField& nuw,
    scalarField& uTau,
    scalarField& y,
    scalarField& magUt,
    label* nClamped,
    scalar* minRatio
)
{
    const fvPatch& p = Up.patch();
    const vectorField n(p.nf());
    const vectorField d(p.Cf() - p.Cn());
    const vectorField dU(Up.patchInternalField() - Up);

    uTau.resize(p.size());
    y.resize(p.size());
    magUt.resize(p.size());

    forAll(n, facei)
    {
        // Normal distance, bounded below on distorted cells (cell centre
        // close to the face plane), where it would give omega_vis -> inf
        const scalar yn = n[facei] & d[facei];
        const scalar magD = mag(d[facei]);
        y[facei] = max(yn, yMinFraction*magD);

        if (nClamped && yn < yMinFraction*magD)
        {
            ++(*nClamped);
        }
        if (minRatio && magD > VSMALL)
        {
            *minRatio = min(*minRatio, yn/magD);
        }
        magUt[facei] = mag(dU[facei] - n[facei]*(n[facei] & dU[facei]));
        uTau[facei] = law.uTau(magUt[facei], y[facei], nuw[facei]);
    }
}


// ************************************************************************* //
