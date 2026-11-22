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

#include "spaldingWallLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(spaldingWallLaw, 0);
    addToRunTimeSelectionTable(ibWallLaw, spaldingWallLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spaldingWallLaw::spaldingWallLaw(const dictionary& dict)
:
    ibWallLaw(dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::spaldingWallLaw::uPlus
(
    const scalar yPlus,
    const scalar ksPlus
) const
{
    const scalar k = kappa_;
    const scalar expMinusKB = 1.0/E_;   // exp(-kappa B) = 1/E

    // y+ as a function of u+ and its derivative
    auto f = [&](const scalar up)
    {
        const scalar ku = k*up;
        return
            up
          + expMinusKB
           *(Foam::exp(ku) - 1 - ku - 0.5*sqr(ku) - pow3(ku)/6.0);
    };
    auto dfdu = [&](const scalar up)
    {
        const scalar ku = k*up;
        return 1 + expMinusKB*k*(Foam::exp(ku) - 1 - ku - 0.5*sqr(ku));
    };

    if (yPlus <= 0)
    {
        return 0;
    }

    // Initial guess: the smaller of the linear and log laws
    scalar up = min(yPlus, Foam::log(max(E_*yPlus, 1.0 + SMALL))/k);

    for (label iter = 0; iter < 50; ++iter)
    {
        const scalar dUp = (f(up) - yPlus)/dfdu(up);
        up = max(up - dUp, 0.5*up);

        if (mag(dUp) < 1e-10*max(up, 1.0))
        {
            break;
        }
    }

    return up;
}

// ************************************************************************* //
