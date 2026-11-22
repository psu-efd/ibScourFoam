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

#include "roughLogLawWallLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(roughLogLawWallLaw, 0);
    addToRunTimeSelectionTable(ibWallLaw, roughLogLawWallLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::roughLogLawWallLaw::roughLogLawWallLaw(const dictionary& dict)
:
    ibWallLaw(dict),
    Cs_(dict.getOrDefault<scalar>("Cs", 0.5))
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::roughLogLawWallLaw::uPlus
(
    const scalar yPlus,
    const scalar ksPlus
) const
{
    const scalar Edash = E_/(1 + Cs_*ksPlus);

    // Intersection of u+ = y+ and the log law (fixed-point iteration)
    scalar yPlusLam = 11;
    for (label iter = 0; iter < 20; ++iter)
    {
        yPlusLam = Foam::log(max(Edash*yPlusLam, 1.0 + SMALL))/kappa_;
    }

    if (yPlus <= yPlusLam)
    {
        return yPlus;
    }

    // Bounded below: for y < ks/(E/Cs) the log law would be negative
    return max(Foam::log(Edash*yPlus)/kappa_, 0.1);
}

// ************************************************************************* //
