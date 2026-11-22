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

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ibScourFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "seedRahman.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace poreGenerationModels
{
    defineTypeNameAndDebug(seedRahman, 0);
    addToRunTimeSelectionTable(poreGenerationModel, seedRahman, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poreGenerationModels::seedRahman::seedRahman
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    poreGenerationModel(mesh, dict),
    theta_(dict.getOrDefault<scalar>("theta", 0.7))
{
    Info<< "    Seed-Rahman generation: theta " << theta_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::poreGenerationModels::seedRahman::generated
(
    const scalar sigmaV0,
    const scalar Nl,
    const scalar N
) const
{
    if (N >= Nl)
    {
        return sigmaV0;
    }

    const scalar x = Foam::pow(N/Nl, 1/(2*theta_));

    return sigmaV0*2/constant::mathematical::pi*Foam::asin(min(x, scalar(1)));
}


Foam::scalar Foam::poreGenerationModels::seedRahman::cyclesFor
(
    const scalar sigmaV0,
    const scalar Nl,
    const scalar u
) const
{
    // u/sigma'_v0 = (2/pi) asin(x), x = (N/Nl)^(1/(2 theta))
    const scalar x = Foam::sin(0.5*constant::mathematical::pi*min(max(u/sigmaV0, scalar(0)), scalar(1)));
    return Nl*Foam::pow(x, 2*theta_);
}


// ************************************************************************* //
