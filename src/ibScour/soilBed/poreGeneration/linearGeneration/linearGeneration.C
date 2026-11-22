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

#include "linearGeneration.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace poreGenerationModels
{
    defineTypeNameAndDebug(linearGeneration, 0);
    addToRunTimeSelectionTable(poreGenerationModel, linearGeneration, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poreGenerationModels::linearGeneration::linearGeneration
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    poreGenerationModel(mesh, dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::poreGenerationModels::linearGeneration::generated
(
    const scalar sigmaV0,
    const scalar Nl,
    const scalar N
) const
{
    return sigmaV0*min(N/Nl, scalar(1));
}


Foam::scalar Foam::poreGenerationModels::linearGeneration::cyclesFor
(
    const scalar sigmaV0,
    const scalar Nl,
    const scalar u
) const
{
    return Nl*min(max(u/sigmaV0, scalar(0)), scalar(1));
}


// ************************************************************************* //
