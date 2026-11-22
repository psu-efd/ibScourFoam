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

Class
    Foam::poreGenerationModels::seedRahman

Description
    Seed et al. (1976) / Seed & Rahman (1978) generation law,

        u_g/sigma'_v0 = (2/pi) arcsin[(N/N_l)^(1/(2 theta))],

    with theta = 0.7 by default (u_g = sigma'_v0 at N = N_l).

SourceFiles
    seedRahman.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_poreGenerationModels_seedRahman_H
#define Foam_poreGenerationModels_seedRahman_H

#include "poreGenerationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace poreGenerationModels
{

/*---------------------------------------------------------------------------*\
                         Class seedRahman Declaration
\*---------------------------------------------------------------------------*/

class seedRahman
:
    public poreGenerationModel
{
    // Private Data

        //- Shape parameter theta
        const scalar theta_;


public:

    //- Runtime type information
    TypeName("seedRahman");


    // Constructors

        //- Construct from the soil mesh and dictionary
        seedRahman(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~seedRahman() = default;


    // Member Functions

        virtual scalar generated
        (
            const scalar sigmaV0,
            const scalar Nl,
            const scalar N
        ) const;

        virtual scalar cyclesFor
        (
            const scalar sigmaV0,
            const scalar Nl,
            const scalar u
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace poreGenerationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
