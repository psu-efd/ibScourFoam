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

Class
    Foam::seepageModels::chengChiew

Description
    Seepage reduces the effective weight of the grains (Cheng & Chiew
    1999; Lu & Chiew 2007): the flat-bed critical Shields number becomes
    theta_c0 (1 - i/i_c) for an upward gradient i (downward seepage,
    i < 0, stabilises), clipped at zero, and, if 'reposeEffect' is true
    (default), the friction coefficient tan(phi_r) of the slope
    correction and of the sand slide is reduced by the same factor.

    References:
        Cheng, N.-S., & Chiew, Y.-M. (1999). Incipient sediment motion
        with upward seepage. J. Hydraul. Res., 37(5), 665-681.
        Lu, Y., & Chiew, Y.-M. (2007). Seepage effects on dune
        dimensions. J. Hydraul. Eng., 133(5), 560-563.

SourceFiles
    chengChiew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_chengChiew_H
#define Foam_chengChiew_H

#include "seepageModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace seepageModels
{

/*---------------------------------------------------------------------------*\
                          Class chengChiew Declaration
\*---------------------------------------------------------------------------*/

class chengChiew
:
    public seepageModel
{
protected:

    // Protected Data

        //- Apply the factor to the friction coefficient as well
        const bool reposeEffect_;

        //- Lower bound of the repose factor (default 0.05)
        const scalar minReposeFactor_;


    // Protected Member Functions

        //- Effective-weight factor 1 - i/i_c, clipped at zero
        scalar weightFactor(const scalar i) const
        {
            return max(1 - i/iCrit_, scalar(0));
        }



public:

    //- Runtime type information
    TypeName("chengChiew");


    // Constructors

        //- Construct from dictionary and the critical gradient
        chengChiew(const dictionary& dict, const scalar iCrit);


    //- Destructor
    virtual ~chengChiew() = default;


    // Member Functions

        //- 1 - i/i_c, clipped at zero
        virtual scalar thetaCFactor(const scalar i) const
        {
            return weightFactor(i);
        }

        //- The same factor if reposeEffect (floored at minReposeFactor so
        //- that the slope correction stays defined), else 1
        virtual scalar reposeFactor(const scalar i) const
        {
            return
            (
                reposeEffect_
              ? max(weightFactor(i), minReposeFactor_)
              : scalar(1)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace seepageModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
