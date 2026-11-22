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
    Foam::seepageModels::francalanci

Description
    Seepage model with the generalised Shields number of Francalanci,
    Parker & Solari (2008): the seepage-induced non-hydrostatic pressure
    gradient reduces the effective weight of the grains, so the Shields
    number of the transport closure is

        theta_g = theta / f(i),   f(i) = max(1 - i/i_c, f_min),

    with the flat-bed critical Shields number unchanged (the incipient
    condition theta_g >= theta_c0 is the Cheng-Chiew one, theta >=
    theta_c0 f, but the transport rate above threshold grows with the
    reduced weight). The friction coefficient and the horizontal seepage
    force (apparent slope) are those of the pressureGradient model.

    Reference:
        Francalanci, S., Parker, G., & Solari, L. (2008). Effect of
        seepage-induced nonhydrostatic pressure distribution on bed-load
        transport and bed morphodynamics. J. Hydraul. Eng., 134(4),
        378-389.

Usage
    \verbatim
    seepage
    {
        type            francalanci;
        minReposeFactor 0.05;       // optional, floor of f
        horizontalForce true;       // optional
    }
    \endverbatim

SourceFiles
    francalanci.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_seepageModels_francalanci_H
#define Foam_seepageModels_francalanci_H

#include "pressureGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace seepageModels
{

/*---------------------------------------------------------------------------*\
                         Class francalanci Declaration
\*---------------------------------------------------------------------------*/

class francalanci
:
    public pressureGradient
{
public:

    //- Runtime type information
    TypeName("francalanci");


    // Constructors

        //- Construct from dictionary and the critical gradient
        francalanci(const dictionary& dict, const scalar iCrit);


    //- Destructor
    virtual ~francalanci() = default;


    // Member Functions

        //- The flat-bed critical Shields number is unchanged
        virtual scalar thetaCFactor(const scalar i) const
        {
            return 1;
        }

        //- Generalised Shields number: theta / f(i), f floored at
        //- minReposeFactor
        virtual scalar thetaFactor(const scalar i) const
        {
            return 1/max(weightFactor(i), minReposeFactor_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace seepageModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
