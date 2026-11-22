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
    Foam::seepageModels::pressureGradient

Description
    Seepage force as a body force on the grains: the pore-pressure
    gradient just below the bed, grad p = (grad_h p, dp/dz), acts on the
    surface layer in addition to the submerged weight gamma'. Its vertical
    (upward) part reduces the effective weight as in chengChiew, theta_c0
    (1 - i/i_c), and its horizontal part tilts the effective gravity: the
    bed-load slope correction and the sand slide see an apparent slope
    grad_h p / (gamma' (1 - i/i_c)) added to grad z_b (Sleath 1999; Yeh &
    Mason 2014; the non-hydrostatic pressure effect of Francalanci et al.
    2008 in its simplest form). The horizontal gradient is that of the
    pore pressure at the gradient-evaluation depth of the soil region
    (which includes the flow pressure on the bed).

    References:
        Sleath, J. F. A. (1999). Conditions for plug formation in
        oscillatory flow. Cont. Shelf Res., 19(13), 1643-1664.
        Yeh, H., & Mason, H. B. (2014). Sediment response to tsunami
        loading: mechanisms and estimates. Geotechnique, 64(2), 131-143.
        Francalanci, S., Parker, G., & Solari, L. (2008). Effect of
        seepage-induced nonhydrostatic pressure distribution on bed-load
        transport and bed morphodynamics. J. Hydraul. Eng., 134(4),
        378-389.

Usage
    \verbatim
    seepage
    {
        type            pressureGradient;
        reposeEffect    true;
        minReposeFactor 0.05;
        horizontalForce true;           // apparent slope from grad_h p
    }
    \endverbatim

SourceFiles
    pressureGradient.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_seepagePressureGradient_H
#define Foam_seepagePressureGradient_H

#include "chengChiew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace seepageModels
{

/*---------------------------------------------------------------------------*\
                       Class pressureGradient Declaration
\*---------------------------------------------------------------------------*/

class pressureGradient
:
    public chengChiew
{
    // Private Data

        //- Apply the horizontal seepage force as an apparent slope
        const bool horizontalForce_;


public:

    //- Runtime type information
    TypeName("pressureGradient");


    // Constructors

        //- Construct from dictionary and the critical gradient
        pressureGradient(const dictionary& dict, const scalar iCrit);


    //- Destructor
    virtual ~pressureGradient() = default;


    // Member Functions

        //- True: the sediment core adds grad_h p / (gamma' f) to the bed
        //- slope, f the effective-weight factor
        virtual bool horizontalForce() const
        {
            return horizontalForce_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace seepageModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
