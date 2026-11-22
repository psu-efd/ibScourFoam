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
    Foam::vanDriestWallLaw

Description
    Van Driest (1956) mixing-length velocity profile with the roughness
    coordinate shift of Rotta (1962) / Cebeci & Chang (1978), as used for
    the bed shear stress by Roulund et al. (2005, J. Fluid Mech. 534,
    eqs. 21-22):

        u+ = 2 int_0^y+ dy' / (1 + sqrt(1 + 4 kappa^2 (y' + dy+)^2
                                         [1 - exp(-(y' + dy+)/A)]^2))

        dy+ = 0.9 (sqrt(ks+) - ks+ exp(-ks+/6)),   5 < ks+ < 2000

    with A = 25, y+ measured from the top of the roughness elements,
    dy+ = 0 for ks+ <= 5 (hydraulically smooth) and dy+ held at its
    ks+ = 2000 value above. The integral is tabulated once in
    (ln(1 + y+), dy+) and interpolated bilinearly; beyond the table the
    log-law asymptote du+/dy+ = 1/(kappa (y+ + dy+)) is used.

    For the turbulence values of the IB cells the law also provides the
    Wilcox (1993) rough-wall coefficient of Roulund et al. (eqs. 19-20),
    omega_w = S_r u_tau^2/nu with

        S_r = (40/ks+)^3      ks+ < 20.2
        S_r = 100/ks+^0.85    ks+ > 20.2

    and the coordinate shift dy+ for the mixing length.

Usage
    \verbatim
    wallLaw
    {
        type            vanDriest;
        ks              0.0007;     // equivalent sand roughness [m]
        A               25;         // optional, damping constant
    }
    \endverbatim

SourceFiles
    vanDriestWallLaw.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vanDriestWallLaw_H
#define Foam_vanDriestWallLaw_H

#include "ibWallLaw.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class vanDriestWallLaw Declaration
\*---------------------------------------------------------------------------*/

class vanDriestWallLaw
:
    public ibWallLaw
{
    // Private Data

        //- Van Driest damping constant
        scalar A_;

        //- Table abscissa t = ln(1 + y+): number of points and spacing
        label nT_;
        scalar dT_;

        //- Table shift values dy+: number of points and spacing
        label nShift_;
        scalar dShift_;

        //- Tabulated u+(t_i, dy+_j), row-major in j
        scalarList table_;


    // Private Member Functions

        //- Integrand du+/dy+ at y+ for shift dy+
        scalar dUPlus(const scalar yPlus, const scalar shift) const;

        //- Fill the table by Simpson integration
        void tabulate();

        //- Tabulated u+ at y+ for shift dy+
        scalar lookup(const scalar yPlus, const scalar shift) const;


public:

    //- Runtime type information
    TypeName("vanDriest");


    // Constructors

        //- Construct from dictionary
        explicit vanDriestWallLaw(const dictionary& dict);


    //- Destructor
    virtual ~vanDriestWallLaw() = default;


    // Member Functions

        //- Dimensionless velocity u+ at y+ for roughness ks+
        virtual scalar uPlus(const scalar yPlus, const scalar ksPlus) const;

        //- Roughness coordinate shift dy+ for ks+
        virtual scalar yPlusShift(const scalar ksPlus) const;

        //- Wilcox rough-wall coefficient S_r for ks+
        virtual scalar omegaWallCoeff(const scalar ksPlus) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
