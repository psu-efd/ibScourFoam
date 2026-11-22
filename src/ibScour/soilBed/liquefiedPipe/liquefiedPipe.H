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
    Foam::liquefiedPipe

Description
    Reduced-order vertical motion of a pipeline (or other long horizontal
    body) in a wave-liquefied seabed, driven by the soil region:
    sinking of heavy pipes, floatation of light ones (Sumer et
    al. 1999, 2006; Teh et al. 2006).

    Per unit length, a rigid pipe of diameter D and specific gravity s_p at
    one soil column moves vertically with the quasi-static velocity at
    which the drag balances its net weight
        f = rho_w g (A_w + s_l A_s) - s_p rho_w g A,
    A_w and A_s its areas above and below the bed (water and liquefied soil
    of specific gravity s_l). Drag laws:
      - sumer1999: the empirical fit of Sumer et al. (1999, Fig. 15),
        C_D = a/sqrt(Re_p), Re_p = w D/nu_w, F = 0.5 s_l rho_w C_D D w^2
        (a = 8e7);
      - linear: F = c w.
    The pipe moves only where the soil is mobilised: the period mean of the
    excess pore pressure of the column at the pipe centre (the top cell for
    a pipe on the bed) has reached onsetRatio sigma'_0, with
    sigma'_0 = sigma'_v0 (1 + 2 k0)/3 (Sumer et al. 1999: motion starts at
    0.4-0.7 sigma'_0). A stable (re-solidified) layer grows from the soil
    bottom at stableLayerSpeed once the bottom cell has reached sigma'_0
    (Teh et al. 2006); the pipe comes to rest on it. A floating pipe stops
    with its bottom at the bed.

    Benchmark: Sumer et al. (1999), with
    stableLayerSpeed 1e-4 m/s calibrated on their Fig. 9, predicts the
    burial depths of their Figs. 12-13 within 0.3-0.6 D; the depth limit of
    very heavy pipes (s_p > 4-5) is not reproduced.

    The pipe is not part of the soil mesh (no structure mask): a reduced
    model on top of the soil region's liquefaction. Driving an immersed
    solid (heightFieldBed 'liquefiedPipes'), the burial is set from the
    solid and the current bed before each update (setBurial), so that scour
    or deposition at the pipe is accounted for; a pipe heavier than water
    whose bottom is above the bed (scoured from under it) is unsupported and
    drops onto the bed at most at freeFallSpeed.

Usage
    \verbatim
    liquefiedPipes          // one sub-dictionary per pipe
    {
      pipe1
      {
        x               0.45;       // position of the column [m]
        diameter        0.04;       // D [m]
        specificGravity 3.1;        // s_p
        initialBurial   0.07;       // pipe-bottom depth [m] (0: on the bed)
        liquefiedSpecificGravity 1.3;   // s_l
        period          1.6;        // wave period for the means [s]
        onsetRatio      0.5;
        stableLayerSpeed 1e-4;      // [m/s] (0: none)
        drag            sumer1999;  // | linear
        dragCoefficient 8e7;        // a (sumer1999) or c [N s/m2] (linear)
        freeFallSpeed   0.1;        // [m/s] (default 0.1)
      }
    }
    \endverbatim

SourceFiles
    liquefiedPipe.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_liquefiedPipe_H
#define Foam_liquefiedPipe_H

#include "soilRegion.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class liquefiedPipe Declaration
\*---------------------------------------------------------------------------*/

class liquefiedPipe
{
    // Private Data

        //- The soil region
        const soilRegion& soil_;

        //- Diameter [m], specific gravities of the pipe and the liquefied
        //- soil, wave period [s], onset ratio, stable-layer speed [m/s]
        const scalar D_;
        const scalar sp_;
        const scalar sl_;
        const scalar period_;
        const scalar onsetRatio_;
        const scalar stableSpeed_;

        //- Drag law and its coefficient
        const word drag_;
        const scalar dragCoeff_;

        //- Largest speed of an unsupported pipe dropping onto the bed [m/s]
        const scalar freeFall_;

        //- Water density [kg/m3], gravity [m/s2], water viscosity [m2/s]
        const scalar rhoW_;
        const scalar g_;
        const scalar nuW_;

        //- Column of the pipe on this processor (-1: elsewhere)
        label col_;

        //- Depth of the pipe bottom below the bed [m] and velocity
        //- (positive down) [m/s]
        scalar e_;
        scalar w_;

        //- Running time integral of the column excess pore pressure over
        //- the current period, and the mean of the last complete period
        scalarField pSum_;
        scalarField pMean_;
        scalar tPeriod_;
        bool havePeriod_;

        //- Start of the stable layer [s] (GREAT: not yet)
        scalar tStable_;

        //- Vertical external force per unit length [N/m] (up positive)
        scalar fExt_;

        //- Last diagnostics: ratio at the centre, stable-layer top depth
        scalar ratio_;
        scalar stableTop_;

        //- Trajectory file (master)
        autoPtr<OFstream> file_;


    // Private Member Functions

        //- Net upward force per unit length [N/m] for a bottom depth e
        scalar netForce(const scalar e) const;

        //- Speed [m/s] at which the drag balances the force magnitude f
        scalar speed(const scalar f) const;


public:

    // Constructors

        //- Construct from the soil region, a name (of the trajectory file
        //- postProcessing/liquefiedPipe_<name>.dat) and the dictionary;
        //- a positive diameter or initial burial replaces the dictionary
        //- entry (taken from the geometry of an immersed solid)
        liquefiedPipe
        (
            const soilRegion& soil,
            const word& name,
            const dictionary& dict,
            const scalar diameter = -1,
            const scalar initialBurial = -1
        );


    // Member Functions

        //- Set the depth of the pipe bottom below the current bed [m]
        //- (negative: above the bed) from the geometry of the structure
        void setBurial(const scalar e) noexcept { e_ = e; }

        //- Set the vertical external force per unit length [N/m] (up
        //- positive), e.g. the hydrodynamic lift of the flow, added to the
        //- net weight
        void setExternalForce(const scalar f) noexcept { fExt_ = f; }

        //- Advance over the time step dt [s] after the soil has evolved
        void update(const scalar dt);

        //- Depth of the pipe bottom below the bed [m]
        scalar bottomDepth() const noexcept { return e_; }

        //- Velocity [m/s] (positive down)
        scalar velocity() const noexcept { return w_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
