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
    Foam::ibSuspendedLoad

Description
    Suspended sediment over an immersed-boundary bed (formulation of
    libsuspendedSediment): the volume concentration C is carried by the
    drift flux phiC = phi + ws.Sf (settling only across faces between
    fluid cells) and diffuses with DC = nu/ScC + nut/sigmaC,

        ddt(C) + div(phiC, C) - laplacian(DC, C) = (E - D) A_w/V  in IB cells,

    with C = 0 in solid cells. The bed exchange of every IB cell uses the
    wall area A_w, wall distance y_c, normal n and shear velocity u* of
    the IB wall model, w_n = -ws.n the settling velocity through the bed,
    Gamma_1 = DC of the IB cell and the reference concentration C_b* of
    referenceConcentrationModel (sedimentProperties). The net flux into
    the fluid, with C_P the IB-cell value, is the two-point wall flux of
    the sedimentBedFlux condition with the face value eliminated:

        referenceLevel : E = w_n C_b*,                    net = E - w_n C_P
        HRN (firstCell): E = w_n C_b* (1 + P2) [or exp(P2)],
                         P2 = w_n delta/Gamma_1,
                         net = E (1 - Pe1) - w_n C_P,     Pe1 = w_n y_c/Gamma_1
                         (Liu 2014 Eqs. 18/17; stops when Pe1 >= 1 as sedimentBedFlux)
        logLaw         : net = w_n [C_b* (delta/y_c)^Z' - C_P],
                         Z' = sigma_c w_n/(kappa u*)  (extension of Liu 2014)
        ripple         : suspension over wave ripples (Nielsen 1992): the
                         ripple vortices mix the sand over a layer of a few
                         L_s = 1.4 eta (eta the ripple height) with the
                         diffusivity eps_s = w_s L_s, so the period-mean
                         concentration is C0 exp(-z/L_s) with
                         C0 = 0.005 theta_r^3 at the mean bed level,
                         theta_r = theta'_max/(1 - pi eta/lambda)^2 the
                         ripple-crest Shields number from the peak skin
                         friction of the previous wave period (wall-model
                         u*). E = w_n C0 and the deposition is that of the
                         bed value of the exponential profile through C_P,
                         net = w_n [C0 - C_P exp(y_c/L_s)]. The
                         referenceConcentrationModel is not used.

    The deposition w_n C_P is implicit. The net deposition D - E per IB cell
    (grain volume per wall area and time) goes to the bed model for the
    Exner equation. Only the total shear velocity of the wall model is used
    for the Shields number (shieldsNumber total).

Usage
    constant/sedimentProperties (entries of suspendedSedimentFoam) with
    \verbatim
    solveSediment   yes;
    ScC             1;
    sigmaC          1;
    d50             2e-4;
    s               2.65;
    settlingVelocityModel  vanRijn1984;     // constant (ws) | stokes | soulsby
    referenceConcentrationModel zysermanFredsoe;
    shieldsNumber   total;
    ibBedExchange
    {
        mode        HRN;            // referenceLevel | HRN | logLaw | ripple
        form        linear;         // HRN: linear | exponential
        kappa       0.41;           // logLaw
    }
    \endverbatim

    Ripple mode (ibBedExchange):
    \verbatim
        mode            ripple;
        period          3.5;        // wave period: window of theta'_max [s]
        mixingLayers    5;          // eps_s up to mixingLayers*L_s above bed
        ripples
        {
            model       fixed;      // fixed | nielsen1981
            height      0.01;       // fixed [m]
            length      0.05;       // fixed [m]
            // nielsen1981 (regular waves, mobility number psi < 156):
            // Um 0.174; (undisturbed near-bed velocity amplitude [m/s];
            // T is the period above)
        }
    \endverbatim

SourceFiles
    ibSuspendedLoad.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ibSuspendedLoad_H
#define Foam_ibSuspendedLoad_H

#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "settlingVelocityModel.H"
#include "referenceConcentrationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class immersedBoundary;

/*---------------------------------------------------------------------------*\
                        Class ibSuspendedLoad Declaration
\*---------------------------------------------------------------------------*/

class ibSuspendedLoad
{
    // Private Data

        const fvMesh& mesh_;
        const immersedBoundary& ib_;

        //- Schmidt numbers
        const scalar ScC_;
        const scalar sigmaC_;

        //- Bed-exchange mode and options
        const word mode_;
        const bool exponential_;
        const scalar kappa_;

        //- Gravity magnitude [m/s2]
        const scalar g_;

        //- Settling velocity [m/s] (vertical, downward)
        scalar wsMag_;
        vector ws_;

        //- Reference concentration model
        autoPtr<suspendedSediment::referenceConcentrationModel> refModel_;

        //- Concentration and diffusivity
        volScalarField C_;
        volScalarField DC_;

        //- Per IB cell: entrainment E, deposition D [m/s], C_b*
        scalarField E_;
        scalarField D_;
        scalarField Cstar_;

        //- Erodibility per IB cell (1 erodible, 0 fixed: no pickup);
        //- empty: all erodible
        scalarField erodible_;

        // Ripple mode (Nielsen 1992)

            //- Ripple height eta and length lambda [m]
            scalar rippleHeight_;
            scalar rippleLength_;

            //- Concentration decay length L_s = 1.4 eta [m]
            scalar Ls_;

            //- Thickness of the ripple mixing layer in units of L_s
            scalar mixingLayers_;

            //- Window of the peak Shields number (the wave period) [s]
            scalar ripplePeriod_;

            //- Grain diameter and specific gravity (Shields number)
            scalar d50_;
            scalar s_;

            //- Running maximum of the skin-friction Shields number in the
            //- current window, and that of the last completed window
            //- (per cell: the IB cell list changes as the bed moves)
            autoPtr<volScalarField> thetaWindowPtr_;
            autoPtr<volScalarField> thetaPeakPtr_;

            //- Time accumulated in the current window
            scalar rippleWindowTime_;

            //- A completed window exists (read back on a restart)
            bool haveThetaPeak_;


    // Private Member Functions

        //- Ripple mode: read the ripple geometry and allocate the fields
        void readRipples(const dictionary& dict);


public:

    // Private Member Functions

        //- The equation; alpha is null for the single-phase solver
        void solveImpl
        (
            const surfaceScalarField& phi,
            const scalar nu,
            const volScalarField* alpha
        );


public:

    // Constructors

        //- Construct from the mesh, sedimentProperties, the immersed
        //- boundary and the kinematic viscosity
        ibSuspendedLoad
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const immersedBoundary& ib,
            const scalar nu
        );


    // Member Functions

        //- Diffusivity from the turbulence fields
        void correctDiffusivity
        (
            const volScalarField& nu,
            const volScalarField& nut
        );

        //- Solve the concentration equation with the bed exchange
        void solve(const surfaceScalarField& phi, const scalar nu);

        //- Two-phase (VoF) variant: the sediment lives in the water only.
        //- The equation is for alpha C with the water flux phiWater (the
        //- flux that advected alpha), the settling and diffusive fluxes
        //- are weighted by the face water fraction, the bed exchange acts
        //- in wet IB cells only and C is zero in air
        void solve
        (
            const surfaceScalarField& phiWater,
            const scalar nu,
            const volScalarField& alpha
        );

        //- Erodibility per IB cell (order of ib.ibCells()): C_b* and E
        //- are zero where it is 0 (a fixed strip, a rigid apron)
        void setErodible(const scalarField& erodible)
        {
            erodible_ = erodible;
        }

        //- Net deposition D - E per IB cell [m/s] (for the bed model)
        tmp<scalarField> netDeposition() const
        {
            return D_ - E_;
        }

        const volScalarField& C() const noexcept { return C_; }

        scalar ws() const noexcept { return wsMag_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
