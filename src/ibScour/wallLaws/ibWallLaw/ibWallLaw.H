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
    Foam::ibWallLaw

Description
    Abstract wall law u+ = f(y+, ks+) used by the immersed-boundary wall
    model, with a robust solver for the friction velocity.

    Given the tangential velocity magnitude Ut at wall distance y, the
    friction velocity u_tau solves

        Ut = u_tau u+(y u_tau/nu, ks u_tau/nu)

    The left-hand side minus the right-hand side decreases monotonically in
    u_tau for the laws implemented, so the root is bracketed and found by a
    safeguarded secant (regula falsi) iteration.

Usage
    \verbatim
    wallLaw
    {
        type            spalding;   // spalding, roughLogLaw, vanDriest, laminar
        kappa           0.41;       // optional
        E               9.8;        // optional
    }
    \endverbatim

SourceFiles
    ibWallLaw.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ibWallLaw_H
#define Foam_ibWallLaw_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class ibWallLaw Declaration
\*---------------------------------------------------------------------------*/

class ibWallLaw
{
protected:

    // Protected Data

        //- von Karman constant
        scalar kappa_;

        //- Log-law constant E (B = ln(E)/kappa)
        scalar E_;

        //- Equivalent sand roughness ks [m]
        scalar ks_;


public:

    //- Runtime type information
    TypeName("ibWallLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            ibWallLaw,
            dictionary,
            (const dictionary& dict),
            (dict)
        );


    // Constructors

        //- Construct from dictionary
        explicit ibWallLaw(const dictionary& dict);


    // Selectors

        //- Return the selected wall law
        static autoPtr<ibWallLaw> New(const dictionary& dict);


    //- Destructor
    virtual ~ibWallLaw() = default;


    // Member Functions

        //- von Karman constant
        scalar kappa() const noexcept
        {
            return kappa_;
        }

        //- Log-law constant E
        scalar E() const noexcept
        {
            return E_;
        }

        //- Equivalent sand roughness [m]
        scalar ks() const noexcept
        {
            return ks_;
        }

        //- Dimensionless velocity u+ at y+ for roughness ks+
        virtual scalar uPlus(const scalar yPlus, const scalar ksPlus) const
            = 0;

        //- Roughness coordinate shift dy+ of the mixing length (zero for
        //- laws without a shifted origin)
        virtual scalar yPlusShift(const scalar ksPlus) const
        {
            return 0;
        }

        //- Coefficient S_r of the wall value omega_w = S_r u_tau^2/nu
        //- (Wilcox rough-wall condition); VGREAT for a smooth wall
        virtual scalar omegaWallCoeff(const scalar ksPlus) const
        {
            return VGREAT;
        }

        //- Friction velocity [m/s] from the tangential velocity Ut [m/s]
        //- at wall distance y [m], kinematic viscosity nu [m2/s]
        virtual scalar uTau
        (
            const scalar Ut,
            const scalar y,
            const scalar nu
        ) const;

        //- Turbulence of a wall cell at distance y [m] for friction
        //- velocity uTau [m/s] (shared by the immersed-boundary wall model
        //- and the body-fitted wall conditions):
        //- omega = sqrt(omega_vis^2 + omega_log^2), omega_vis =
        //- 6 nu/(beta1 (y + ys)^2), omega_log = u_tau/(sqrt(Cmu) kappa
        //- (y + dy)), with the roughness shift dy of the mixing length and
        //- the offset ys that gives omega_w = S_r u_tau^2/nu at the wall;
        //- nutLog = kappa u_tau (y + dy)(1 - exp(-(y + dy)+/A+))^2
        void wallTurbulence
        (
            const scalar uTau,
            const scalar y,
            const scalar nu,
            const scalar Cmu,
            const scalar beta1,
            const scalar Aplus,
            scalar& omega,
            scalar& nutLog
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
