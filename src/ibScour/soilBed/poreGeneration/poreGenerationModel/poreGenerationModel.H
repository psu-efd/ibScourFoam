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
    Foam::poreGenerationModel

Description
    Residual pore-pressure generation (build-up) of the soil region under
    cyclic shearing (option 1): the undrained generation rate
    f = du_g/dt [Pa/s] of the pore pressure u_g that cyclic shearing of
    amplitude tau_cyc would produce without drainage, as a source of the
    pressure equation of the soil model (Seed & Rahman 1978; McDougal et
    al. 1989; Jeng, Seymour & Li 2007; Sumer 2014),

        du/dt = c_v laplacian(u) + f.

    The number of cycles to liquefaction N_l of a cell follows from its
    cyclic stress ratio, tau_cyc/sigma'_v0 = alpha N_l^(-beta) (no
    liquefaction below 'thresholdRatio': the asymptote of the cyclic
    strength curves, Alba et al. 1976, behind the critical wave height of
    Sumer et al. 1999, Fig. 5), and the
    generated pressure u_g(N) of the law (seedRahman, linear) from the
    number of cycles N. With equivalentCycles (default), N is the
    equivalent number of cycles that corresponds to the current excess
    pore pressure of the cell through the inverse of the law (Booker,
    Rahman & Seed 1976), bounded by the elapsed cycles: drainage lowers N
    and the generation goes on until the cell actually liquefies, while
    pressure fed by diffusion from a liquefied neighbour does not count
    as loading (unbounded, the liquefaction front runs away with the
    grid). When the shear amplitude is computed from the solved stress,
    the cyclic loading is also present in the solved pressure, and the
    equivalent cycles are counted from the period mean of the excess pore
    pressure (the residual pressure; updated when a period completes)
    rather than from its instantaneous value, which the momentary
    component would corrupt within every cycle. Without it, N = (t - t0)/T
    elapsed since the loading started (the form of the analytical solutions of
    McDougal et al. 1989 and Jeng et al. 2007) and the generation stops
    at N_l whatever the drainage. The rate of a step is the increment of
    u_g over the step (no singularity at N = 0 or N = N_l), capped so
    that the pressure does not exceed sigma'_v0; the cumulative generated
    pressure is kept as the field poreGenerated.

    Cyclic shear amplitude tau_cyc(z) per cell ('shearAmplitude'):
    - wave: tau = p0 k z exp(-k z), the progressive-wave shear stress in
      a saturated poro-elastic bed (amplitude p0, wavelength 2 pi/k; z
      the depth below the bed);
    - computed: the amplitude of the solved effective shear stress
      sigma'_xz over the last completed period (models that solve the
      displacement), from a running min/max envelope;
    - constant: a uniform value.

Usage
    \verbatim
    soil
    {
        ...
        poreGeneration
        {
            type            seedRahman;     // | linear
            period          8;              // T [s]
            startTime       0;              // t0 [s]
            equivalentCycles yes;           // N from the current pressure
            alpha           0.246;          // tau/sigma'_v0 = alpha N_l^-beta
            beta            0.165;
            thresholdRatio  0;              // no generation below (default 0)
            theta           0.7;            // seedRahman shape
            shearAmplitude
            {
                type        wave;           // | computed | constant
                amplitude   5000;           // p0 [Pa]
                wavelength  20;             // [m]
            }
        }
    }
    \endverbatim

SourceFiles
    poreGenerationModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_poreGenerationModel_H
#define Foam_poreGenerationModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class poreGenerationModel Declaration
\*---------------------------------------------------------------------------*/

class poreGenerationModel
{
protected:

    // Protected Data

        //- The soil mesh
        const fvMesh& mesh_;

        //- Settings
        dictionary dict_;

        //- Loading period T [s] and start time t0 [s]
        const scalar T_;
        const scalar t0_;

        //- Cyclic strength: tau/sigma'_v0 = alpha N_l^(-beta)
        const scalar alpha_;
        const scalar beta_;

        //- Threshold stress ratio: no generation for tau/sigma'_v0 below it
        const scalar thresholdRatio_;

        //- Equivalent cycles from the current pressure (true) or
        //- elapsed cycles (false)
        const bool equivalentCycles_;

        //- Shear-amplitude type: wave, computed, constant
        const word amplitudeType_;

        //- wave: pressure amplitude p0 [Pa] and wavenumber k [1/m];
        //- constant: the value [Pa]
        scalar p0_;
        scalar k_;
        scalar tauConst_;

        //- Generated pore pressure u_g [Pa] (state, written)
        volScalarField ug_;

        //- Cyclic shear amplitude of the last evaluation [Pa] (written)
        volScalarField tau_;

        //- computed: running envelope of sigma'_xz over the current
        //- period and the start time of the period
        scalarField sMax_;
        scalarField sMin_;
        scalar tPeriod_;

        //- computed: the excess pore pressure integrated over the current
        //- period and its duration, and the period mean of the last
        //- completed period (the residual pressure that counts the
        //- equivalent cycles; empty before the first period)
        scalarField pSum_;
        scalar tSum_;
        scalarField pBar_;


public:

    //- Runtime type information
    TypeName("poreGenerationModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            poreGenerationModel,
            dictionary,
            (const fvMesh& mesh, const dictionary& dict),
            (mesh, dict)
        );


    // Constructors

        //- Construct from the soil mesh and the poreGeneration dictionary
        poreGenerationModel(const fvMesh& mesh, const dictionary& dict);


    // Selectors

        //- Return the selected model from the 'poreGeneration'
        //- sub-dictionary of soilDict, or nullptr if absent
        static autoPtr<poreGenerationModel> New
        (
            const fvMesh& mesh,
            const dictionary& soilDict
        );


    //- Destructor
    virtual ~poreGenerationModel() = default;


    // Member Functions

        //- Generated pore pressure u_g [Pa] after N cycles of a point
        //- with the vertical effective stress sigmaV0 and the cycles to
        //- liquefaction Nl (the law)
        virtual scalar generated
        (
            const scalar sigmaV0,
            const scalar Nl,
            const scalar N
        ) const = 0;

        //- Inverse of the law: the cycles N at which the generated
        //- pressure equals u (0 <= u <= sigmaV0)
        virtual scalar cyclesFor
        (
            const scalar sigmaV0,
            const scalar Nl,
            const scalar u
        ) const = 0;

        //- Whether N is the equivalent number of cycles of the current
        //- pressure
        bool equivalentCycles() const noexcept { return equivalentCycles_; }

        //- Generation increment [Pa] of a point over the interval
        //- (tOld, t] for the current excess pore pressure u: the rate is
        //- the increment over dt (shared with the finite-difference
        //- reference of the verification)
        scalar increment
        (
            const scalar sigmaV0,
            const scalar tau,
            const scalar u,
            const scalar tOld,
            const scalar t
        ) const;

        //- Cycles to liquefaction for the cyclic shear amplitude tau and
        //- the vertical effective stress sigmaV0 (GREAT when no shear)
        scalar cyclesToLiquefaction
        (
            const scalar sigmaV0,
            const scalar tau
        ) const;

        //- Cyclic shear amplitude [Pa] of the wave/constant types at the
        //- depth z below the bed
        scalar amplitude(const scalar depth) const;

        //- Loading period [s]
        scalar period() const noexcept { return T_; }

        //- Start time of the loading [s]
        scalar startTime() const noexcept { return t0_; }

        //- Cycles elapsed at time t
        scalar cycles(const scalar t) const
        {
            return max(t - t0_, scalar(0))/T_;
        }

        //- Generated pore pressure field [Pa]
        const volScalarField& ug() const noexcept { return ug_; }

        //- Update the envelope of the solved shear stress sigma'_xz at
        //- time t (computed amplitude)
        void updateEnvelope(const scalarField& sigmaXZ, const scalar t);

        //- Generation rate f [Pa/s] per cell over the physical interval
        //- dt ending at time t, for the vertical effective stress sigmaV0,
        //- depth below the bed and current excess pore pressure of the
        //- cells; accumulates u_g
        tmp<scalarField> rate
        (
            const scalar t,
            const scalar dt,
            const scalarField& sigmaV0,
            const scalarField& depth,
            const scalarField& pExcess
        );

        //- Remap the state after the soil columns moved (see
        //- soilModel::remap)
        void remap
        (
            const labelList& c0,
            const labelList& c1,
            const scalarField& w
        );

        //- Write the state fields
        void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
