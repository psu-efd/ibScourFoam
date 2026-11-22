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
    Foam::soilConstitutiveModels::pastorZienkiewicz

Description
    Generalised-plasticity law for sand of Pastor, Zienkiewicz & Chan
    (1990), the "Mark III" model, for the Biot soil models:
    contractive plastic shearing below the characteristic
    state generates pore pressure under poor drainage, the mechanism of
    residual liquefaction, without a yield surface.

    Triaxial invariants (compression positive): p' = -tr(sigma')/3,
    q = sqrt(3/2 s:s), eta = q/p'; volumetric strain epsilon_v = -tr(eps),
    deviatoric strain epsilon_s = sqrt(2/3 e:e).
    - elastic moduli K = K0 p'/p_ref, G = G0 p'/p_ref, with p' capped at
      modulusCap x the geostatic p' (the implicit stiffness of the
      segregated scheme is set at the cap);
    - plastic flow direction n_g = (d_g, 1)/sqrt(1 + d_g^2) in (p', q),
      dilatancy d_g = (1 + alpha_g)(M_g - eta); in unloading the
      volumetric component is -|d_g| (contractive);
    - loading direction n_f likewise with M_f, alpha_f; loading when
      n_f : dsigma'^e > 0 (elastic trial of the strain increment);
    - plastic modulus in loading H_L = H0 p' H_f max(H_v + H_s, 0) H_DM
      (the softening branch beyond the peak is floored at zero) with
      H_f = (1 - eta/eta_f)^4, eta_f = (1 + 1/alpha_f) M_f, H_v = 1 - eta/M_g,
      H_s = beta0 beta1 exp(-beta0 xi) (xi the accumulated deviatoric
      plastic strain), H_DM = (zeta_max/zeta)^gammaDM with
      zeta = p' (1 - alpha_f eta/((1 + alpha_f) M_f))^(-1/alpha_f);
    - plastic modulus in unloading H_U = Hu0 (M_g/eta_u)^gammaU for
      |M_g/eta_u| > 1, else Hu0 (eta_u the stress ratio at the reversal);
    - strain-driven explicit integration from the committed state of the
      step with sub-steps (dsigma' limited to 'substepFraction' of p'),
      dsigma' = D^e (deps - dLambda n_g), dLambda = (n_f : D^e deps)/
      (H + n_f : D^e n_g); p' is kept above 'minPressure'.

    The stress increment of the interface is the total effective stress
    minus the geostatic one; the elastic moduli of the implicit operator
    are those of the geostatic p' (set by geostatic()).

    Default parameters: loose Banding sand of Pastor et al. (1990).

Usage
    \verbatim
    constitutive
    {
        type            pastorZienkiewicz;
        K0              35e6;       // bulk modulus at p_ref [Pa]
        G0              52.5e6;     // shear modulus at p_ref [Pa]
        pRef            100e3;      // reference pressure [Pa]
        Mg              1.5;
        Mf              1.1;
        alphaF          0.45;
        alphaG          0.45;
        H0              350;
        beta0           4.2;
        beta1           0.2;
        gammaDM         2;
        Hu0             2e6;
        gammaU          2;
        minPressure     1000;       // p' floor [Pa]
        tensileStrength 0;          // p_t [Pa]: the law sees p' + p_t
        floorWidth      0;          // C1 transition to the floor [Pa]
                                    // (0: hard floor)
        modulusCap      1.5;        // moduli capped at 1.5 x geostatic p'
        substepFraction 0.02;       // |dsigma'| per sub-step / p'
        substepStrain   1e-4;       // |deps| per sub-step
        writeState      true;       // pPrime, qDev, xi fields
    }
    \endverbatim

SourceFiles
    pastorZienkiewicz.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_soilConstitutiveModels_pastorZienkiewicz_H
#define Foam_soilConstitutiveModels_pastorZienkiewicz_H

#include "soilConstitutiveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace soilConstitutiveModels
{

/*---------------------------------------------------------------------------*\
                     Class pastorZienkiewicz Declaration
\*---------------------------------------------------------------------------*/

class pastorZienkiewicz
:
    public soilConstitutiveModel
{
public:

    //- Point state of the law
    struct state
    {
        symmTensor sigma;   //- total effective stress (tension positive)
        symmTensor eps;     //- total strain
        symmTensor epsP;    //- plastic strain
        scalar xi;          //- accumulated deviatoric plastic strain
        scalar zetaMax;     //- maximum of the mobilised stress function
        scalar etaU;        //- stress ratio at the last reversal
        scalar loading;     //- 1 loading, -1 unloading, 0 at start
        scalar pCap;        //- p' at which the elastic moduli are capped

        state()
        :
            sigma(Zero), eps(Zero), epsP(Zero), xi(0), zetaMax(0),
            etaU(0), loading(0), pCap(GREAT)
        {}
    };


private:

    // Private Data

        //- Parameters
        const scalar K0_;
        const scalar G0_;
        const scalar pRef_;
        const scalar Mg_;
        const scalar Mf_;
        const scalar alphaF_;
        const scalar alphaG_;
        const scalar H0_;
        const scalar beta0_;
        const scalar beta1_;
        const scalar gammaDM_;
        const scalar Hu0_;
        const scalar gammaU_;
        const scalar pMin_;

        //- Tensile strength p_t [Pa]: the law acts on p' + p_t (a small
        //- apparent cohesion), so that the moduli, the stress ratio and
        //- the flow and loading directions stay defined and smooth at
        //- the free surface of a seabed, where p' -> 0; the floor
        //- minPressure then applies to p' + p_t
        const scalar pT_;

        //- Width [Pa] of the C1 transition of p' + p_t to the floor
        //- minPressure (the floor and the tension cut-off are applied
        //- through f(x) = (x + w)^2/(4w) for |x| < w, x = p - p_min,
        //- instead of max(x, 0): no kink for the Newton solve); 0 keeps
        //- the hard floor
        const scalar floorWidth_;

        //- Floored p' + p_t of a raw value (hard or C1 floor)
        scalar floored(const scalar pRaw) const;

        //- The elastic moduli grow with p' only up to modulusCap x the
        //- geostatic p' of the point, where the implicit stiffness of
        //- the segregated scheme is set: a compressive overshoot of a
        //- corrector iterate cannot stiffen the law beyond the operator
        //- (the fixed point then diverges)
        const scalar modulusCap_;

        const scalar substepFraction_;
        const scalar substepStrain_;
        const label maxSubsteps_;

        //- Freeze the loading direction of each cell over the correctors
        //- of a step (decided at the first one)
        const bool freezeDirection_;

        //- Moduli of the implicit operator (geostatic p')
        volScalarField mu_;
        volScalarField lambda_;
        volScalarField impK_;

        //- Committed states (cells and boundary faces) and the working
        //- ones of the current step
        List<state> cellOld_;
        List<state> cell_;
        List<List<state>> patchOld_;
        List<List<state>> patch_;

        //- Written diagnostics: p', q, xi
        volScalarField pPrime_;
        volScalarField qDev_;
        volScalarField xiField_;

        //- Loading direction of the step per cell and boundary face
        //- (0 undecided, 1 loading, -1 unloading): decided at the first
        //- corrector from the elastic trial and kept over the correctors
        //- (the segregated iteration otherwise flips between the loading
        //- and unloading branches)
        labelList cellDir_;
        List<labelList> patchDir_;

        //- Cells that yielded (plastic multiplier > 0) at the last call
        label nPlastic_;

        //- States initialised (from whatever sigma0 the first correct
        //- saw) and initialised from the geostatic stress (once: the soil
        //- region sets the geostatic state after every solve)
        bool initialised_;
        bool geostaticSeen_;


    // Private Member Functions

        //- Mobilised-stress function zeta of (p', eta)
        scalar zeta(const scalar p, const scalar eta) const;

        //- Integrate the committed state s0 to the strain eps (total)
        //- with the loading direction dir (decided from the first
        //- plastic sub-step when 0); returns the number of sub-steps
        //- that were plastic
        label integrate
        (
            const state& s0,
            const symmTensor& eps,
            state& s,
            label& dir
        ) const;

        //- Initialise the states from the geostatic stress
        void initialise(const volSymmTensorField& sigma0);


public:

    //- Runtime type information
    TypeName("pastorZienkiewicz");


    // Constructors

        //- Construct from the soil mesh and dictionary
        pastorZienkiewicz(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~pastorZienkiewicz() = default;


    // Member Functions

        virtual const volScalarField& mu() const { return mu_; }
        virtual const volScalarField& lambda() const { return lambda_; }
        virtual const volScalarField& impK() const { return impK_; }

        virtual scalar constrainedModulus() const
        {
            return K0_ + 4*G0_/3;
        }

        //- Bulk modulus K at the current p' of each cell (capped, floored)
        virtual tmp<volScalarField> stabilisationModulusField() const;

        virtual void geostatic(const volSymmTensorField& sigma0);

        virtual void correct
        (
            const volTensorField& gradD,
            const volSymmTensorField& sigma0,
            volSymmTensorField& sigma
        );

        virtual void updateState();

        virtual void beginStep();

        //- Reset the frozen loading directions
        virtual void unfreeze();

        virtual void remap
        (
            const labelList& c0,
            const labelList& c1,
            const scalarField& w
        );

        //- Integrate one point from a state along a strain path step
        //- (public for the point-level verification)
        label integratePoint
        (
            const state& s0,
            const symmTensor& eps,
            state& s
        ) const
        {
            label dir = 0;
            return integrate(s0, eps, s, dir);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace soilConstitutiveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
