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
    Foam::soilConstitutiveModels::linearElastic

Description
    Isotropic linear elasticity: sigma' = 2 mu eps + lambda tr(eps) I with
    mu = E / (2 (1 + nu)), lambda = nu E / ((1 + nu)(1 - 2 nu)).

    Optional liquefaction envelope ('tensionCutoff'), stateless (a
    nonlinear elastic envelope, not plasticity):
      - p' < cutoffStress (liquefied): the total effective stress
        sigma0 + sigma is the apex, -cutoffStress I, plus the residual
        fraction f ('cutoffStiffness') of the trial stress beyond it. Pure
        shear does not change p', so a liquefied cell carries no shear
        (a cap of the principal stresses alone leaves a no-tension solid
        that still resists shear through the compression in the other
        direction). The residual fraction keeps a liquefied layer from
        being singular in the quasi-static problem (f = 1e-4 to 1e-2).
      - optional friction limit, 'frictionAngle' phi: q <= M p'
        (M = 6 sin(phi)/(3 - sin(phi))), the deviator scaled at constant
        p' to M p' plus the fraction f of the excess. Experimental: on the
        cone the shear stiffness falls to f mu while the bulk modulus stays,
        the near-incompressible limit in which the segregated displacement
        scheme diverges (seen with f = 1e-4); the coupled solve is needed.

    Tangent operator: the segregated Biot scheme contracts with the ratio
    of the true stiffness to the implicit one, so a liquefied cell under the
    elastic implicit operator needs ~1/f correctors, and the fixed-stress
    split diverges when its compliance estimate is below the true one. With
    the envelope the stabilisation moduli (mu_, lambda_) follow the branch
    of each cell after every stress evaluation, and the implicit stiffness
    is set at the start of each step (beginStep) from the branches of the
    last correct. The tangent moduli: elastic; on the cone
    the bulk modulus with the shear modulus f mu; at the apex f times the
    elastic moduli plus the viscous eta/dt in the shear modulus. The
    implicit stiffness 2 mu + lambda is floored at 'tangentFloor' times the
    elastic one (default 0), as a cell that re-stiffens within a step sees
    an implicit operator softer than its stiffness.

    Viscous liquefied soil (stage B): in the liquefied cells the deviatoric
    stress 2 eta dev(eps_dot) is added, with the strain rate from the
    displacement gradient of the last committed step,
    eps_dot = symm(gradD - gradD_old)/(t - t_old). A structure in liquefied
    soil then moves at the rate its net weight drives against the viscous
    drag. The viscous stress is explicit in the segregated scheme (stable
    while eta/dt is well below the shear modulus). Apparent viscosities of
    liquefied sand are to be calibrated (sphere-sinking and dam-break tests,
    Sumer et al. 1999; Della Vecchia et al. 2019).

Usage
    \verbatim
    constitutive
    {
        type            linearElastic;
        E               5e7;            // [Pa]
        nu              0.3;
        tensionCutoff   yes;            // liquefaction envelope (default no)
        frictionAngle   30;             // [deg] q <= M p' (experimental)
        cutoffStress    0;              // [Pa] apex p' (default 0)
        cutoffStiffness 1e-3;           // residual fraction (default 0)
        tangentFloor    0;              // of the elastic impK (default 0)
        viscosity       100;            // liquefied [Pa s] (default 0)
    }
    \endverbatim

SourceFiles
    linearElastic.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_soilLinearElastic_H
#define Foam_soilLinearElastic_H

#include "soilConstitutiveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace soilConstitutiveModels
{

/*---------------------------------------------------------------------------*\
                        Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class linearElastic
:
    public soilConstitutiveModel
{
protected:

    // Protected Data

        //- Young's modulus [Pa]
        const scalar E_;

        //- Poisson's ratio
        const scalar nu_;

        //- Liquefaction envelope of the total effective stress and the
        //- mean effective stress of its apex, cutoffStress [Pa]
        const bool tensionCutoff_;
        const scalar cutoffStress_;

        //- Fraction of the stress beyond the cap that is kept
        const scalar cutoffStiffness_;

        //- Slope M of the friction limit q <= M p' (0: none)
        const scalar frictionM_;

        //- Viscosity of the liquefied (cut-off) soil [Pa s]
        const scalar viscosity_;

        //- Floor of the tangent implicit stiffness, fraction of the elastic
        const scalar tangentFloor_;


    // Protected Member Functions

        //- Tangent moduli mu_ and lambda_ (the stabilisation of the
        //- fixed-stress split) from the branches of the last correct, and
        //- with withImpK the implicit stiffness (only between steps)
        void setTangent(const bool withImpK);

        //- Limit sigma0 + sigma by the liquefaction envelope (sigma is the
        //- increment, changed in place) and add the viscous stress
        //- 2 eta devRate in the liquefied entries (devRate: deviatoric
        //- strain rate, empty without viscosity); the branch of each entry
        //- is returned in branch when it is not null
        void cutTension
        (
            const symmTensorField& sigma0,
            const symmTensorField& devRate,
            symmTensorField& sigma,
            labelList* branch
        ) const;


private:

    // Private Data

        //- Lame's parameters and the implicit stiffness (uniform)
        volScalarField mu_;
        volScalarField lambda_;
        volScalarField impK_;

        //- Displacement gradient of the last correct and of the last
        //- committed step, and the time of that step (viscosity only)
        volTensorField gradDLast_;
        volTensorField gradDOld_;
        scalar tOld_;

        //- The columns were remapped since the last commit
        bool rateReset_;

        //- Branch of the envelope per cell at the last correct (0 elastic,
        //- 1 friction limit, 2 liquefied)
        labelList branch_;

        //- The geostatic state has been set once
        bool geostaticSeen_;


public:

    //- Runtime type information
    TypeName("linearElastic");


    // Constructors

        //- Construct from the soil mesh and dictionary
        linearElastic(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~linearElastic() = default;


    // Member Functions

        virtual const volScalarField& mu() const { return mu_; }

        virtual const volScalarField& lambda() const { return lambda_; }

        virtual const volScalarField& impK() const { return impK_; }

        virtual scalar constrainedModulus() const
        {
            return E_*(1 - nu_)/((1 + nu_)*(1 - 2*nu_));
        }

        virtual void correct
        (
            const volTensorField& gradD,
            const volSymmTensorField& sigma0,
            volSymmTensorField& sigma
        );

        //- Commit the displacement gradient for the viscous strain rate
        virtual void updateState();

        //- Tangent moduli of the step from the branches of the envelope
        virtual void beginStep();

        //- The geostatic stress is set: elastic branch everywhere (first
        //- call)
        virtual void geostatic(const volSymmTensorField& sigma0);

        //- The columns moved: no strain rate for the next step
        virtual void remap
        (
            const labelList& c0,
            const labelList& c1,
            const scalarField& w
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace soilConstitutiveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
