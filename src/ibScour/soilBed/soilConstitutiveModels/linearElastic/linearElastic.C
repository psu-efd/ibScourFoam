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
\*---------------------------------------------------------------------------*/

#include "linearElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace soilConstitutiveModels
{
    defineTypeNameAndDebug(linearElastic, 0);
    addToRunTimeSelectionTable
    (
        soilConstitutiveModel,
        linearElastic,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::soilConstitutiveModels::linearElastic::linearElastic
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    soilConstitutiveModel(mesh, dict),
    E_(dict.get<scalar>("E")),
    nu_(dict.get<scalar>("nu")),
    tensionCutoff_(dict.getOrDefault<bool>("tensionCutoff", false)),
    cutoffStress_(dict.getOrDefault<scalar>("cutoffStress", 0)),
    cutoffStiffness_(dict.getOrDefault<scalar>("cutoffStiffness", 0)),
    frictionM_
    (
        tensionCutoff_ && dict.found("frictionAngle")
      ? 6*Foam::sin(degToRad(dict.get<scalar>("frictionAngle")))
       /(3 - Foam::sin(degToRad(dict.get<scalar>("frictionAngle"))))
      : 0
    ),
    viscosity_(dict.getOrDefault<scalar>("viscosity", 0)),
    tangentFloor_(dict.getOrDefault<scalar>("tangentFloor", 0)),
    mu_
    (
        IOobject("soilMu", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(dimPressure, E_/(2*(1 + nu_)))
    ),
    lambda_
    (
        IOobject("soilLambda", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(dimPressure, nu_*E_/((1 + nu_)*(1 - 2*nu_)))
    ),
    impK_
    (
        IOobject("soilImpK", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar
        (
            dimPressure,
            E_/(1 + nu_) + nu_*E_/((1 + nu_)*(1 - 2*nu_))
        )
    ),
    gradDLast_
    (
        IOobject("soilGradDLast", mesh.time().timeName(), mesh),
        mesh,
        dimensionedTensor(dimless, Zero)
    ),
    gradDOld_
    (
        IOobject("soilGradDOld", mesh.time().timeName(), mesh),
        mesh,
        dimensionedTensor(dimless, Zero)
    ),
    tOld_(mesh.time().value()),
    rateReset_(false),
    branch_(mesh.nCells(), 0),
    geostaticSeen_(false)
{
    if (viscosity_ > 0 && !tensionCutoff_)
    {
        FatalIOErrorInFunction(dict)
            << "viscosity acts in the cut-off cells: set tensionCutoff yes"
            << exit(FatalIOError);
    }

    Info<< "    Linear elastic soil: E " << E_ << " Pa, nu " << nu_
        << ", constrained modulus " << constrainedModulus() << " Pa";
    if (tensionCutoff_)
    {
        Info<< ", liquefaction envelope: apex p' " << cutoffStress_
            << " Pa, M " << frictionM_ << ", residual " << cutoffStiffness_
            << ", liquefied viscosity " << viscosity_ << " Pa s";
    }
    Info<< endl;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::soilConstitutiveModels::linearElastic::cutTension
(
    const symmTensorField& sigma0,
    const symmTensorField& devRate,
    symmTensorField& sigma,
    labelList* branch
) const
{
    const bool viscous = (devRate.size() == sigma.size());

    const scalar pApex = cutoffStress_;
    const scalar f = cutoffStiffness_;

    forAll(sigma, i)
    {
        const symmTensor sT(sigma0[i] + sigma[i]);

        // Mean effective stress (compression positive) and deviator
        const scalar p = -tr(sT)/3;
        const symmTensor sd(dev(sT));
        const scalar q = Foam::sqrt(1.5*magSqr(sd));
        const scalar qMax = frictionM_*(p - pApex);

        label b = 0;

        if (p < pApex)
        {
            b = 2;

            // Liquefied: the apex plus the residual fraction of the trial
            // stress beyond it, and the viscous stress
            sigma[i] = -pApex*I + f*(sT + pApex*I) - sigma0[i];

            if (viscous)
            {
                sigma[i] += 2*viscosity_*devRate[i];
            }
        }
        else if (frictionM_ > 0 && q > qMax)
        {
            b = 1;

            // Friction limit at constant p'
            sigma[i] = sd*(qMax + f*(q - qMax))/q - p*I - sigma0[i];
        }

        if (branch)
        {
            (*branch)[i] = b;
        }

        if (debug && !std::isfinite(mag(sigma[i])))
        {
            Pout<< "linearElastic: non-finite stress at " << i << " (of "
                << sigma.size() << "), branch " << b << ", trial " << sT
                << ", sigma0 " << sigma0[i] << ", p " << p << ", q " << q
                << endl;
        }
    }
}


void Foam::soilConstitutiveModels::linearElastic::setTangent
(
    const bool withImpK
)
{

    const scalar mu0 = E_/(2*(1 + nu_));
    const scalar lambda0 = nu_*E_/((1 + nu_)*(1 - 2*nu_));
    const scalar K0 = lambda0 + 2*mu0/3;
    const scalar f = cutoffStiffness_;
    const scalar dt = mesh_.time().value() - tOld_;
    const scalar muVisc =
        (viscosity_ > 0 && dt > SMALL) ? viscosity_/dt : scalar(0);
    const scalar impKMin = tangentFloor_*(lambda0 + 2*mu0);

    scalarField& mu = mu_.primitiveFieldRef();
    scalarField& lambda = lambda_.primitiveFieldRef();
    scalarField& impK = impK_.primitiveFieldRef();

    forAll(branch_, celli)
    {
        if (branch_[celli] == 1)
        {
            mu[celli] = f*mu0;
            lambda[celli] = K0 - 2*mu[celli]/3;
        }
        else if (branch_[celli] == 2)
        {
            mu[celli] = f*mu0 + muVisc;
            lambda[celli] = f*lambda0;
        }
        else
        {
            mu[celli] = mu0;
            lambda[celli] = lambda0;
        }
        if (withImpK)
        {
            impK[celli] = max(2*mu[celli] + lambda[celli], impKMin);
        }
    }

    mu_.correctBoundaryConditions();
    lambda_.correctBoundaryConditions();
    if (withImpK)
    {
        impK_.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::soilConstitutiveModels::linearElastic::correct
(
    const volTensorField& gradD,
    const volSymmTensorField& sigma0,
    volSymmTensorField& sigma
)
{
    // Elastic trial stress with the elastic moduli (the mu_ and lambda_
    // fields hold the tangent of the envelope for the stabilisation)
    const dimensionedScalar mu0(dimPressure, E_/(2*(1 + nu_)));
    const dimensionedScalar lambda0
    (
        dimPressure,
        nu_*E_/((1 + nu_)*(1 - 2*nu_))
    );
    sigma = mu0*twoSymm(gradD) + lambda0*(I*tr(gradD));

    if (tensionCutoff_)
    {
        // Deviatoric strain rate over the step (none in the step after a
        // remap)
        const scalar dt = mesh_.time().value() - tOld_;
        const bool viscous = (viscosity_ > 0 && dt > SMALL && !rateReset_);

        symmTensorField rate;
        if (viscosity_ > 0)
        {
            gradDLast_ = gradD;
        }
        if (viscous)
        {
            rate =
                dev
                (
                    symm(gradD.primitiveField() - gradDOld_.primitiveField())
                )/dt;
        }
        cutTension
        (
            sigma0.primitiveField(),
            rate,
            sigma.primitiveFieldRef(),
            &branch_
        );

        forAll(sigma.boundaryField(), patchi)
        {
            symmTensorField pRate;
            if (viscous)
            {
                pRate =
                    dev
                    (
                        symm
                        (
                            gradD.boundaryField()[patchi]
                          - gradDOld_.boundaryField()[patchi]
                        )
                    )/dt;
            }
            cutTension
            (
                sigma0.boundaryField()[patchi],
                pRate,
                sigma.boundaryFieldRef()[patchi],
                nullptr
            );
        }

        // The stabilisation of the fixed-stress split must not
        // under-estimate the compliance: its moduli follow the branches of
        // this evaluation (the implicit stiffness changes only between
        // steps, so that the implicit and explicit parts agree)
        setTangent(false);

        if (debug)
        {
            label nApex = 0;
            forAll(branch_, celli)
            {
                nApex += (branch_[celli] == 2);
            }
            Info<< "    linearElastic: " << returnReduce(nApex, sumOp<label>())
                << " liquefied cells, max |gradD| "
                << gMax(mag(gradD.primitiveField())) << ", max |sigma| "
                << gMax(mag(sigma.primitiveField())) << endl;
        }
    }
}


void Foam::soilConstitutiveModels::linearElastic::updateState()
{
    if (viscosity_ > 0)
    {
        gradDOld_ = gradDLast_;
        tOld_ = mesh_.time().value();
        rateReset_ = false;
    }
}


void Foam::soilConstitutiveModels::linearElastic::beginStep()
{
    if (tensionCutoff_)
    {
        setTangent(true);
    }
}


void Foam::soilConstitutiveModels::linearElastic::geostatic
(
    const volSymmTensorField& sigma0
)
{
    // Evaluations before the geostatic state (sigma0 = 0) put every cell at
    // the apex: start from the elastic branch (first call only; the soil
    // region sets the geostatic state again at every update)
    if (!geostaticSeen_)
    {
        geostaticSeen_ = true;
        branch_ = 0;

        if (tensionCutoff_)
        {
            setTangent(true);
        }
    }
}


void Foam::soilConstitutiveModels::linearElastic::remap
(
    const labelList& c0,
    const labelList& c1,
    const scalarField& w
)
{
    // The displacement was remapped: the committed gradient is stale, so
    // the next step has no viscous stress and restarts the strain rate
    rateReset_ = true;
}


// ************************************************************************* //
