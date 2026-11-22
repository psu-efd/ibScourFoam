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

\*---------------------------------------------------------------------------*/

#include "poreGenerationModel.H"
#include "soilModel.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poreGenerationModel, 0);
    defineRunTimeSelectionTable(poreGenerationModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poreGenerationModel::poreGenerationModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    T_(dict.get<scalar>("period")),
    t0_(dict.getOrDefault<scalar>("startTime", 0)),
    alpha_(dict.getOrDefault<scalar>("alpha", 0.246)),
    beta_(dict.getOrDefault<scalar>("beta", 0.165)),
    thresholdRatio_(dict.getOrDefault<scalar>("thresholdRatio", 0)),
    equivalentCycles_(dict.getOrDefault<bool>("equivalentCycles", true)),
    amplitudeType_(dict.subDict("shearAmplitude").get<word>("type")),
    p0_(0),
    k_(0),
    tauConst_(0),
    ug_
    (
        IOobject
        (
            "poreGenerated",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPressure, Zero)
    ),
    tau_
    (
        IOobject
        (
            "cyclicShear",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPressure, Zero)
    ),
    sMax_(),
    sMin_(),
    tPeriod_(t0_),
    pSum_(),
    tSum_(0),
    pBar_()
{
    const dictionary& ampDict = dict.subDict("shearAmplitude");

    if (amplitudeType_ == "wave")
    {
        p0_ = ampDict.get<scalar>("amplitude");
        k_ = 2*constant::mathematical::pi/ampDict.get<scalar>("wavelength");
    }
    else if (amplitudeType_ == "constant")
    {
        tauConst_ = ampDict.get<scalar>("value");
    }
    else if (amplitudeType_ == "computed")
    {
        sMax_.resize(mesh.nCells(), -GREAT);
        sMin_.resize(mesh.nCells(), GREAT);
        pSum_.resize(mesh.nCells(), Zero);
    }
    else
    {
        FatalIOErrorInFunction(ampDict)
            << "Unknown shearAmplitude type " << amplitudeType_
            << ": wave, computed or constant" << exit(FatalIOError);
    }

    Info<< "    Pore-pressure generation " << type() << ": period " << T_
        << " s from t = " << t0_ << " s, alpha " << alpha_ << ", beta "
        << beta_ << ", threshold ratio " << thresholdRatio_
        << (equivalentCycles_ ? ", equivalent" : ", elapsed")
        << " cycles, shear amplitude " << amplitudeType_;
    if (amplitudeType_ == "wave")
    {
        Info<< " (p0 " << p0_ << " Pa, k " << k_ << " 1/m)";
    }
    Info<< endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::poreGenerationModel> Foam::poreGenerationModel::New
(
    const fvMesh& mesh,
    const dictionary& soilDict
)
{
    if (!soilDict.found("poreGeneration"))
    {
        return nullptr;
    }

    const dictionary& dict = soilDict.subDict("poreGeneration");
    const word modelType(dict.get<word>("type"));

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "poreGeneration",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<poreGenerationModel>(ctorPtr(mesh, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::poreGenerationModel::cyclesToLiquefaction
(
    const scalar sigmaV0,
    const scalar tau
) const
{
    if (tau <= 0 || sigmaV0 <= 0 || tau < thresholdRatio_*sigmaV0)
    {
        return GREAT;
    }

    // tau/sigma'_v0 = alpha N_l^(-beta)
    return Foam::pow(tau/(alpha_*sigmaV0), -1/beta_);
}


Foam::scalar Foam::poreGenerationModel::amplitude(const scalar depth) const
{
    if (amplitudeType_ == "wave")
    {
        return p0_*k_*depth*Foam::exp(-k_*depth);
    }
    else if (amplitudeType_ == "constant")
    {
        return tauConst_;
    }

    return 0;
}


void Foam::poreGenerationModel::updateEnvelope
(
    const scalarField& sigmaXZ,
    const scalar t
)
{
    if (amplitudeType_ != "computed")
    {
        return;
    }

    if (t - tPeriod_ >= T_*(1 - SMALL))
    {
        // Period completed: its amplitude, then a new envelope
        scalarField& tau = tau_.primitiveFieldRef();
        forAll(tau, celli)
        {
            tau[celli] =
                (sMax_[celli] > sMin_[celli])
              ? 0.5*(sMax_[celli] - sMin_[celli])
              : 0;
        }
        sMax_ = -GREAT;
        sMin_ = GREAT;
        tPeriod_ += T_*std::floor((t - tPeriod_)/T_ + SMALL);
    }

    sMax_ = max(sMax_, sigmaXZ);
    sMin_ = min(sMin_, sigmaXZ);
}


Foam::scalar Foam::poreGenerationModel::increment
(
    const scalar sigmaV0,
    const scalar tau,
    const scalar u,
    const scalar tOld,
    const scalar t
) const
{
    if (t <= t0_ || sigmaV0 <= 0)
    {
        return 0;
    }

    const scalar Nl = cyclesToLiquefaction(sigmaV0, tau);

    if (Nl >= GREAT)
    {
        return 0;
    }

    const scalar uc = min(max(u, scalar(0)), sigmaV0);
    const scalar dN = (t - max(tOld, t0_))/T_;

    scalar N0, uOld;
    if (equivalentCycles_)
    {
        // Cycles that correspond to the current pressure, at most the
        // cycles elapsed: drainage lowers the count, but pressure fed by
        // diffusion from a liquefied neighbour is not loading (without
        // the bound the front runs away with the grid)
        N0 = min(cyclesFor(sigmaV0, Nl, uc), cycles(max(tOld, t0_)));
        uOld = uc;
    }
    else
    {
        N0 = cycles(max(tOld, t0_));
        uOld = generated(sigmaV0, Nl, N0);
    }

    const scalar du = generated(sigmaV0, Nl, N0 + dN) - uOld;

    // No generation past liquefaction
    return min(max(du, scalar(0)), sigmaV0 - uc);
}


Foam::tmp<Foam::scalarField> Foam::poreGenerationModel::rate
(
    const scalar t,
    const scalar dt,
    const scalarField& sigmaV0,
    const scalarField& depth,
    const scalarField& pExcess
)
{
    scalarField& ug = ug_.primitiveFieldRef();
    scalarField& tau = tau_.primitiveFieldRef();

    auto tf = tmp<scalarField>::New(ug.size(), Zero);
    scalarField& f = tf.ref();

    if (t <= t0_ || dt <= 0)
    {
        return tf;
    }

    // Computed amplitude: the period mean of the excess pore pressure
    // counts the equivalent cycles (the instantaneous value carries the
    // momentary component of the cyclic loading)
    const bool periodMean = (amplitudeType_ == "computed");
    if (periodMean)
    {
        pSum_ += pExcess*dt;
        tSum_ += dt;
        if (tSum_ >= T_*(1 - SMALL))
        {
            pBar_ = pSum_/tSum_;
            pSum_ = Zero;
            tSum_ = 0;
        }
    }
    const scalarField& u = (periodMean && pBar_.size()) ? pBar_ : pExcess;

    forAll(f, celli)
    {
        if (amplitudeType_ != "computed")
        {
            tau[celli] = amplitude(depth[celli]);
        }

        const scalar du =
            increment(sigmaV0[celli], tau[celli], u[celli], t - dt, t);

        f[celli] = du/dt;
        ug[celli] += du;
    }

    if (debug)
    {
        label nGen = 0;
        forAll(f, celli)
        {
            if (f[celli] > 0) ++nGen;
        }
        Info<< "    Pore-pressure generation: "
            << returnReduce(nGen, sumOp<label>()) << " cells, max f "
            << gMax(f) << " Pa/s, max p_e " << gMax(pExcess) << " Pa"
            << endl;
    }

    return tf;
}


void Foam::poreGenerationModel::remap
(
    const labelList& c0,
    const labelList& c1,
    const scalarField& w
)
{
    ug_.primitiveFieldRef() =
        soilModel::remapField(ug_.primitiveField(), c0, c1, w);
    tau_.primitiveFieldRef() =
        soilModel::remapField(tau_.primitiveField(), c0, c1, w);

    if (amplitudeType_ == "computed")
    {
        sMax_ = soilModel::remapField(sMax_, c0, c1, w);
        sMin_ = soilModel::remapField(sMin_, c0, c1, w);
        pSum_ = soilModel::remapField(pSum_, c0, c1, w);
        if (pBar_.size())
        {
            pBar_ = soilModel::remapField(pBar_, c0, c1, w);
        }
    }
}


void Foam::poreGenerationModel::write() const
{
    ug_.write();
    tau_.write();
}


// ************************************************************************* //
