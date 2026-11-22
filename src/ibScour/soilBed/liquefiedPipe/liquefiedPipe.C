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

#include "liquefiedPipe.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::liquefiedPipe::netForce(const scalar e) const
{
    // Area of the section below the bed: the circle below a chord at
    // height -zc above the centre, zc = r - e the centre elevation
    const scalar r = 0.5*D_;
    const scalar h = min(max(e - r, -r), r);
    const scalar below =
        sqr(r)*Foam::acos(-h/r) + h*Foam::sqrt(max(sqr(r) - sqr(h), 0));
    const scalar A = constant::mathematical::pi*sqr(r);

    return rhoW_*g_*((A - below) + sl_*below - sp_*A) + fExt_;
}


Foam::scalar Foam::liquefiedPipe::speed(const scalar f) const
{
    if (drag_ == "linear")
    {
        return f/dragCoeff_;
    }

    // sumer1999: F = 0.5 s_l rho_w a sqrt(nu_w D) w^1.5
    const scalar k = 0.5*sl_*rhoW_*dragCoeff_*Foam::sqrt(nuW_*D_);
    return Foam::pow(f/k, 2.0/3.0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquefiedPipe::liquefiedPipe
(
    const soilRegion& soil,
    const word& name,
    const dictionary& dict,
    const scalar diameter,
    const scalar initialBurial
)
:
    soil_(soil),
    D_(diameter > 0 ? diameter : dict.get<scalar>("diameter")),
    sp_(dict.get<scalar>("specificGravity")),
    sl_(dict.get<scalar>("liquefiedSpecificGravity")),
    period_(dict.get<scalar>("period")),
    onsetRatio_(dict.getOrDefault<scalar>("onsetRatio", 0.5)),
    stableSpeed_(dict.getOrDefault<scalar>("stableLayerSpeed", 0)),
    drag_(dict.getOrDefault<word>("drag", "sumer1999")),
    dragCoeff_(dict.getOrDefault<scalar>("dragCoefficient", 8e7)),
    freeFall_(dict.getOrDefault<scalar>("freeFallSpeed", 0.1)),
    rhoW_(dict.getOrDefault<scalar>("waterDensity", 1000)),
    g_(dict.getOrDefault<scalar>("gravity", 9.81)),
    nuW_(dict.getOrDefault<scalar>("waterViscosity", 1e-6)),
    col_(-1),
    e_
    (
        initialBurial >= 0
      ? initialBurial
      : dict.getOrDefault<scalar>("initialBurial", 0)
    ),
    w_(0),
    pSum_(),
    pMean_(),
    tPeriod_(0),
    havePeriod_(false),
    tStable_(GREAT),
    fExt_(0),
    ratio_(0),
    stableTop_(0)
{
    if (drag_ != "sumer1999" && drag_ != "linear")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown drag " << drag_ << ": sumer1999 or linear"
            << exit(FatalIOError);
    }

    // The column nearest x, found over all processors
    const scalar x = dict.get<scalar>("x");
    const vectorField Cf(soil_.floorPatch().faceCentres());
    scalar dLoc = GREAT;
    label iLoc = -1;
    forAll(Cf, i)
    {
        if (mag(Cf[i].x() - x) < dLoc)
        {
            dLoc = mag(Cf[i].x() - x);
            iLoc = i;
        }
    }
    const scalar dMin = returnReduce(dLoc, minOp<scalar>());
    const label owner = returnReduce
    (
        (iLoc >= 0 && dLoc <= dMin) ? Pstream::myProcNo() : Pstream::nProcs(),
        minOp<label>()
    );
    if (Pstream::myProcNo() == owner)
    {
        col_ = iLoc;
    }

    if (Pstream::master())
    {
        const fileName dir(soil_.mesh().time().globalPath()/"postProcessing");
        mkDir(dir);
        file_.reset(new OFstream(dir/("liquefiedPipe_" + name + ".dat")));
        file_() << "# t [s]  bottom depth e [m]  velocity w [m/s] (down)"
            << "  p~/sigma'_0 at the centre  stable-layer top depth [m]"
            << "  external force [N/m] (up)" << nl;
    }

    Info<< "    Liquefied pipe " << name << ": D " << D_ << " m, s_p " << sp_ << ", s_l "
        << sl_ << ", initial burial " << e_ << " m, onset ratio "
        << onsetRatio_ << ", stable layer " << stableSpeed_ << " m/s, drag "
        << drag_ << " " << dragCoeff_ << ", column at x = " << x
        << " on processor " << owner << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::liquefiedPipe::update(const scalar dt)
{
    const scalar t = soil_.mesh().time().value();

    // Period mean of the excess pore pressure of the column
    if (col_ >= 0)
    {
        const scalarField p(soil_.columnPorePressure(col_));
        if (pSum_.size() != p.size())
        {
            pSum_.setSize(p.size(), Zero);
        }
        pSum_ += p*dt;
    }
    tPeriod_ += dt;
    bool periodEnd = false;
    if (tPeriod_ >= period_*(1 - 1e-6))
    {
        if (col_ >= 0)
        {
            pMean_ = pSum_/tPeriod_;
            pSum_ = Zero;
        }
        tPeriod_ = 0;
        havePeriod_ = true;
        periodEnd = true;
    }

    // Mobilisation at the pipe centre and at the soil bottom, and the
    // layer thickness (on the owner; summed over the processors)
    scalar ratio = 0;
    scalar ratioBottom = 0;
    scalar layer = 0;
    if (col_ >= 0 && havePeriod_)
    {
        const scalar zBed = soil_.bedElevation(col_);
        const scalarField depth(zBed - soil_.columnZ(col_)());
        const scalar f0 = (1 + 2*soil_.k0())/3;
        scalarField r(depth.size());
        forAll(r, j)
        {
            const scalar s0 = (soil_.overburden() + soil_.gammaS()*depth[j])*f0;
            r[j] = pMean_[j]/max(s0, SMALL);
        }

        // Centre depth, at least the top cell (a pipe on the bed)
        const scalar dc = max(e_ - 0.5*D_, depth[0]);
        label j = 0;
        while (j < depth.size() - 2 && depth[j + 1] < dc)
        {
            ++j;
        }
        const scalar a =
            min(max((dc - depth[j])/max(depth[j + 1] - depth[j], SMALL), 0), 1);
        ratio = (1 - a)*r[j] + a*r[j + 1];
        ratioBottom = r.last();
        layer = zBed - soil_.bottom();
    }
    reduce(ratio, sumOp<scalar>());
    reduce(ratioBottom, sumOp<scalar>());
    reduce(layer, sumOp<scalar>());
    ratio_ = ratio;

    // Stable layer from the bottom once it has liquefied
    if (havePeriod_ && tStable_ >= GREAT && ratioBottom >= 1)
    {
        tStable_ = t;
    }
    stableTop_ =
        (tStable_ < GREAT)
      ? layer - max(stableSpeed_*(t - tStable_), scalar(0))
      : layer;

    // An unsupported pipe heavier than water (its bottom above the bed,
    // e.g. scoured from under it) drops onto the bed
    w_ = 0;
    if (e_ < -SMALL && sp_ > 1)
    {
        w_ = min(-e_/dt, freeFall_);
        e_ += dt*w_;
    }
    // Quasi-static motion of the mobilised pipe
    else if (havePeriod_ && ratio >= onsetRatio_ && e_ < stableTop_ - SMALL)
    {
        const scalar f = netForce(e_);
        w_ = (f < 0 ? 1 : -1)*speed(mag(f));
        scalar eNew = e_ + dt*w_;
        if (eNew > stableTop_)
        {
            eNew = stableTop_;
            w_ = 0;
        }
        if (eNew < 0 && w_ < 0)
        {
            // Floats onto the bed
            eNew = 0;
            w_ = 0;
        }
        e_ = eNew;
    }

    if (periodEnd && file_)
    {
        file_() << t << ' ' << e_ << ' ' << w_ << ' ' << ratio_ << ' '
            << stableTop_ << ' ' << fExt_ << endl;
    }
}


// ************************************************************************* //
