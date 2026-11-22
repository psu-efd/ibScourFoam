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

#include "vanDriestWallLaw.H"
#include "addToRunTimeSelectionTable.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vanDriestWallLaw, 0);
    addToRunTimeSelectionTable(ibWallLaw, vanDriestWallLaw, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::vanDriestWallLaw::dUPlus
(
    const scalar yPlus,
    const scalar shift
) const
{
    const scalar ys = yPlus + shift;
    const scalar lPlus = kappa_*ys*(1 - Foam::exp(-ys/A_));

    return 2/(1 + Foam::sqrt(1 + 4*sqr(lPlus)));
}


void Foam::vanDriestWallLaw::tabulate()
{
    // u+ is integrated in t = ln(1 + y+) (dy+ = e^t dt), where the
    // integrand is smooth from the viscous to the log region
    const label nSub = 8;

    table_.setSize(nT_*nShift_);

    for (label j = 0; j < nShift_; ++j)
    {
        const scalar shift = j*dShift_;

        auto f = [&](const scalar t)
        {
            return dUPlus(Foam::exp(t) - 1, shift)*Foam::exp(t);
        };

        scalar u = 0;
        table_[j*nT_] = 0;

        for (label i = 1; i < nT_; ++i)
        {
            const scalar t0 = (i - 1)*dT_;
            const scalar h = dT_/nSub;

            // Composite Simpson over [t0, t0 + dT]
            scalar s = f(t0) + f(t0 + dT_);
            for (label k = 1; k < nSub; ++k)
            {
                s += (k % 2 ? 4 : 2)*f(t0 + k*h);
            }
            u += s*h/3;

            table_[j*nT_ + i] = u;
        }
    }
}


Foam::scalar Foam::vanDriestWallLaw::lookup
(
    const scalar yPlus,
    const scalar shift
) const
{
    const scalar tMax = (nT_ - 1)*dT_;
    const scalar sMax = (nShift_ - 1)*dShift_;

    const scalar s = min(max(shift, scalar(0)), sMax);
    const label j0 = min(label(s/dShift_), nShift_ - 2);
    const scalar ws = s/dShift_ - j0;

    auto row = [&](const label j, const scalar t)
    {
        const label i0 = min(label(t/dT_), nT_ - 2);
        const scalar wt = t/dT_ - i0;
        const scalar* r = &table_[j*nT_];

        return (1 - wt)*r[i0] + wt*r[i0 + 1];
    };

    const scalar t = std::log1p(max(yPlus, scalar(0)));

    if (t <= tMax)
    {
        return (1 - ws)*row(j0, t) + ws*row(j0 + 1, t);
    }

    // Beyond the table: log-law asymptote from the last tabulated point
    const scalar yMax = std::expm1(tMax);
    const scalar uMax = (1 - ws)*row(j0, tMax) + ws*row(j0 + 1, tMax);

    return uMax + Foam::log((yPlus + s)/(yMax + s))/kappa_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vanDriestWallLaw::vanDriestWallLaw(const dictionary& dict)
:
    ibWallLaw(dict),
    A_(dict.getOrDefault<scalar>("A", 25)),
    nT_(801),
    dT_(std::log1p(1e7)/(nT_ - 1)),
    nShift_(0),
    dShift_(0.25),
    table_()
{
    nShift_ = label(yPlusShift(2000)/dShift_) + 2;

    tabulate();

    Info<< "    van Driest wall law: A " << A_ << ", ks " << ks_
        << " m, table " << nT_ << " x " << nShift_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::vanDriestWallLaw::uPlus
(
    const scalar yPlus,
    const scalar ksPlus
) const
{
    return lookup(yPlus, yPlusShift(ksPlus));
}


Foam::scalar Foam::vanDriestWallLaw::yPlusShift(const scalar ksPlus) const
{
    if (ksPlus <= 5)
    {
        return 0;
    }

    const scalar k = min(ksPlus, scalar(2000));

    return max(0.9*(Foam::sqrt(k) - k*Foam::exp(-k/6)), scalar(0));
}


Foam::scalar Foam::vanDriestWallLaw::omegaWallCoeff
(
    const scalar ksPlus
) const
{
    if (ksPlus < SMALL)
    {
        return VGREAT;
    }

    if (ksPlus < 20.2)
    {
        return Foam::pow3(40/ksPlus);
    }

    return 100/Foam::pow(ksPlus, 0.85);
}


// ************************************************************************* //
