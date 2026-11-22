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

#include "ibWallLaw.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ibWallLaw, 0);
    defineRunTimeSelectionTable(ibWallLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ibWallLaw::ibWallLaw(const dictionary& dict)
:
    kappa_(dict.getOrDefault<scalar>("kappa", 0.41)),
    E_(dict.getOrDefault<scalar>("E", 9.8)),
    ks_(dict.getOrDefault<scalar>("ks", 0))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ibWallLaw> Foam::ibWallLaw::New(const dictionary& dict)
{
    const word modelType(dict.get<word>("type"));

    Info<< "    Selecting immersed-boundary wall law " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "ibWallLaw",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<ibWallLaw>(ctorPtr(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::ibWallLaw::uTau
(
    const scalar Ut,
    const scalar y,
    const scalar nu
) const
{
    if (Ut < VSMALL || y < VSMALL)
    {
        return 0;
    }

    // Residual F(u) = u*u+(y u/nu, ks u/nu) - Ut, increasing in u
    auto residual = [&](const scalar u)
    {
        return u*uPlus(y*u/nu, ks_*u/nu) - Ut;
    };

    // Bracket: the laminar estimate is an upper bound for u+ >= y+ laws
    // only, so expand the bracket until the sign changes
    scalar uLo = 0;
    scalar fLo = -Ut;
    scalar uHi = max(Foam::sqrt(nu*Ut/y), Ut/30);
    scalar fHi = residual(uHi);

    for (label iter = 0; fHi < 0 && iter < 60; ++iter)
    {
        uLo = uHi;
        fLo = fHi;
        uHi *= 2;
        fHi = residual(uHi);
    }

    // Safeguarded regula falsi (Illinois variant)
    scalar u = uHi;
    int side = 0;

    for (label iter = 0; iter < 100; ++iter)
    {
        u = (uLo*fHi - uHi*fLo)/(fHi - fLo);
        const scalar f = residual(u);

        if (mag(f) < 1e-10*Ut || (uHi - uLo) < 1e-12*uHi)
        {
            break;
        }

        if (f*fHi > 0)
        {
            uHi = u;
            fHi = f;
            if (side == -1)
            {
                fLo *= 0.5;
            }
            side = -1;
        }
        else
        {
            uLo = u;
            fLo = f;
            if (side == 1)
            {
                fHi *= 0.5;
            }
            side = 1;
        }
    }

    return u;
}


void Foam::ibWallLaw::wallTurbulence
(
    const scalar uTau,
    const scalar y,
    const scalar nu,
    const scalar Cmu,
    const scalar beta1,
    const scalar Aplus,
    scalar& omega,
    scalar& nutLog
) const
{
    const scalar kappa = kappa_;
    const scalar sqrtCmu = Foam::sqrt(Cmu);

    // Roughness: shift of the mixing length and offset of the near-wall
    // omega solution that gives omega_w = S_r u_tau^2/nu at the wall
    scalar dy = 0;
    scalar ys = 0;

    if (uTau > VSMALL)
    {
        const scalar ksPlus = ks_*uTau/nu;
        dy = yPlusShift(ksPlus)*nu/uTau;
        ys = nu*Foam::sqrt(6/(beta1*omegaWallCoeff(ksPlus)))/uTau;
    }

    const scalar omegaVis = 6*nu/(beta1*sqr(y + ys));
    // Log layer: omega = sqrt(k)/(Cmu^0.25 kappa y), sqrt(k) = u_tau/Cmu^0.25
    const scalar omegaLog = uTau/(sqrtCmu*kappa*(y + dy));
    const scalar damping = sqr(1 - Foam::exp(-(y + dy)*uTau/(nu*Aplus)));

    omega = Foam::sqrt(sqr(omegaVis) + sqr(omegaLog));
    nutLog = kappa*uTau*(y + dy)*damping;
}


// ************************************************************************* //
