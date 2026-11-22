/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Xiaofeng Liu, Penn State University
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "referenceConcentrationModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * Base class * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace suspendedSediment
{
    defineTypeNameAndDebug(referenceConcentrationModel, 0);
    defineRunTimeSelectionTable(referenceConcentrationModel, dictionary);
}
}


Foam::suspendedSediment::referenceConcentrationModel::
referenceConcentrationModel(const dictionary& dict)
:
    d50_(dict.get<scalar>("d50")),
    s_(dict.get<scalar>("s")),
    delta_(dict.getOrDefault<scalar>("delta", -1)),
    skinFriction_(true),
    kN_(dict.getOrDefault<scalar>("kN", -1))
{
    const word shields
    (
        dict.getOrDefault<word>("shieldsNumber", "skinFriction")
    );

    if (shields == "total")
    {
        skinFriction_ = false;
    }
    else if (shields != "skinFriction")
    {
        FatalIOErrorInFunction(dict)
            << "shieldsNumber " << shields
            << ": use skinFriction (default) or total" << exit(FatalIOError);
    }
}


Foam::autoPtr<Foam::suspendedSediment::referenceConcentrationModel>
Foam::suspendedSediment::referenceConcentrationModel::New
(
    const dictionary& dict
)
{
    const word modelType
    (
        dict.getOrDefault<word>("referenceConcentrationModel", "zysermanFredsoe")
    );

    Info<< "Selecting reference concentration model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "referenceConcentrationModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<referenceConcentrationModel>(ctorPtr(dict));
}


Foam::tmp<Foam::scalarField>
Foam::suspendedSediment::referenceConcentrationModel::uStarSkin
(
    const scalarField& uStar,
    const scalarField&,
    const scalarField&,
    const scalar
) const
{
    return tmp<scalarField>::New(uStar);
}


Foam::tmp<Foam::scalarField>
Foam::suspendedSediment::referenceConcentrationModel::theta
(
    const scalarField& uStar,
    const scalar g
) const
{
    return sqr(uStar)/((s_ - 1)*g*d50_);
}


// * * * * * * * * * * * * * * * * * Models * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace suspendedSediment
{
namespace referenceConcentrationModels
{

//- Zyserman & Fredsoe (1994), at 2 d50 unless delta is given
class zysermanFredsoe : public referenceConcentrationModel
{
public:

    TypeName("zysermanFredsoe");

    explicit zysermanFredsoe(const dictionary& dict)
    :
        referenceConcentrationModel(dict)
    {
        if (delta_ <= 0)
        {
            delta_ = 2*d50_;
        }
        if (kN_ <= 0)
        {
            kN_ = 2.5*d50_;     // Engelund & Hansen (1967), ZF p. 1023
        }
        if (skinFriction_)
        {
            Info<< "    Shields number from the skin friction of"
                << " Zyserman & Fredsoe (1994) Eqs. 9-10, kN = " << kN_
                << " m" << endl;
        }
    }

    //- u*' of ZF Eqs. 9-10 with g D I = u*^2: solves
    //  V = u*' [6 + 2.5 ln(D u*'^2/(u*^2 kN))] by bisection; u*' <= u*
    tmp<scalarField> uStarSkin
    (
        const scalarField& uStar,
        const scalarField& V,
        const scalarField& D,
        const scalar
    ) const override
    {
        auto tus = tmp<scalarField>::New(uStar);
        scalarField& us = tus.ref();

        forAll(us, facei)
        {
            const scalar u = uStar[facei];
            const scalar Vf = V[facei];
            const scalar Df = D[facei];

            if (u <= SMALL || Df <= kN_)
            {
                continue;       // no flow or no depth: total u*
            }
            if (Vf <= 0)
            {
                us[facei] = 0;
                continue;
            }

            // F(x) = x [6 + 2.5 ln(D x^2/(u^2 kN))] - V, increasing where
            // the bracket is positive
            auto F = [&](const scalar x)
            {
                return x*(6 + 2.5*log(Df*sqr(x)/(sqr(u)*kN_))) - Vf;
            };

            if (F(u) <= 0)
            {
                continue;       // bed not rougher than kN: u*' = u*
            }

            // bracket = 0 at lo (F(lo) = -V < 0)
            scalar lo = u*sqrt(kN_/Df)*exp(-1.2);
            scalar hi = u;
            for (label iter = 0; iter < 60; ++iter)
            {
                const scalar mid = 0.5*(lo + hi);
                (F(mid) > 0 ? hi : lo) = mid;
            }
            us[facei] = 0.5*(lo + hi);
        }

        return tus;
    }

    tmp<scalarField> Cb
    (
        const scalarField& uStar,
        const scalar,
        const scalar g
    ) const override
    {
        const scalarField x(pow(max(theta(uStar, g) - 0.045, scalar(0)), 1.75));
        return 0.331*x/(1 + 0.72*x);
    }
};


//- van Rijn (1984b) Eq. 38 at the reference level a = delta
class vanRijn1984 : public referenceConcentrationModel
{
    //- Critical Shields parameter, van Rijn (1984a) Fig. 1
    static scalar thetaCr(const scalar Ds)
    {
        if (Ds <= 4) return 0.24/Ds;
        if (Ds <= 10) return 0.14*pow(Ds, -0.64);
        if (Ds <= 20) return 0.04*pow(Ds, -0.10);
        if (Ds <= 150) return 0.013*pow(Ds, 0.29);
        return 0.055;
    }

public:

    TypeName("vanRijn1984");

    explicit vanRijn1984(const dictionary& dict)
    :
        referenceConcentrationModel(dict)
    {
        if (delta_ <= 0)
        {
            FatalIOErrorInFunction(dict)
                << "vanRijn1984 needs the reference level delta (a = ks or"
                << " half the bed-form height, a >= 0.01 h)" << exit(FatalIOError);
        }
        if (skinFriction_ && kN_ <= 0)
        {
            if (!dict.found("d90"))
            {
                FatalIOErrorInFunction(dict)
                    << "vanRijn1984 with shieldsNumber skinFriction needs"
                    << " d90 (kN = 3 d90) or kN" << exit(FatalIOError);
            }
            kN_ = 3*dict.get<scalar>("d90");
        }
        if (skinFriction_)
        {
            Info<< "    Shields number from the grain shear velocity of"
                << " van Rijn (1984), C' = 18 log10(12 D/kN), kN = " << kN_
                << " m" << endl;
        }
    }

    //- u*' = g^(1/2) V/C', C' = 18 log10(12 D/kN) (R_b = D); u*' <= u*
    tmp<scalarField> uStarSkin
    (
        const scalarField& uStar,
        const scalarField& V,
        const scalarField& D,
        const scalar g
    ) const override
    {
        auto tus = tmp<scalarField>::New(uStar);
        scalarField& us = tus.ref();

        forAll(us, facei)
        {
            const scalar Cp = 18*log10(12*D[facei]/kN_);
            if (Cp > 0)
            {
                us[facei] = min(sqrt(g)*max(V[facei], scalar(0))/Cp, uStar[facei]);
            }
        }

        return tus;
    }

    tmp<scalarField> Cb
    (
        const scalarField& uStar,
        const scalar nu,
        const scalar g
    ) const override
    {
        const scalar Ds = d50_*cbrt((s_ - 1)*g/sqr(nu));
        const scalar uCr2 = thetaCr(Ds)*(s_ - 1)*g*d50_;
        const scalarField T(max((sqr(uStar) - uCr2)/uCr2, scalar(0)));
        return 0.015*d50_/delta_*pow(T, 1.5)/pow(Ds, 0.3);
    }
};


//- Reference concentration given (entry Cb); delta required
class constant : public referenceConcentrationModel
{
    const scalar Cb_;

public:

    TypeName("constant");

    explicit constant(const dictionary& dict)
    :
        referenceConcentrationModel(dict),
        Cb_(dict.get<scalar>("Cb"))
    {
        skinFriction_ = false;
        if (delta_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "constant reference concentration needs delta"
                << exit(FatalIOError);
        }
    }

    bool needsShear() const override
    {
        return false;
    }

    tmp<scalarField> Cb
    (
        const scalarField& uStar,
        const scalar,
        const scalar
    ) const override
    {
        return tmp<scalarField>::New(uStar.size(), Cb_);
    }
};


defineTypeNameAndDebug(zysermanFredsoe, 0);
defineTypeNameAndDebug(vanRijn1984, 0);
defineTypeNameAndDebug(constant, 0);

addToRunTimeSelectionTable(referenceConcentrationModel, zysermanFredsoe, dictionary);
addToRunTimeSelectionTable(referenceConcentrationModel, vanRijn1984, dictionary);
addToRunTimeSelectionTable(referenceConcentrationModel, constant, dictionary);

} // End namespace referenceConcentrationModels
} // End namespace suspendedSediment
} // End namespace Foam


// ************************************************************************* //
