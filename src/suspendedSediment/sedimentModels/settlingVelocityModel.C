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

#include "settlingVelocityModel.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * Base class * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace suspendedSediment
{
    defineTypeNameAndDebug(settlingVelocityModel, 0);
    defineRunTimeSelectionTable(settlingVelocityModel, dictionary);
}
}


Foam::suspendedSediment::settlingVelocityModel::settlingVelocityModel
(
    const dictionary& dict
)
:
    d50_(dict.get<scalar>("d50")),
    s_(dict.get<scalar>("s"))
{}


Foam::autoPtr<Foam::suspendedSediment::settlingVelocityModel>
Foam::suspendedSediment::settlingVelocityModel::New(const dictionary& dict)
{
    const word modelType(dict.get<word>("settlingVelocityModel"));

    Info<< "Selecting settling velocity model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "settlingVelocityModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<settlingVelocityModel>(ctorPtr(dict));
}


Foam::scalar Foam::suspendedSediment::settlingVelocityModel::Dstar
(
    const scalar nu,
    const scalar g
) const
{
    return d50_*cbrt((s_ - 1)*g/sqr(nu));
}


// * * * * * * * * * * * * * * * * * Models * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace suspendedSediment
{
namespace settlingVelocityModels
{

//- Settling velocity given in sedimentProperties (entry ws)
class constant : public settlingVelocityModel
{
    const scalar ws_;

public:

    TypeName("constant");

    explicit constant(const dictionary& dict)
    :
        settlingVelocityModel(dict),
        ws_(dict.get<scalar>("ws"))
    {}

    scalar ws(const scalar, const scalar) const override
    {
        return ws_;
    }
};


//- Stokes law (van Rijn 1984b, Eq. 11)
class stokes : public settlingVelocityModel
{
public:

    TypeName("stokes");

    explicit stokes(const dictionary& dict)
    :
        settlingVelocityModel(dict)
    {}

    scalar ws(const scalar nu, const scalar g) const override
    {
        return (s_ - 1)*g*sqr(d50_)/(18*nu);
    }
};


//- van Rijn (1984b) Eqs. 11-13: Stokes, Zanke, and the square-root law
class vanRijn1984 : public settlingVelocityModel
{
public:

    TypeName("vanRijn1984");

    explicit vanRijn1984(const dictionary& dict)
    :
        settlingVelocityModel(dict)
    {}

    scalar ws(const scalar nu, const scalar g) const override
    {
        const scalar d = d50_;

        if (d <= 1e-4)
        {
            return (s_ - 1)*g*sqr(d)/(18*nu);
        }
        else if (d <= 1e-3)
        {
            return 10*nu/d*(sqrt(1 + 0.01*(s_ - 1)*g*pow3(d)/sqr(nu)) - 1);
        }

        return 1.1*sqrt((s_ - 1)*g*d);
    }
};


//- Soulsby (1997), as quoted by Liang, Cheng & Li (2005), Eq. 8
class soulsby : public settlingVelocityModel
{
public:

    TypeName("soulsby");

    explicit soulsby(const dictionary& dict)
    :
        settlingVelocityModel(dict)
    {}

    scalar ws(const scalar nu, const scalar g) const override
    {
        return nu/d50_*(sqrt(sqr(10.36) + 1.049*pow3(Dstar(nu, g))) - 10.36);
    }
};


defineTypeNameAndDebug(constant, 0);
defineTypeNameAndDebug(stokes, 0);
defineTypeNameAndDebug(vanRijn1984, 0);
defineTypeNameAndDebug(soulsby, 0);

addToRunTimeSelectionTable(settlingVelocityModel, constant, dictionary);
addToRunTimeSelectionTable(settlingVelocityModel, stokes, dictionary);
addToRunTimeSelectionTable(settlingVelocityModel, vanRijn1984, dictionary);
addToRunTimeSelectionTable(settlingVelocityModel, soulsby, dictionary);

} // End namespace settlingVelocityModels
} // End namespace suspendedSediment
} // End namespace Foam


// ************************************************************************* //
