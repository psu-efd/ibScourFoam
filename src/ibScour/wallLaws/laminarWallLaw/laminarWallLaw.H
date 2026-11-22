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
    Foam::laminarWallLaw

Description
    Linear (viscous sublayer / laminar) wall law, u+ = y+. Exact for
    laminar flow over a plane wall when the sampling point lies in the
    linear part of the profile.

Usage
    \verbatim
    wallLaw
    {
        type            laminar;
    }
    \endverbatim

SourceFiles
    laminarWallLaw.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_laminarWallLaw_H
#define Foam_laminarWallLaw_H

#include "ibWallLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class laminarWallLaw Declaration
\*---------------------------------------------------------------------------*/

class laminarWallLaw
:
    public ibWallLaw
{

public:

    //- Runtime type information
    TypeName("laminar");


    // Constructors

        //- Construct from dictionary
        explicit laminarWallLaw(const dictionary& dict);


    //- Destructor
    virtual ~laminarWallLaw() = default;


    // Member Functions

        //- Dimensionless velocity u+ at y+ for roughness ks+
        virtual scalar uPlus(const scalar yPlus, const scalar ksPlus) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
