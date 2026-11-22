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
    Foam::spaldingWallLaw

Description
    Spalding's (1961) single formula for a smooth wall, valid from the
    viscous sublayer through the log layer:

        y+ = u+ + exp(-kappa B) [exp(kappa u+) - 1 - kappa u+
             - (kappa u+)^2/2 - (kappa u+)^3/6],    B = ln(E)/kappa

    It adapts continuously to the y+ of the sampling point, as the
    y+-adaptive wall function of Xu & Liu (2021, IJNMF 93:1929-1946).
    u+(y+) is obtained by Newton iteration on the formula.

Usage
    \verbatim
    wallLaw
    {
        type            spalding;
        kappa           0.41;       // optional
        E               9.8;        // optional
    }
    \endverbatim

SourceFiles
    spaldingWallLaw.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_spaldingWallLaw_H
#define Foam_spaldingWallLaw_H

#include "ibWallLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class spaldingWallLaw Declaration
\*---------------------------------------------------------------------------*/

class spaldingWallLaw
:
    public ibWallLaw
{

public:

    //- Runtime type information
    TypeName("spalding");


    // Constructors

        //- Construct from dictionary
        explicit spaldingWallLaw(const dictionary& dict);


    //- Destructor
    virtual ~spaldingWallLaw() = default;


    // Member Functions

        //- Dimensionless velocity u+ at y+ for roughness ks+
        virtual scalar uPlus(const scalar yPlus, const scalar ksPlus) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
