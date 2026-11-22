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
    Foam::roughLogLawWallLaw

Description
    Log law with equivalent sand roughness ks, using the roughness
    function of OpenFOAM's nutkRoughWallFunction (Cebeci & Bradshaw):

        u+ = ln(E y+ / (1 + Cs ks+))/kappa   for y+ > y+_lam
        u+ = y+                              otherwise

    where y+_lam is the intersection of the two laws. For ks = 0 it reduces
    to the standard two-layer law.

Usage
    \verbatim
    wallLaw
    {
        type            roughLogLaw;
        ks              0.0005;     // equivalent sand roughness [m]
        Cs              0.5;        // roughness constant, optional
    }
    \endverbatim

SourceFiles
    roughLogLawWallLaw.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_roughLogLawWallLaw_H
#define Foam_roughLogLawWallLaw_H

#include "ibWallLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class roughLogLawWallLaw Declaration
\*---------------------------------------------------------------------------*/

class roughLogLawWallLaw
:
    public ibWallLaw
{
    // Private Data

        //- Roughness constant
        scalar Cs_;


public:

    //- Runtime type information
    TypeName("roughLogLaw");


    // Constructors

        //- Construct from dictionary
        explicit roughLogLawWallLaw(const dictionary& dict);


    //- Destructor
    virtual ~roughLogLawWallLaw() = default;


    // Member Functions

        //- Dimensionless velocity u+ at y+ for roughness ks+
        virtual scalar uPlus(const scalar yPlus, const scalar ksPlus) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
