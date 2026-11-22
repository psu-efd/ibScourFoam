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
    Foam::sphereSolid

Description
    Sphere defined by its centre and radius (e.g. a gravel particle or a
    verification body).

Usage
    \verbatim
    stone
    {
        type            sphere;
        centre          (0 0 0.01); // [m]
        radius          0.01;       // [m]
    }
    \endverbatim

SourceFiles
    sphereSolid.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sphereSolid_H
#define Foam_sphereSolid_H

#include "ibSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sphereSolid Declaration
\*---------------------------------------------------------------------------*/

class sphereSolid
:
    public ibSolid
{
    // Private Data

        //- Centre (body frame) [m]
        point centre_;

        //- Radius [m]
        scalar radius_;


protected:

    // Protected Member Functions

        //- Signed distance and unit normal (into the fluid), body frame
        virtual void bodyDistance
        (
            const pointField& bodyPoints,
            scalarField& distance,
            vectorField& normal
        ) const;


public:

    //- Runtime type information
    TypeName("sphere");


    // Constructors

        //- Construct from name, dictionary and registry
        sphereSolid
        (
            const word& name,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~sphereSolid() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
