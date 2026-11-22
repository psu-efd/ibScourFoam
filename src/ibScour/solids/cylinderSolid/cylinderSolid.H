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
    Foam::cylinderSolid

Description
    Circular cylinder defined by two axis points and a radius. With
    'finite false' (default) the cylinder extends infinitely along its axis
    (e.g. a pipeline); with 'finite true' it is capped at the two points.

Usage
    \verbatim
    pipe
    {
        type            cylinder;
        point1          (0 -1 0);   // axis start [m]
        point2          (0  1 0);   // axis end [m]
        radius          0.025;      // [m]
        finite          false;      // optional, default false
    }
    \endverbatim

SourceFiles
    cylinderSolid.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cylinderSolid_H
#define Foam_cylinderSolid_H

#include "ibSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cylinderSolid Declaration
\*---------------------------------------------------------------------------*/

class cylinderSolid
:
    public ibSolid
{
    // Private Data

        //- Axis start point (body frame) [m]
        point point1_;

        //- Unit axis direction
        vector axis_;

        //- Axis length [m]
        scalar length_;

        //- Radius [m]
        scalar radius_;

        //- Capped (true) or infinite (false) cylinder
        bool finite_;


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
    TypeName("cylinder");


    // Constructors

        //- Construct from name, dictionary and registry
        cylinderSolid
        (
            const word& name,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~cylinderSolid() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
