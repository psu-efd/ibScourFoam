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
    Foam::planeSolid

Description
    Half-space solid bounded by a plane: the fluid is on the side the
    normal points to. Used for flat beds and verification cases.

Usage
    \verbatim
    bed
    {
        type            plane;
        point           (0 0 0);    // point on the plane [m]
        normal          (0 0 1);    // points into the fluid
    }
    \endverbatim

SourceFiles
    planeSolid.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_planeSolid_H
#define Foam_planeSolid_H

#include "ibSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class planeSolid Declaration
\*---------------------------------------------------------------------------*/

class planeSolid
:
    public ibSolid
{
    // Private Data

        //- Point on the plane (body frame) [m]
        point point_;

        //- Unit normal, pointing into the fluid
        vector normal_;


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
    TypeName("plane");


    // Constructors

        //- Construct from name, dictionary and registry
        planeSolid
        (
            const word& name,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~planeSolid() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
