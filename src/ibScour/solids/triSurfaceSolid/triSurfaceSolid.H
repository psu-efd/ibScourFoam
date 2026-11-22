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
    Foam::triSurfaceSolid

Description
    Solid described by a triangulated surface in constant/triSurface.

    For a closed surface the solid is its interior (inside/outside test).
    For an open surface, such as a bed surface spanning the domain, the
    fluid is on the side the face normals point to (or the opposite side
    with 'flipNormals true'). The normal returned is the face normal of the
    nearest triangle, oriented into the fluid; it does not depend on the
    position of the query point, which keeps it smooth near the surface.

    The surface is held on every processor; it is meant for structures of
    moderate size. Large evolving beds use a bed model instead (Stage 1D).

Usage
    \verbatim
    pier
    {
        type            triSurface;
        file            "pier.stl";
        closed          true;       // optional; default: from the surface
        flipNormals     false;      // open surfaces only
    }
    \endverbatim

SourceFiles
    triSurfaceSolid.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_triSurfaceSolid_H
#define Foam_triSurfaceSolid_H

#include "ibSolid.H"
#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class triSurfaceSolid Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceSolid
:
    public ibSolid
{
    // Private Data

        //- The surface with its search tree
        triSurfaceMesh surface_;

        //- Treat the surface as closed (inside = solid)
        bool closed_;

        //- Orientation of the face normals relative to the fluid (+1/-1)
        scalar normalSign_;


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
    TypeName("triSurface");


    // Constructors

        //- Construct from name, dictionary and registry
        triSurfaceSolid
        (
            const word& name,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~triSurfaceSolid() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
