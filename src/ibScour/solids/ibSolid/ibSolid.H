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
    Foam::ibSolid

Description
    Abstract base class of a solid immersed in the fluid mesh.

    A solid provides, for any set of points, the signed distance to its
    surface (positive in the fluid, negative inside the solid) and the unit
    surface normal at the nearest surface point, pointing into the fluid.

    Every solid carries a rigid-body state: the position of its reference
    point, its orientation, and its linear and angular velocity. Derived
    classes describe the geometry in the body frame; the base class maps
    points and normals between the body and the global frame. A static
    solid has the identity transform and zero velocity. Motion solvers
    (later stages) update the state through setState().

Usage
    \verbatim
    <name>
    {
        type            plane;      // plane, cylinder, sphere, triSurface

        // Optional rigid-body state (default: at rest, identity)
        origin          (0 0 0);    // body reference point [m]
        velocity        (0 0 0);    // [m/s]
        omega           (0 0 0);    // angular velocity [rad/s]

        // Optional prescribed translation: displacement of the reference
        // point from 'origin' as a Function1 of time (e.g. a lifted gate);
        // the velocity follows from the displacement over the time step,
        // and a restart reproduces the position
        displacement    table ((0 (0 0 0)) (0.23 (0 0 0.4)));

        // type-specific entries
    }
    \endverbatim

SourceFiles
    ibSolid.C
    ibSolidNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ibSolid_H
#define Foam_ibSolid_H

#include "dictionary.H"
#include "pointField.H"
#include "tensor.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "objectRegistry.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ibSolid Declaration
\*---------------------------------------------------------------------------*/

class ibSolid
{
    // Private Data

        //- Name of the solid
        word name_;

        //- Body reference point in the global frame [m]
        point origin_;

        //- Rotation from the body frame to the global frame
        tensor R_;

        //- Linear velocity of the reference point [m/s]
        vector velocity_;

        //- Angular velocity [rad/s]
        vector omega_;

        //- Reference point at zero displacement [m]
        point origin0_;

        //- Optional prescribed displacement of the reference point [m]
        autoPtr<Function1<vector>> displacement_;


protected:

    // Protected Data

        //- Construction dictionary
        dictionary dict_;


    // Protected Member Functions

        //- Signed distance and outward (into fluid) unit normal of points
        //- given in the body frame
        virtual void bodyDistance
        (
            const pointField& bodyPoints,
            scalarField& distance,
            vectorField& normal
        ) const = 0;


public:

    //- Runtime type information
    TypeName("ibSolid");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            ibSolid,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const objectRegistry& obr
            ),
            (name, dict, obr)
        );


    // Constructors

        //- Construct from name and dictionary
        ibSolid(const word& name, const dictionary& dict);


    // Selectors

        //- Return a reference to the selected solid
        static autoPtr<ibSolid> New
        (
            const word& name,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~ibSolid() = default;


    // Member Functions

        // Access

            //- Name of the solid
            const word& name() const noexcept
            {
                return name_;
            }

            //- Body reference point [m]
            const point& origin() const noexcept
            {
                return origin_;
            }

            //- Rotation from the body frame to the global frame
            const tensor& R() const noexcept
            {
                return R_;
            }

            //- Linear velocity [m/s]
            const vector& velocity() const noexcept
            {
                return velocity_;
            }

            //- Angular velocity [rad/s]
            const vector& omega() const noexcept
            {
                return omega_;
            }

            //- True if the solid moves (non-zero velocity or rotation rate)
            bool moving() const;

            //- True if the solid follows a prescribed displacement
            bool prescribed() const noexcept
            {
                return bool(displacement_);
            }

            //- Apply the prescribed displacement at time t (step dt): set
            //- the origin and the velocity; true if the origin changed
            bool applyPrescribedMotion(const scalar t, const scalar dt);


        // Geometry

            //- Signed distance (> 0 in fluid) and unit normal into the
            //- fluid, for points in the global frame
            void distance
            (
                const pointField& points,
                scalarField& distance,
                vectorField& normal
            ) const;

            //- Signed distance and normal for the solid at the given
            //- rigid-body state (origin, R) instead of its current one
            void distance
            (
                const pointField& points,
                scalarField& distance,
                vectorField& normal,
                const point& origin,
                const tensor& R
            ) const;

            //- Rigid-body velocity of points on the solid [m/s]
            tmp<vectorField> wallVelocity(const pointField& points) const;


        // Motion

            //- Set the rigid-body state
            void setState
            (
                const point& origin,
                const tensor& R,
                const vector& velocity,
                const vector& omega
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
