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
    Foam::bedModel

Description
    Abstract base class of an evolving sediment bed for the immersed
    boundary.

    A bed model provides the signed distance (positive in the fluid) and
    the unit normal (into the fluid) of every cell of the mesh, and evolves
    the bed from the wall shear stress of the IB cells (Exner equation).
    It is held by immersedBoundary next to the immersed solids; the cell
    classification uses the minimum distance over solids and bed.

Usage
    In constant/immersedBoundaryDict:
    \verbatim
    bed
    {
        type            heightField;
        // type-specific entries
    }
    \endverbatim

SourceFiles
    bedModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_bedModel_H
#define Foam_bedModel_H

#include "fvMesh.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "HashTable.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class immersedBoundary;

/*---------------------------------------------------------------------------*\
                          Class bedModel Declaration
\*---------------------------------------------------------------------------*/

class bedModel
{
protected:

    // Protected Data

        //- The mesh
        const fvMesh& mesh_;

        //- Settings
        dictionary dict_;


public:

    //- Runtime type information
    TypeName("bedModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            bedModel,
            dictionary,
            (const fvMesh& mesh, const dictionary& dict),
            (mesh, dict)
        );


    // Constructors

        //- Construct from mesh and dictionary
        bedModel(const fvMesh& mesh, const dictionary& dict);


    // Selectors

        //- Return the selected bed model
        static autoPtr<bedModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~bedModel() = default;


    // Member Functions

        //- Signed distance (> 0 in fluid) and unit normal (into the fluid)
        //- of all cell centres
        virtual void cellDistance
        (
            scalarField& distance,
            vectorField& normal
        ) const = 0;

        //- Level set of the bed at the mesh points, phi = z - z_b(x, y)
        //- (> 0 in the fluid), consistent across processors; used for the
        //- face apertures and volume fractions of the cut-face immersed
        //- boundary. Not implemented by default.
        virtual tmp<scalarField> pointLevel() const;

        //- Advance the bed over the flow time step dt [s] using the wall
        //- shear stress of the IB cells. Returns true when the bed has moved
        //- enough since the last classification to require a new one.
        virtual bool evolve
        (
            const immersedBoundary& ib,
            const scalar dt
        ) = 0;

        //- Net deposition rate of suspended sediment per IB cell [m/s]
        //- (grain volume per wall area and time, D - E, in the order of
        //- ib.ibCells()) for the Exner equation of the next advance;
        //- default: ignored
        virtual void setSuspendedExchange
        (
            const immersedBoundary& ib,
            const scalarField& netDeposition
        )
        {}

        //- Erodibility (1 erodible, 0 fixed) per IB cell in the order of
        //- ib.ibCells(), for the suspended-load pickup: a fixed bed (a
        //- rigid apron) entrains no sediment; default: invalid tmp, all
        //- cells erodible
        virtual tmp<scalarField> erodibleIbCells
        (
            const immersedBoundary& ib
        ) const
        {
            return tmp<scalarField>();
        }

        //- Called after the immersed boundary has re-classified the cells
        virtual void classified() = 0;

        //- Called when an immersed solid has moved (prescribed motion or a
        //- structure); the bed refreshes what it caches about the solids
        virtual void solidsMoved()
        {}

        //- Velocities [m/s] the bed prescribes to immersed solids, by
        //- solid name, after the last evolve (e.g. a pipe sinking in a
        //- liquefied seabed); default: none
        virtual void solidVelocities(HashTable<vector>& velocities) const
        {}

        //- Angular velocities [rad/s] the bed prescribes to immersed
        //- solids, by solid name (structures that rotate); default: none
        virtual void solidAngularVelocities(HashTable<vector>& omegas) const
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
