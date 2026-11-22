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
    Foam::columnStretchMotionSolver

Description
    Algebraic mesh motion for a quasi-two-dimensional body-fitted bed next
    to a rigid apron with a vertical downstream face (the wall-jet layout
    of Liu and Garcia 2008): every column of points above the erodible bed
    is stretched vertically between the bed and the top of the domain, the
    points on the vertical face slide down it as the scour hole deepens,
    and the apron side of the face does not move.

    Points are moved in z only, as a function of their original x and z.
    Upstream of the face (x < xFace) the displacement is zero. Downstream
    the full column stretch of a point at height z0 above a bed that moved
    by dz_b is dz_b (H - z0)/(H - z_b0), with H the top of the domain and
    z_b0 the original bed level. Points above the apron level take the
    fraction w = (x - xFace)/blendLength of it (capped at one), so that the
    shear between the fixed column of the face and the moving columns is
    spread over blendLength instead of one cell; points below the apron
    level (the gap block between the bed and the apron level) are placed
    linearly between the moved bed and the moved apron-level point of
    their column, which at the face itself is the fixed apron corner, so
    they slide along the vertical face. The gap block keeps at least
    minGapFraction of its original height where sand rises against the
    face.

    The bed displacement is read from the fixedValue pointDisplacement
    condition of the bed patch, which the deformingBed function object
    maintains. The bed table is gathered from all processors.

Usage
    \verbatim
    dynamicFvMesh   dynamicMotionSolverFvMesh;
    motionSolverLibs (ibScour);
    motionSolver    columnStretch;
    columnStretchCoeffs
    {
        bedPatch        bed;
        xFace           0.66;   // vertical face of the apron
        apronLevel      0;      // z of the apron surface
        blendLength     0.1;    // shear spread above the apron level [m]
        minGapFraction  0.25;   // optional
    }
    \endverbatim

SourceFiles
    columnStretchMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_columnStretchMotionSolver_H
#define Foam_columnStretchMotionSolver_H

#include "displacementMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class columnStretchMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class columnStretchMotionSolver
:
    public displacementMotionSolver
{
    // Private Data

        //- Bed patch
        word bedPatchName_;
        label bedPatchID_;

        //- x of the vertical face of the apron
        scalar xFace_;

        //- z of the apron surface
        scalar apronLevel_;

        //- Horizontal length over which the stretch above the apron level
        //- grows from zero at the face to the full column stretch
        scalar blendLength_;

        //- Minimum height of the gap block as a fraction of its original
        scalar minGapFraction_;

        //- Top of the domain (largest original z)
        scalar zTop_;

        //- Tolerance in x for the bed table and the face
        scalar xTol_;


    // Private Member Functions

        //- Gather the bed table (x, original z, displacement z) sorted by x
        void bedTable
        (
            scalarField& x,
            scalarField& zb0,
            scalarField& dz
        ) const;

        //- Linear interpolation in a sorted table, clamped at the ends
        static scalar interpolate
        (
            const scalar x,
            const scalarField& xs,
            const scalarField& ys
        );

        //- No copy construct
        columnStretchMotionSolver(const columnStretchMotionSolver&) = delete;

        //- No copy assignment
        void operator=(const columnStretchMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("columnStretch");


    // Constructors

        //- Construct from mesh and dictionary
        columnStretchMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict
        );

        //- Construct from components
        columnStretchMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    //- Destructor
    virtual ~columnStretchMotionSolver() = default;


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
