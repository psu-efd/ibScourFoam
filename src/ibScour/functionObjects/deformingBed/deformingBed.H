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
    Foam::functionObjects::deformingBed

Group
    grpFieldFunctionObjects

Description
    Erodible bed of a body-fitted flow mesh that deforms with the bed
    (mesh deformation approach of Liu & Garcia 2008, J. Waterway Port
    Coastal Ocean Eng. 134:203-217, and Roulund et al. 2005, J. Fluid Mech.
    534:351-401).

    The bed is a patch of the flow mesh. Its plan view is copied once into
    a static bed mesh (bedRegionMesh), on which the sediment core
    (bedSediment) solves the bed load, slope effects, the Exner equation
    and the sand slide for z_b(x, y), exactly as for the immersed-boundary
    bed (heightFieldBed). After every flow time step:

    - the kinematic wall shear stress of each bed face,
      tau = -nu_eff (I - n n) & snGrad(U), with nu_eff = nu + nu_t of the
      wall boundary condition (the stress the flow actually transmits),
      projected on the horizontal plane;
    - z_b is advanced over the step (times morphFactor);
    - the bed points are moved vertically to the bed elevation (the
      inverse-distance average of the surrounding faces): the fixedValue
      condition of the bed patch of pointDisplacementz (vertical component
      solver displacementComponentLaplacian, component z; recommended for
      meshes with vertical columns, with 'diffusivity directional
      (1e-3 1e-3 1)' each column stretches linearly between bed and lid)
      or the vertical component of pointDisplacement (displacementLaplacian)
      is set, and the motion solver of dynamicMotionSolverFvMesh moves the
      interior points at the next mesh update.

    Without a z_b field (a new run) the initial bathymetry is that of the
    mesh: z_b is the elevation of the bed-face centres. The mesh follows
    z_b from the first step on, also before morphStartTime (the bed
    points are set to the inverse-distance average of the face values).
    The bed moves vertically only, so the static plan view stays valid.
    Use movingWallVelocity for U on the bed patch so that the flow sees
    the bed velocity.

Usage
    In system/controlDict:
    \verbatim
    functions
    {
        bed
        {
            type            deformingBed;
            libs            (ibScour);
            bedPatch        bed;          // patch of the flow mesh
            C               C;              // optional: C field whose bed
                                            // condition gives E and D
            U               U;            // optional
            moveMesh        true;         // optional

            // Sediment entries as for the heightField bed model
            // (grainDiameter, criticalShields, bedLoad, morphFactor,
            // morphStartTime, openBoundaries, bedSlope, fixedBed,
            // sandSlide, ...), see bedSediment
            grainDiameter   0.00026;
            openBoundaries  (inlet outlet);
        }
    }
    \endverbatim

    The finite-area schemes of the bed mesh are read from
    system/bedRegion/finite-area/{faSchemes,faSolution}; z_b is written to
    <time>/bedRegion/finite-area/zb (processor directories in parallel)
    and read back on restart.

See also
    Foam::bedSediment
    Foam::bedRegionMesh
    Foam::heightFieldBed

SourceFiles
    deformingBed.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_deformingBed_H
#define Foam_functionObjects_deformingBed_H

#include "fvMeshFunctionObject.H"
#include "bedSediment.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class deformingBed Declaration
\*---------------------------------------------------------------------------*/

class deformingBed
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Bed patch of the flow mesh
        word patchName_;

        //- Index of the bed patch
        label patchID_;

        //- Name of the velocity field
        word UName_;

        //- Name of the suspended-sediment concentration field whose bed
        //- patch condition (sedimentBedFlux) supplies E and D; empty or
        //- absent: no suspended load
        word CName_;

        //- Move the mesh with the bed
        bool moveMesh_;

        //- Kinematic viscosity [m2/s]
        scalar nu_;

        //- Static plan-view mesh of the bed
        autoPtr<polyMesh> bedMeshPtr_;

        //- Sediment core on the bed mesh
        autoPtr<bedSediment> sedimentPtr_;


    // Private Member Functions

        //- Wall shear stress of the bed faces into the sediment core
        void calcShearStress();

        //- Set the bed-patch point displacement to the bed elevation
        void setBedDisplacement();


public:

    //- Runtime type information
    TypeName("deformingBed");


    // Constructors

        //- Construct from name, Time and dictionary
        deformingBed
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        deformingBed(const deformingBed&) = delete;

        //- No copy assignment
        void operator=(const deformingBed&) = delete;


    //- Destructor
    virtual ~deformingBed() = default;


    // Member Functions

        //- Read the settings
        virtual bool read(const dictionary& dict);

        //- Advance the bed over the last flow time step and move the bed
        //- points
        virtual bool execute();

        //- No separate output (the bed surface is written by execute)
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
