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
    Foam::heightFieldBed

Description
    Single-valued sediment bed z_b(x, y) for the immersed boundary, on a
    finite-area mesh (faMesh) built at run time from a flat floor patch
    below the sediment. The sediment physics (bed load, slope effects,
    Exner equation, sand slide, fixed strips, volume balance, bed surface
    output) is in bedSediment; this class couples it to the immersed
    boundary:

    - Distance of a cell: (z_c - z_b)/sqrt(1 + |grad z_b|^2), normal
      (-grad z_b, 1)/sqrt(1 + |grad z_b|^2), with z_b of the floor face
      below the cell centre (column map built once; values exchanged with
      a mapDistribute). Gravity is along -z.
    - Bed shear stress per floor face: area-weighted average of the
      wall shear stress of the IB cells in its column.
    - Floor faces without bed IB cells in their column (covered by a
      body-fitted structure or an immersed solid: IB cells whose nearest
      surface is a solid do not count) are inactive: no sediment flux
      across their edges. They become active when erosion exposes them.
    - The height field counts no sediment inside an immersed solid: the
      vertical spans of the solids in each column are traced once and
      passed to the sediment core, which moves z_b through a span at no
      sediment cost and stops deposition at a span's underside.
    - The cells are re-classified when z_b has changed by more than
      updateDistance since the last classification.

Usage
    \verbatim
    bed
    {
        type            heightField;
        floorPatch      bottom;         // flat patch below the sediment
        initialElevation 0;             // z_b at t = 0 if there is no zb
                                        // field file: a constant, or any
                                        // PatchFunction1 on floorPatch:
                                        // { type expression; expression
                                        // #{ 0.01*pos().x() #}; } or
                                        // { type mappedFile; mapMethod
                                        // planar; } for survey points in
                                        // constant/boundaryData/<floorPatch>
        initialSurface  "bed.stl";      // optional, instead: z_b from an
                                        // STL/OBJ in constant/triSurface
                                        // (vertical line through each face)
        rigidElevation  0;              // optional non-erodible floor, as
        rigidSurface    "rock.stl";     // initialElevation/initialSurface
        openBoundaries  (inlet outlet); // edges with zero-gradient flux

        grainDiameter   0.00026;        // d50 [m]
        specificGravity 2.65;
        porosity        0.4;
        gravity         9.81;           // [m/s2]
        criticalShields 0.05;           // flat-bed theta_c0
        frictionCoefficient 0.63;       // mu_s = tan(angle of repose)
        bedLoad         engelundFredsoe;// | meyerPeterMuller

        morphFactor     1;
        morphStartTime  0;              // bed frozen before this time [s]
        maxBedChange    1e-4;           // per Exner sub-step [m]
        updateDistance  5e-4;           // re-classify after this change [m]

        // Testing only: either replaces the IB wall shear stress
        // by a PatchFunction1<vector> on the floor patch ...
        // prescribedShear { type expression; expression #{ vector(1e-3, 0, 0) #}; }
        // ... or the Exner advance by a prescribed elevation z_b(x, y, t)
        // prescribedElevation { type expression; expression #{ -0.02*time() #}; }
        // ... and the soil region's liquefied depth [m] (liquefied bed)
        // prescribedLiquefiedDepth { type uniformValue; uniformValue 0.01; }

        fixedBed                        // optional
        {
            inlet       0.05;           // fixed strip width [m] per patch
        }
        fixedBedLoad    capacity;       // bed load on the fixed strips:
                                        // capacity (equilibrium feed) |
                                        // none (rigid apron, no supply);
                                        // they never entrain suspended
                                        // sediment

        bedSlope                        // optional
        {
            model       kochFlokstra;   // | talmon | none
            alpha       1.5;            // kochFlokstra
            flowDepth   0.2;            // talmon: reference depth h [m]
        }

        sandSlide       true;
        slideDiffusivity 1e-3;          // K0 [m2/s] (pseudo-time units)
        slideTolerance  0.02;           // relative excess of tan(phi_r)
        slideMaxIter    500;
        liquefiedSlide                  // optional: seepage-weakened
        {                               // faces lower at most at this
            velocity    0.05;           // rate [m/s] (liquefied: also
        }                               // at most the liquefied depth)
        seepage                         // optional (seepageModel):
        {                               // chengChiew | pressureGradient |
            type        chengChiew;     // francalanci | none
        }
        // Two-phase solver (ibScourInterFoam): still-water level and air
        // pressure that the soil and the solids see as zero load
        referenceLevel  0.3;            // default: the solver's hRef
        atmosphericPressure 0;          // default 0 (the pRefValue)
        reportSolidForces (cylinder);   // optional: force and moment of
                                        // the flow on these solids every
                                        // step, postProcessing/solidForces
        soil                            // optional (soilRegion)
        {
            nLayers     40;
            bottom      -2;
            grading     20;
            type        darcyStorage;
            conductivity 1e-4;
            saturation  0.99;
        }
        structure                       // optional, with a soil region
        {                               // one rigid body
            solids          (pile);     // its immersed solids (default:
                                        // all solids except the pipes)
            hydrodynamicForce yes;      // flow force and moment (default)
            density         2400;       // and the rigidStructure entries
            freeTranslations (0 0 1);   // (free dofs, stiffness
            freeRotations   (0 1 0);    // estimates, external loads,
            ...                         // referencePoint, volume)
        }
        liquefiedPipes                  // optional, with a soil region
        {
            pipe1
            {
                solid       pipe;       // immersed solid moved vertically
                x           0.45;       // column of the pipe [m]
                specificGravity 3.1;    // and the other liquefiedPipe
                ...                     // entries (diameter and initial
            }                           // burial from the solid geometry)
        }
    }
    \endverbatim

    Liquefied pipes: each moves the named immersed solid
    vertically with the velocity of a liquefiedPipe driven by the soil
    region at its column (sinking or floatation in a wave-liquefied
    seabed); its diameter and initial burial are measured from the solid
    (signed-distance probes above and below the column). With
    'hydrodynamicForce' (default yes) the vertical force of the flow on the
    solid (immersedBoundary::solidForce) per unit length ('length',
    default the y extent of the mesh) is added to its net weight.

See also
    Foam::bedSediment

SourceFiles
    heightFieldBed.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_heightFieldBed_H
#define Foam_heightFieldBed_H

#include "bedModel.H"
#include "bedSediment.H"
#include "soilRegion.H"
#include "liquefiedPipe.H"
#include "rigidStructure.H"
#include "mapDistribute.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class heightFieldBed Declaration
\*---------------------------------------------------------------------------*/

class heightFieldBed
:
    public bedModel
{
    // Private Data

        //- Sediment core on the floor patch
        bedSediment sediment_;

        //- Column (compact floor-face index) of every cell
        labelList cellColumn_;

        //- Map fetching floor-face data into compact column addressing
        autoPtr<mapDistribute> columnMap_;

        //- Bed elevation at the last classification [m]
        scalarField zbClassified_;

        //- Re-classify after this change of z_b [m]
        scalar updateDistance_;

        //- Floor patch (the finite-area faces are its faces, in order)
        const polyPatch& floorPatch_;

        //- Optional prescribed bed shear stress [m2/s2] on the floor patch
        //- (testing: replaces the IB wall shear; the active flags still
        //- come from the classification)
        autoPtr<PatchFunction1<vector>> prescribedShear_;

        //- Optional prescribed bed elevation [m] on the floor patch
        //- (testing: replaces the Exner advance)
        autoPtr<PatchFunction1<scalar>> prescribedElevation_;

        //- Optional prescribed liquefied-layer depth [m] on the floor
        //- patch (testing of the liquefied bed without a soil region;
        //- replaces the soil region's value)
        autoPtr<PatchFunction1<scalar>> prescribedLiquefiedDepth_;

        //- Optional soil region under the bed ('soil' sub-dictionary):
        //- pore pressure driven by the bed pressure, seepage gradient fed
        //- back to the sediment core
        autoPtr<soilRegion> soil_;

        //- Water density [kg/m3] converting the kinematic flow pressure
        //- to the soil pressure (from the soil dictionary)
        scalar rhoW_;

        //- Two-phase solver (stages 2 and 3 of the free-surface plan): the
        //- soil takes the excess of the total bed pressure over the
        //- still-water hydrostatic pressure, p_b - p_atm - rho_w g
        //- (z_ref - z_b), and the solids the excess of the IB-cell
        //- pressure over the same, with z_ref the still-water level (entry
        //- referenceLevel of the bed or soil dictionary, default the
        //- solver's hRef) and p_atm the pressure of the air (entry
        //- atmosphericPressure, default 0, the pRefValue of the solver)
        scalar fsReferenceLevel_;
        bool fsReferenceLevelSet_;
        scalar fsAtmosphericPressure_;

        //- Solids whose hydrodynamic force and moment (about the origin)
        //- are reported every step (entry reportSolidForces) and written
        //- to postProcessing/solidForces/<name>.dat
        wordList reportSolids_;

        //- Liquefied pipes: settings, models (built at the first evolve,
        //- from the solid geometry) and the names of their solids
        dictionary pipesDict_;
        PtrList<liquefiedPipe> pipes_;
        wordList pipeSolids_;

        //- Probe line of each pipe: x, y of its column
        List<vector2D> pipeLines_;

        //- Length of each pipe in the domain [m] for the hydrodynamic
        //- force per unit length (0: no hydrodynamic force)
        scalarList pipeLengths_;

        //- Rigid structure of the 'structure' sub-dictionary:
        //- the body formed by its immersed solids, masked in the soil
        //- region with the solids' initial geometry, moved in the flow
        //- by the soil balance, loaded by the hydrodynamic force and
        //- moment of the flow
        autoPtr<rigidStructure> structure_;
        wordList structureSolids_;
        List<point> structureOrigins_;
        List<tensor> structureRs_;
        bool structureHydroForce_;
        bool structureRefSet_;

        //- Structure state at the previous evolve (for the velocities)
        vector structureDPrev_;
        vector structureThetaPrev_;
        vector structureV_;
        vector structureOmega_;


    // Private Member Functions

        //- Build the column map (cells -> floor faces)
        void makeColumnMap();

        //- Values of a floor-patch function on the finite-area faces
        template<class Type>
        tmp<Field<Type>> floorValues
        (
            const PatchFunction1<Type>& func
        ) const;

        //- Floor-face values in compact column addressing
        template<class Type>
        List<Type> columnValues(const UList<Type>& faceValues) const;

        //- Bed shear stress and active flags per floor face from the IB
        //- cells
        void collectShearStress(const immersedBoundary& ib);

        //- Bed pressure [Pa] per floor face: area-weighted average of the
        //- kinematic pressure of the IB cells of its column times rhoW
        //- (two-phase solver: the excess of the total bed pressure over
        //- the still-water hydrostatic pressure)
        tmp<scalarField> collectPressure(const immersedBoundary& ib) const;

        //- Hydrodynamic pressure [Pa] of the IB cells for the solid
        //- forces: rhoW p in the single-phase solver; in the two-phase
        //- solver the total pressure at the wall point of the cell
        //- (cell value extrapolated with the hydrostatic gradient of its
        //- phase) less the still-water hydrostatic pressure there,
        //- p_w - p_atm - rho_w g max(z_ref - z_w, 0): zero in water and
        //- air at rest, rho_w g eta under a wave, minus the lost
        //- hydrostatic pressure on a drained wall
        tmp<scalarField> ibPressure(const immersedBoundary& ib) const;

        //- Density [kg/m3] of the IB cells for the wall shear stress:
        //- rhoW, or the mixture density in the two-phase solver
        tmp<scalarField> ibDensity(const immersedBoundary& ib) const;

        //- Vertical spans of the immersed solids in every floor-face
        //- column (traced along the column line) for the sediment core;
        //- recomputed when the solids change
        void updateSolidSpans(const immersedBoundary& ib);

        //- Structure mask of the soil cells from the immersed solids
        //- (cell centre inside a solid), after a reclassification
        void updateStructureMask(const immersedBoundary& ib);

        //- The soil columns moved since the last mask update
        bool maskDirty_;

        //- True until the solid spans have been traced
        bool spansDirty_;

        //- Bed elevation at the local points of the floor patch, from the
        //- point level set of the soil surface [m]
        tmp<scalarField> zbAtFloorPoints() const;

        //- Surface of the soil region [m]: z_b, except in the columns
        //- covered by an immersed solid (z_b inside a solid span: the
        //- solid replaces the sand there and z_b is the solid's bottom),
        //- which take the level of the sand around the solid (average of
        //- the uncovered neighbours, outward passes), so that the soil
        //- surface along the solid's perimeter does not dip to the
        //- solid's bottom and the soil cells inside the solid above it
        //- exist for the structure mask (lateral contact down from the
        //- sand surface)
        tmp<areaScalarField> soilSurface() const;

        //- Build the liquefied pipes from the geometry of their solids
        void buildPipes(const immersedBoundary& ib);

        //- Bed elevation and lowest and highest point of the solid on the
        //- probe line (x, y) at column x [m]
        void probePipe
        (
            const immersedBoundary& ib,
            const word& solidName,
            const scalar x,
            const scalar L,
            scalar& zBed,
            scalar& y,
            scalar& zLow,
            scalar& zHigh
        ) const;


public:

    //- Runtime type information
    TypeName("heightField");


    // Constructors

        //- Construct from mesh and dictionary
        heightFieldBed(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~heightFieldBed() = default;


    // Member Functions

        //- Signed distance and normal of all cell centres
        //- Level set phi = z - z_b at the mesh points: average over the
        //- cells of each point of the column plane z_b + grad(z_b).(x - x_f)
        //- of the cell, synchronised across processors
        virtual tmp<scalarField> pointLevel() const;

        //- The same for a given surface field
        tmp<scalarField> pointLevel(const areaScalarField& zb) const;

        virtual void cellDistance
        (
            scalarField& distance,
            vectorField& normal
        ) const;

        //- Advance the bed over the flow time step
        virtual bool evolve(const immersedBoundary& ib, const scalar dt);

        //- Record the bed at the latest classification
        virtual void classified();

        //- The solids moved: trace the solid spans of the columns again
        virtual void solidsMoved()
        {
            spansDirty_ = true;
        }

        virtual void setSuspendedExchange
        (
            const immersedBoundary& ib,
            const scalarField& netDeposition
        );

        //- Erodibility per IB cell: 0 in the columns of the fixed strips
        virtual tmp<scalarField> erodibleIbCells
        (
            const immersedBoundary& ib
        ) const;

        //- Velocities of the solids moved by liquefied pipes and by the
        //- structure
        virtual void solidVelocities(HashTable<vector>& velocities) const;

        //- Angular velocities of the solids of the structure
        virtual void solidAngularVelocities(HashTable<vector>& omegas) const;

        //- The rigid structure (nullptr: none)
        const rigidStructure* structure() const noexcept
        {
            return structure_.get();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
