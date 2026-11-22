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
    Foam::bedSediment

Description
    Sediment core of a single-valued bed z_b(x, y), independent of how the
    bed is coupled to the flow (immersed boundary or body-fitted mesh).

    The bed lives on a finite-area mesh (faMesh) built at run time from a
    flat, static patch whose plan view is the bed: a floor patch of the
    flow mesh (immersed boundary) or the bed patch of a separate bed mesh
    (body-fitted, deforming flow mesh; see bedRegionMesh). z_b is a field
    on it, so all geometry is that of the horizontal projection. The coupling supplies the
    bed shear stress and the active flags per face; the core computes the
    bed load and advances z_b:

    - Exner equation in flux form, (1 - n) dz_b/dt + div(q) = 0, explicit
      with upwind edge fluxes and sub-cycling limited by maxBedChange;
      morphFactor multiplies the flow time step.
    - Adaptive morphological factor ('bedCourant' entry, optional, 0 =
      off): the explicit flow-bed coupling (shear frozen over the
      morphological step of one flow step) develops a growing even-odd
      oscillation of the bed once the bed change per flow step is a
      sizeable fraction of the bed-face size (pier runs at a fixed factor
      failed at MF dt above 0.015-0.02 s). With bedCourant C > 0 the
      morphological step of the flow step is limited so that
          max over faces |dz_b| / sqrt(S) <= C
      from the bed-change rate of the first Exner sub-step; the factor in
      force is min(morphFactor, that limit), and may grow by at most
      'morphFactorGrowth' (default 1.2) per flow step. morphFactor is then
      the maximum factor. The morphological time actually advanced,
      t_m, is kept in <time>/uniform/bedMorphology (restart-safe), printed
      with the bed report ("Bed: MF <in force>, t_m <s>"), and appended to
      postProcessing/bed/morphTime.dat at every bed write.
    - Bed load: engelundFredsoe (Roulund et al. 2005 slope-corrected
      critical Shields number, q* = 18.74(theta - theta_c)(sqrt(theta)
      - 0.7 sqrt(theta_c))) or meyerPeterMuller (q* = 8(theta -
      theta_c)^1.5); direction of the bed shear stress.
    - Gravitational (down-slope) bed load ('bedSlope' entry, optional):
      q = q0 (t - gamma grad z_b) with gamma = 1/f(theta), as an edge flux
      -gamma |q0| dz_b/dn with a compact two-point gradient (which also
      damps grid-scale bed modes):
        - kochFlokstra: gamma = alpha sqrt(theta_c0/theta) (Koch &
          Flokstra 1981; Struiksma et al. 1985, J. Hydraul. Res.
          23:57-79), alpha 1.5 by default;
        - talmon: gamma = 1/(9 (d50/h)^0.3 sqrt(theta)) (Talmon et al.
          1995, J. Hydraul. Res. 33:495-517) with a reference flow
          depth h;
        - none (default).
      The formulations were calibrated for transverse slopes; here the
      term acts isotropically, in addition to the slope correction of
      the critical Shields number. The coefficient is poorly constrained
      (Baar et al. 2019, Nat. Commun. 10:4903). The explicit diffusion
      adds a positivity limit to the Exner sub-step.
    - Bed-load adaptation (saturation) length L ('adaptationLength',
      optional, 0 = local equilibrium by default): the magnitude of the
      bed load relaxes to its equilibrium value along the transport
      direction e (direction of the bed shear stress),
          L (e . grad) q = q_eq - q,
      the steady non-equilibrium bed-load model (adaptation length:
      Phillips & Sutherland 1989, J. Hydraul. Res. 27:105-127; Wu 2007,
      Computational River Dynamics; saturation length: Charru 2006, Phys.
      Fluids 18:121508; Fourriere et al. 2010, J. Fluid Mech. 649:287-328).
      With a grain-scale length (10-20 d50) it restores the short-wave
      cut-off of the bedform instability, which otherwise grows fastest at
      the grid scale. Discretised with upwind inflow edges and solved by
      Jacobi sweeps; sediment enters from open boundaries at equilibrium
      and not from closed ones (walls, structures).
    - Sand slide: conservative, gradient-limited diffusive flux on edges
      steeper than the angle of repose (Song, Xu & Liu 2020), iterated to
      tolerance with an iteration cap.
    - Inactive faces (flag set by the coupling, e.g. covered by a
      structure) exchange no sediment.
    - Optional fixed (non-erodible) strips next to open boundaries
      (fixedBed): faces whose centre lies within the given horizontal
      distance of the patch keep their elevation. Bed load still crosses
      them (a rigid apron); their exchange is booked as boundary flux in
      the volume balance.
    - Initial bed: the area field zb, else 'initialSurface' (an STL/OBJ
      surface sampled along vertical lines) or 'initialElevation' (any
      PatchFunction1 on the bed patch: a constant, an expression, or
      mappedFile for scattered survey points with 2-D Delaunay
      interpolation).
    - The finite-area mesh is read from constant/finite-area when it exists
      (makeFaMesh with system/finite-area/faMeshDefinition; decomposePar and
      reconstructPar then handle the bed fields), else built from the bed
      patch at run time.
    - Optional non-erodible floor z_r below the sediment ('rigidElevation'
      or 'rigidSurface' as for the initial bed, or the area field zRigid, e.g. a rigid flume
      bed, a sill or a scour protection): a face gives at most the sand
      it holds above the floor, (z_b - z_r)(1 - n) A per sub-step. Each
      Exner sub-step takes two passes: the first finds the outgoing flux
      of every face, the second scales the edge fluxes of a face whose
      outgoing volume exceeds its sand (a non-erodible layer as in
      Delft3D and TELEMAC, applied per donor so the update stays
      conservative). The sand slide skips empty donors. A
      residual dip below z_r (suspended pickup from a bare floor) is
      clamped and its volume reported.
    - Global sediment-volume balance reported every step.
    - Bed surface written to postProcessing/bed/<time>/ (VTK, merged in
      parallel) with z_b, the bed shear stress, the bed-load flux, the
      active and fixed flags.

    The dictionary entries are those of the bed model that owns the core
    (see heightFieldBed).

SourceFiles
    bedSediment.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_bedSediment_H
#define Foam_bedSediment_H

#include "fvMesh.H"
#include "faMesh.H"
#include "Pair.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "Function1.H"
#include "seepageModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class bedSediment Declaration
\*---------------------------------------------------------------------------*/

class bedSediment
{
    // Private Data

        //- The flow mesh (time, open-boundary patches for fixedBed)
        const fvMesh& mesh_;

        //- The mesh carrying the faMesh patch (the flow mesh or a static
        //- bed mesh)
        const polyMesh& bedMesh_;

        //- Settings
        dictionary dict_;

        //- Finite-area mesh of the bed patch
        autoPtr<faMesh> aMeshPtr_;

        //- Bed elevation [m]
        autoPtr<areaScalarField> zbPtr_;

        //- Bed shear stress (kinematic) [m2/s2]
        autoPtr<areaVectorField> tauPtr_;

        //- Bed-load flux [m2/s]
        autoPtr<areaVectorField> qPtr_;

        //- Active bed faces (1); inactive faces (0) exchange no sediment
        autoPtr<areaScalarField> activePtr_;

        //- Fixed (non-erodible) bed faces (1), e.g. a strip at the inlet
        autoPtr<areaScalarField> fixedPtr_;

        //- Vertical spans (z_bottom, z_top) of the immersed solids in the
        //- column of each face, bottom-up; the height field holds no
        //- sediment inside them (set by the coupling, empty without solids)
        List<List<Pair<scalar>>> solidSpans_;

        //- Faces at or inside a solid span (1): full up to the underside
        //- of the solid, they receive no sediment (zero without solids)
        autoPtr<areaScalarField> cappedPtr_;

        //- Any processor has solid spans (the guard of every branch that
        //- holds a collective; a local span count differs between ranks,
        //- e.g. zero on a rank without floor faces)
        bool haveSolids_ = false;

        //- Elevation of the non-erodible floor below the sediment [m]
        //- ('rigidElevation', 'rigidSurface' or the field file zRigid;
        //- absent: none)
        autoPtr<areaScalarField> zRigidPtr_;

        //- A non-erodible floor is set
        bool rigidFloor_;

        //- Cumulative volume [m3] added by clamping z_b to the rigid floor
        //- (suspended pickup from a bare floor; zero for bed load only)
        scalar cumulativeRigidClamp_;

        //- Fixed-strip width [m] per boundary patch
        dictionary fixedBed_;

        //- Open boundary patches of the faMesh
        wordList openBoundaries_;

        //- Index of the patch of bedMesh_ carrying the faMesh
        label patchID_;


        // Sediment and model parameters

            scalar d50_;
            scalar s_;
            scalar porosity_;
            scalar g_;
            scalar thetaC0_;
            scalar mus_;

            //- Bed-load formula: engelundFredsoe, meyerPeterMuller or
            //- vanRijn (van Rijn 1984, q* = 0.053 T^2.1/D*^0.3)
            word bedLoad_;

            //- Kinematic viscosity of the water [m2/s] for the van Rijn
            //- particle parameter D* (default 1e-6)
            scalar nu_;

            //- Morphological acceleration factor as a function of time
            //- (constant, table, ramp, ...; default 1)
            autoPtr<Function1<scalar>> morphFactor_;

            //- Passes of the three-point (edge-neighbour) smoothing of the bed
            //- shear stress before the bed load is computed (default 1;
            //- 0: none). Removes
            //- the two-cell alternation of the shear on inclined cut-cell beds
            //- while conserving its mean (gain cos^2(pi dx/L) per pass)
            label smoothShear_;
            scalar morphStartTime_;
            scalar maxBedChange_;

            //- Bed Courant number of the adaptive morphological factor
            //- (max |dz_b| per flow step over the face size; 0 = fixed
            //- factor)
            scalar bedCourant_;

            //- Largest bed change per flow step [m] (0 = no cap): limits
            //- the factor in force like bedCourant, with an absolute
            //- length instead of a fraction of the face size; on a
            //- deforming mesh it caps the mesh velocity at the bed
            scalar maxStepChange_;

            //- Growth limit of the factor in force per flow step
            scalar morphFactorGrowth_;

            //- Factor in force at the last advance
            scalar morphFactorInForce_;

            //- Morphological time advanced so far [s] (restart-safe)
            scalar morphTime_;

            //- Flow time of the last advance [s]
            scalar lastAdvanceTime_;
            word slopeModel_;
            scalar slopeAlpha_;
            scalar slopeDepth_;
            bool sandSlide_;
            scalar slideK0_;
            scalar slideTolerance_;
            label slideMaxIter_;

            //- Rate of the seepage-induced flank failure: a face with
            //- seepage-reduced friction lowers at most at this velocity
            //- [m/s] per step, and when liquefied by its liquefied depth
            //- (0: the slide re-grades the whole flank at once)
            scalar liquefiedSlideVelocity_;

            //- Liquefied bed ('liquefiedBed' sub-dictionary): on
            //- faces whose liquefied depth exceeds
            //- minDepth the critical Shields number and the repose
            //- friction are multiplied by criticalShieldsFactor and
            //- reposeFactor (1: no effect, the default), ramped linearly
            //- over transitionDepth so that the transport stays smooth
            scalar liqMinDepth_;
            scalar liqTransition_;
            scalar liqThetaCFactor_;
            scalar liqReposeFactor_;

            //- Bed-load adaptation length [m] (0: local equilibrium)
            scalar adaptLength_;

            //- Bed load on the fixed strips: true = the transport capacity
            //- of their shear (an equilibrium feed through an inlet strip),
            //- false = zero (a rigid apron without sand supply): the strip
            //- then donates no sediment across its edges either, so the
            //- slope correction cannot feed the scour hole next to it
            bool fixedBedLoad_;

            //- Relative tolerance and iteration cap of the adaptation solve
            scalar adaptTolerance_;
            label adaptMaxIter_;

            //- Averaging window of the bed-load flux [s] (0: the
            //- instantaneous flux of every flow step). Under waves the
            //- Exner equation takes the mean flux of the last completed
            //- window (a wave period), so the bed does not oscillate with
            //- the orbital motion at the accelerated rate
            scalar fluxAveragingPeriod_;

            //- Time accumulated in the current averaging window
            scalar fluxWindowTime_;

            //- A completed window exists (the bed is frozen until then)
            bool haveMeanFlux_;

            //- Largest instantaneous |q| seen in the current window and the
            //- face centre where it occurred (diagnostic: a window mean far
            //- above the instantaneous flux of the writes means spikes
            //- between them)
            scalar windowMaxQ_;
            point windowMaxQPoint_;

            //- Upper bound of the Shields number that enters the bed-load
            //- formula (maxShields, 0: none). A single-step garbage shear
            //- on one face (1e4 m2/s2 on the first ring of the wave basin,
            //- 2026-10-31) otherwise dominates the window mean; the faces
            //- clipped in a window are counted and reported
            scalar maxShields_;
            label windowClipped_;


        //- Iterations of the last adaptation solve
        label adaptIter_;

        //- Running sum of q dt over the current averaging window
        autoPtr<areaVectorField> qSumPtr_;

        //- Mean bed-load flux of the last completed window
        autoPtr<areaVectorField> qMeanPtr_;


        //- z_b was read from a field file (restart or initial bathymetry)
        bool zbRead_;

        //- Cumulative sediment-volume error [m3]
        scalar cumulativeMassError_;

        //- Cumulative volume of z_b passage through solid spans [m3] (the
        //- solid overlap the height field released as bed change, not
        //- sediment) and of deposition blocked under solids [m3]
        scalar cumulativeOverlap_;
        scalar cumulativeBlocked_;

        //- Down-slope transport coefficient gamma |q0| per face [m2/s]
        scalarField slopeCoeff_;

        //- Seepage effect on sediment mobility (none by default)
        autoPtr<seepageModel> seepage_;

        //- Upward hydraulic gradient at the bed faces, set by the coupling
        //- with the soil region (0 without one)
        autoPtr<areaScalarField> seepagePtr_;

        //- Depth of the momentarily liquefied layer under the bed faces [m],
        //- set by the coupling with the soil region (diagnostic)
        autoPtr<areaScalarField> liquefiedPtr_;

        //- Net deposition of suspended sediment D - E [m/s] (grain volume
        //- per bed area and time) of the current flow step (unset: none)
        autoPtr<areaScalarField> exchangePtr_;

        //- Capacity index of the soil column under the bed faces [Pa m]
        //- (diagnostic, set by the coupling with the soil region)
        autoPtr<areaScalarField> capacityPtr_;

        //- Apparent bed slope from the horizontal seepage force,
        //- grad_h p / (gamma' f) per face (zero without the
        //- pressureGradient model), set by the coupling
        autoPtr<areaVectorField> seepageSlopePtr_;


    // Private Member Functions

        //- Build the faMesh with one faPatch per neighbouring polyPatch
        void makeFaMesh(const word& patchName);

        //- Is the edge patch an open boundary (bed load in and out)?
        bool isOpenBoundary(const faPatch& p) const;

        //- Bed level [m] on the area faces from the bed dictionary: the
        //- STL/OBJ surface 'surfaceKey' in constant/triSurface (highest
        //- intersection of a vertical line through each face centre), else
        //- the PatchFunction1 'valueKey' on the bed patch (a constant, an
        //- expression, or mappedFile for scattered survey points), else
        //- defaultValue
        tmp<scalarField> bedLevelInput
        (
            const word& patchName,
            const word& valueKey,
            const word& surfaceKey,
            const scalar defaultValue
        ) const;

        //- Mark the faces of the fixed strips (fixedBed)
        void makeFixedMask();

        //- Smooth the bed shear stress over the active edge neighbours
        //- (smoothShear passes of tau <- tau/2 + mean(neighbours)/2)
        void smoothShearStress();

        //- Bed-load flux per face from the shear stress and slope
        void calcBedLoad();

        //- Relax the magnitude of the bed load over the adaptation length
        //- along the transport direction (in place)
        void adaptBedLoad();

        //- Edge fluxes of a face vector field (upwind), with open/closed
        //- boundaries; returns the net outflow per face [m3/s]. With
        //- donorScale, every edge flux is multiplied by the factor of its
        //- donor face (the non-erodible floor limiter); with outgoing, the
        //- sum of the fluxes leaving each face is returned in it
        scalarField netOutflow
        (
            const areaVectorField& q,
            scalar& boundaryOutflow,
            const areaScalarField* donorScale = nullptr,
            scalarField* outgoing = nullptr
        ) const;

        //- Sediment volume per unit area above the non-erodible floor
        //- [m] (GREAT without a rigid floor)
        tmp<scalarField> sedimentAboveRigid() const;

        //- Conservative sand slide; returns the number of iterations
        label sandSlide(const scalar dtMorph);

        //- Update the capped flags from z_b and the solid spans
        void updateCapped();

        //- Room below the underside of the next solid span above z_b of
        //- each face [m] (GREAT without one)
        tmp<scalarField> roomBelowSolid() const;

        //- Correct a change of z_b from zOld to z for the solid spans of
        //- the faces: erosion passes through a span at no sediment cost
        //- (the passage volume is added to overlap), deposition stops at
        //- the underside of a span (the rest is added to blocked); both
        //- volumes are sediment-grain volumes (1 - n) S dz, local sums
        void applySolidSpans
        (
            const scalarField& zOld,
            scalarField& z,
            scalar& overlap,
            scalar& blocked
        ) const;

        //- Down-slope coefficient gamma = 1/f(theta)
        scalar slopeGamma(const scalar theta) const;

        //- Multiplier of a bed property on a liquefied face: 1 off the
        //- liquefied faces, factor f where the liquefied depth exceeds
        //- minDepth + transitionDepth, linear in between
        scalar liquefiedFactor(const label facei, const scalar f) const;

        //- Down-slope coefficient as an area field (coupled neighbours)
        tmp<areaScalarField> slopeCoeffField() const;

        //- Largest stable Exner sub-step of the down-slope diffusion [s]
        scalar slopeTimeStep() const;


public:

    //- Runtime type information
    ClassName("bedSediment");


    // Constructors

        //- Construct on the given (flat, static) patch of bedMesh from
        //- the bed dictionary; mesh is the flow mesh
        bedSediment
        (
            const fvMesh& mesh,
            const polyMesh& bedMesh,
            const dictionary& dict,
            const word& patchName
        );

        //- No copy construct
        bedSediment(const bedSediment&) = delete;

        //- No copy assignment
        void operator=(const bedSediment&) = delete;


    //- Destructor
    ~bedSediment() = default;


    // Member Functions

        // Access

            //- The finite-area mesh of the bed
            const faMesh& aMesh() const
            {
                return *aMeshPtr_;
            }

            //- Index of the patch of the bed mesh carrying the faMesh
            label patchID() const noexcept
            {
                return patchID_;
            }

            //- The mesh carrying the faMesh patch
            const polyMesh& bedMesh() const noexcept
            {
                return bedMesh_;
            }

            //- Bed elevation [m]
            const areaScalarField& zb() const
            {
                return *zbPtr_;
            }

            //- True if z_b was read from a field file
            bool zbRead() const noexcept
            {
                return zbRead_;
            }

            //- Set z_b (e.g. the initial bathymetry from the mesh) [m]
            void setElevation(const scalarField& zb);

            //- Bed shear stress (kinematic) [m2/s2], set by the coupling
            areaVectorField& tau()
            {
                return *tauPtr_;
            }

            //- Fixed (non-erodible) flags (1/0) of the fixedBed strips
            const areaScalarField& fixed() const
            {
                return *fixedPtr_;
            }

            //- Active flags (1/0), set by the coupling
            areaScalarField& active()
            {
                return *activePtr_;
            }

            //- Spans [z_bottom, z_top] of the immersed solids in each
            //- column (empty without solids)
            const List<List<Pair<scalar>>>& solidSpans() const noexcept
            {
                return solidSpans_;
            }

            //- Upward hydraulic gradient at the bed faces (seepage)
            const areaScalarField& seepageGradient() const
            {
                return *seepagePtr_;
            }

            //- Set the upward hydraulic gradient at the bed faces from
            //- the soil region
            void setSeepageGradient(const scalarField& i);

            //- Set the liquefied-layer depth under the bed faces (diagnostic)
            void setLiquefiedDepth(const scalarField& depth);

            //- Set the net deposition of suspended sediment D - E [m/s]
            //- per face for the Exner equation of the next advance
            void setSuspendedExchange(const scalarField& netDeposition);

            //- Set the capacity index of the soil columns (diagnostic)
            void setCapacityIndex(const scalarField& cap);

            //- Set the vertical spans of the immersed solids per face
            //- (bottom-up pairs (z_bottom, z_top)); an empty outer list
            //- removes them
            void setSolidSpans(List<List<Pair<scalar>>>&& spans);

            //- Cumulative solid-overlap volume released as bed change [m3]
            scalar cumulativeOverlap() const noexcept
            {
                return cumulativeOverlap_;
            }

            //- Set the apparent slope of the horizontal seepage force,
            //- gradHp / (gamma' f) with gradHp the horizontal gradient of
            //- the pore pressure just below the bed [Pa/m] and gamma' the
            //- submerged unit weight [N/m3]; f is the effective-weight
            //- factor of the seepage model
            void setSeepageForce
            (
                const vectorField& gradHp,
                const scalar gammaS
            );

            //- The seepage model
            const seepageModel& seepage() const
            {
                return *seepage_;
            }

            //- Bed-load flux [m2/s]
            const areaVectorField& q() const
            {
                return *qPtr_;
            }

            //- Morphological acceleration factor at the current time
            scalar morphFactor() const
            {
                return
                (
                    morphFactor_
                  ? morphFactor_->value(mesh_.time().value())
                  : scalar(1)
                );
            }

            //- Factor in force at the last advance (the adaptive limit,
            //- or morphFactor())
            scalar morphFactorInForce() const
            {
                return morphFactorInForce_;
            }

            //- Morphological time advanced so far [s]
            scalar morphTime() const
            {
                return morphTime_;
            }

            //- True once the bed may move (time >= morphStartTime)
            bool morphActive() const
            {
                return mesh_.time().value() >= morphStartTime_;
            }


        // Evolution

            //- Advance z_b over the flow time step dt [s] (Exner equation
            //- with morphFactor, limited by bedCourant, then sand slide)
            //- from the current shear stress and active flags; reports
            //- the volume balance and advances the morphological time
            void advance(const scalar dt);

            //- Write the morphological state (t_m, factor in force) to
            //- <time>/uniform/bedMorphology; called at bed writes
            void writeMorphologyState() const;


        // Geometry

            //- Bed elevation at the local points of the bed patch (of the
            //- bed mesh): inverse-distance average of the surrounding
            //- faces, consistent across processors [m]
            tmp<scalarField> pointElevation() const;


        // Output

            //- Write the bed surface (patch points at z_b) with its fields
            //- to postProcessing/bed/<time>/bed.vtp (merged in parallel)
            void writeSurface() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
