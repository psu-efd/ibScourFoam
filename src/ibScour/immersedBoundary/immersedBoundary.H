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
    Foam::immersedBoundary

Description
    Conservative, sharp immersed boundary on a fixed fvMesh.

    The immersed solids (ibSolid) give every cell a signed distance and a
    normal. Cells are classified as

        - solid: the cell centre lies inside a solid;
        - IB:    a fluid cell with at least one solid face neighbour;
        - live:  any other fluid cell.

    Faces with a solid cell on either side are blocked: their flux is zero
    and they are removed from the pressure Laplacian, so mass is conserved
    exactly in the fluid. Solid cells are decoupled and pinned. The velocity
    of IB cells is set by the wall model (ibWallModel) through the momentum
    equation and after the pressure correction; the projected flux is never
    rebuilt from the velocity.

    Cut faces ('cutFaces true', height-field beds only): the bed is a level
    set phi = z - z_b at the mesh points. Every face gets its open (fluid)
    area fraction alpha_f and every cell its fluid volume fraction theta
    (cutGeometry, FAVOR-type apertures). Cells with theta < cutThetaMin are
    solid; the other cells cut by the bed are IB cells. Fluxes and the
    pressure Laplacian are weighted by alpha_f (instead of 0/1 blocking), so
    the flow sees the inclined bed instead of a staircase of blocked faces;
    the diffusive coefficients of cut faces (scaleCutDiffusion) and the
    pressure gradient of IB cells are weighted in the same way. The wall
    element of a cut cell follows from closure, A_w n_w = sum_f (1 -
    alpha_f) S_f; its wall distance is that of the cell centre, or, for a
    centre below the bed, that of the fluid part, y = theta V/(2 A_w).
    Other solids keep the blocking treatment.

    Centroid-consistent cut cells (optional 'cutCentroids' sub-dictionary):
    the velocity and pressure of a cut cell are fluid-volume averages, i.e.
    values at the centroid of its fluid part, not at the cell centre. With
    'diffusion' and 'pressure' the diffusive and pressure-Laplacian
    coefficients of a face use the normal distance between the fluid
    centroids, d_c, instead of that between the cell centres, d; d_c is
    bounded below by minDistanceFraction*d (coefficients at most 1/0.5 = 2
    times the uncut value). With 'wallDistance' the wall distance of a cut
    cell is that of its fluid centroid (> 0 for theta > 0); the wall model
    keeps its own bound (sampleBelow).

    Pressure gradient of the momentum equation and the velocity correction
    ('pressureGradient'): 'gauss' (default) is the Gauss gradient with the
    cell value on the blocked part of the faces of IB cells. 'faceIB' and
    'face' reconstruct it in the IB cells (or in all fluid cells) from the
    open-face gradient fluxes of the pressure equation,
        grad p = inv(sum alpha_f |S_f| n_f n_f) sum n_f alpha_f r_f snGrad(p) |S_f|,
    exact for a linear p, so that the cell velocity correction is
    consistent with the face flux correction (no solid-cell pressure).

    The object is a MeshObject: one instance per mesh, constructed on first
    use with immersedBoundary::New(mesh), and shared with fvOptions and
    function objects.

Usage
    Dictionary constant/immersedBoundaryDict:
    \verbatim
    // Optional evolving bed (see heightFieldBed)
    bed
    {
        type        heightField;
        ...
    }

    solids                          // may be empty when a bed is given
    {
        bed
        {
            type        plane;
            point       (0 0 0);
            normal      (0 0 1);
        }
    }

    wallModel
    {
        samplingFactor  1;          // sampling distance / cell extent
        wallLaw
        {
            type        spalding;
        }
    }

    // Optional
    cutFaces        false;          // apertures for the height-field bed
    cutThetaMin     0.05;           // cut cells with theta below are solid
    cutCentroids                    // fluid centroids of cut cells
    {
        diffusion       false;      // diffusive coefficients from d_c
        pressure        false;      // pressure Laplacian from d_c
        wallDistance    false;      // wall distance of the fluid centroid
        minDistanceFraction 0.5;    // d_c >= minDistanceFraction*d
    }
    pressureGradient gauss;         // gauss | faceIB | face (see grad)
    wallSink        isotropic;      // isotropic | tangential (stress mode)
    rAUInterpolation linear;        // linear | harmonic (faces of IB cells)
    solidThreshold  0;              // solid if distance < threshold*h
    pressureSolve   subset;         // subset (default) | penalty
    minFluidRegionFraction 0.01;    // smaller fluid regions become solid
    writeFields     true;           // Gamma, ibDistance, ...
    \endverbatim

SourceFiles
    immersedBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_immersedBoundary_H
#define Foam_immersedBoundary_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "IOdictionary.H"
#include "ibSolid.H"
#include "PtrList.H"
#include "fvMeshSubset.H"
#include "bedModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class ibWallModel;
class pimpleControl;

/*---------------------------------------------------------------------------*\
                       Class immersedBoundary Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundary
:
    public MeshObject<fvMesh, TopologicalMeshObject, immersedBoundary>
{
public:

    // Public Data Types

        //- Cell status
        enum cellStatus : label
        {
            SOLID = 0,
            IB = 1,
            LIVE = 2
        };


private:

    // Private Typedefs

        typedef MeshObject
        <
            fvMesh,
            TopologicalMeshObject,
            immersedBoundary
        > MeshObject_type;


    // Private Data

        //- Settings (constant/immersedBoundaryDict)
        IOdictionary dict_;

        //- The immersed solids
        PtrList<ibSolid> solids_;

        //- Evolving sediment bed (optional)
        autoPtr<bedModel> bed_;

        //- The bed has moved enough to require a new classification
        bool bedMoved_;

        //- Cells that changed from solid to fluid at the last
        //- classification (to be initialised by the solver)
        labelList freshCells_;

        //- Cells that changed from fluid to solid at the last
        //- classification (buried by the bed or entered by a body)
        labelList buriedCells_;

        //- Cells whose centre is closer than solidThreshold*h to a solid
        //- surface are treated as solid (default 0: a cell is solid when
        //- its centre is inside a solid, to a round-off tolerance). The
        //- default wall model samples u_tau at the live neighbours and does
        //- not depend on the IB-cell wall distance; with 'uTauFrom cell' a
        //- threshold of about 0.2 keeps the wall law out of the viscous
        //- sublayer (tau = nu U/y diverges for y -> 0 on a moving bed).
        scalar solidThreshold_;

        //- Fluid regions smaller than this fraction of all fluid cells
        //- (isolated pockets enclosed by solids) are made solid
        scalar minPocketFraction_;

        //- Scale of the diagonal penalty holding the solid-cell pressure,
        //- relative to a typical fluid diagonal
        scalar penaltyFactor_;

        //- Cell status (SOLID, IB, LIVE)
        labelList status_;

        //- Index of the nearest solid per cell
        labelList nearestSolid_;

        //- IB cells
        labelList ibCells_;

        //- Wall area inside each IB cell: sum over its faces to solid
        //- cells of |S_f . n| (projected on the wall normal) [m2]
        scalarField ibWallArea_;

        //- Solid cells
        labelList solidCells_;

        //- Signed distance to the nearest solid surface [m]
        volScalarField distance_;

        //- Wall distance of the mesh walls alone (wallDist before the
        //- immersed solids are included; -GREAT where no mesh wall
        //- reaches), kept so that the distance can be rebuilt after the
        //- solids or the bed moved
        mutable scalarField meshWallDistance_;

        //- Unit normal of the nearest solid surface (into the fluid)
        volVectorField normal_;

        //- Cell indicator for post-processing: 1 live, 0.5 IB, -1 solid
        volScalarField gamma_;

        //- Face mask: 1 open, 0 blocked
        surfaceScalarField faceMask_;

        //- Pressure-Laplacian mask: 0 on faces between a fluid and a solid
        //- cell, 1 elsewhere (solid cells stay coupled to each other)
        surfaceScalarField laplacianMask_;

        //- Wall model for the IB cells
        autoPtr<ibWallModel> wallModel_;

        //- Solve the pressure on a fluid-only mesh subset (true, default)
        //- or on the full mesh with penalised solid cells (false)
        bool pressureSubset_;

        //- Fluid-only mesh subset (non-solid cells) for the pressure
        mutable autoPtr<fvMeshSubset> fluidSubset_;

        //- Base cell -> subset cell (-1 for solid cells)
        mutable labelList subsetCell_;

        //- Cut-face apertures for the height-field bed
        bool cutFaces_;

        //- Cells cut by the bed with a smaller fluid fraction are solid
        scalar cutThetaMin_;

        //- Fluid volume fraction of the cells (1 without cut faces)
        volScalarField theta_;

        //- Geometric open-area fraction of all faces (cut faces)
        scalarField faceAlpha_;

        //- Cut cells: diffusive coefficients from the fluid centroids
        bool centroidDiffusion_;

        //- Cut cells: pressure Laplacian from the fluid centroids
        bool centroidPressure_;

        //- Cut cells: wall distance of the fluid centroid
        bool centroidWallDistance_;

        //- Lower bound of the centroid distance across a face, as a
        //- fraction of the cell-centre distance
        scalar minDistanceFraction_;

        //- Fluid centroid minus cell centre (zero for cells not cut) [m]
        vectorField fluidCentreOffset_;

        //- Wall element of the cut cells from the cut geometry: area vector
        //- A_w n_w pointing into the solid (zero for cells not cut) [m2]
        vectorField cutWallVec_;

        //- Centre of the wall element minus cell centre (cut cells) [m]
        vectorField cutWallOffset_;

        //- Cut cells whose wall-element distance was floored at the last
        //- classification (non-planar bed inside the cell)
        label nFlooredCutDistance_;

        //- Cut cells whose least-squares gradient stencil has an
        //  unresolved direction (last reconstruction)
        mutable label nTruncatedGrad_;

        //- Largest |Sw - closure|/A_w over the cut bed cells at the last
        //- classification (zero for a planar bed inside every cell)
        scalar maxWallClosureError_;

        //- Cut faces: explicit non-orthogonal (skew) correction of the
        //- diffusion flux on the faces of cut cells (cutCentroids
        //- skewCorrection)
        bool skewCorrection_;

        //- Cut faces: interpolation to the centroid of the open part of the
        //- faces of cut cells (cutCentroids faceInterpolation)
        bool faceInterpolation_;

        //- Centroid of the open part of each face minus the face centre
        //- (zero for faces not cut) [m]
        vectorField faceOpenOffset_;

        //- Faces of cut cells between two fluid cells (internal and coupled
        //- boundary faces), the faces the two corrections act on
        surfaceScalarField cutFaceFlag_;

        //- Owner interpolation weights along the fluid-centroid line to the
        //- open-face centroid (mesh weights on other faces)
        surfaceScalarField cutWeights_;

        //- Open-face centroid minus the point of the centroid line on the
        //- face, x_fo - x_line (zero on other faces) [m]
        surfaceVectorField cutInterpOffset_;

        //- Skew-correction vector k_f = n_f - Delta*/d_c of the faces of
        //- cut cells (zero on other faces) [-]
        surfaceVectorField cutSkewVec_;

        //- Cell-centre distance over fluid-centroid distance across each
        //- face (1 unless a cut cell is involved)
        surfaceScalarField cutDeltaRatio_;

        //- Pressure gradient: 0 gauss, 1 face-reconstructed in the IB
        //- cells, 2 face-reconstructed in all fluid cells
        label gradMode_;

        //- Stress-mode wall sink on the tangential velocity only (the
        //- normal part is returned explicitly); false: all components
        bool tangentialSink_;

        //- Harmonic interpolation of rAU on the faces of IB cells in the
        //- pressure Laplacian (false: linear)
        bool harmonicRAU_;

        //- Remove the explicit non-orthogonal correction of the Laplacians
        //- on the faces of IB and solid cells (maskNonOrthCorrection,
        //- default true): built from cell gradients that reach into the
        //- pinned solid cells, it feeds back step by step on non-orthogonal
        //- meshes (IB cells beside freshly buried cells on an O-grid ring
        //- ran to 1.9 m/s in a 0.17 m/s wave, 2026-10-31)
        bool maskNonOrthCorr_;


    // Private Member Functions

        //- Read the solids from the dictionary
        void readSolids();

        //- Distances, normals and status of all cells
        void classify();

        //- Make isolated fluid pockets solid
        void removeFluidPockets();

        //- Apertures and volume fractions from the bed level set
        void computeCutGeometry();

        //- Face mask from the cell status
        void makeFaceMask();

        //- Distance ratio of the faces from the fluid centroids
        void makeCutDeltaRatio();

        //- Interpolation weights, open-face offsets and skew vectors of the
        //- faces of cut cells (after makeCutDeltaRatio)
        void makeCutFaceData();

        //- Gauss gradient, blocked parts of IB-cell faces at the cell value
        tmp<volVectorField> gaussGrad(const volScalarField& p) const;

        //- Gradient reconstructed from the open-face gradient fluxes
        tmp<volVectorField> faceGrad(const volScalarField& p) const;

        //- Least-squares reconstruction of the (masked) face flux F into
        //- the cell vectors g of the IB cells (mode 1) or of all fluid
        //- cells (mode 2); the other cells keep their value
        void openFaceReconstruct
        (
            const surfaceScalarField& F,
            const label mode,
            vectorField& g
        ) const;

        //- constrainU with an optional density per cell for the wall sink
        void constrainUImpl
        (
            fvVectorMatrix& UEqn,
            const bool addWallSink,
            const scalarField* rhoPtr
        ) const;

        //- Report the classification
        void report() const;

        //- No copy construct
        immersedBoundary(const immersedBoundary&) = delete;

        //- No copy assignment
        void operator=(const immersedBoundary&) = delete;


public:

    //- Runtime type information
    TypeName("immersedBoundary");


    // Constructors

        //- Construct from mesh, reading constant/immersedBoundaryDict
        explicit immersedBoundary(const fvMesh& mesh);


    //- Destructor
    virtual ~immersedBoundary();


    // Member Functions

        // Access

            //- The mesh
            const fvMesh& mesh() const noexcept
            {
                return MeshObject_type::mesh();
            }

            //- Settings dictionary
            const dictionary& dict() const noexcept
            {
                return dict_;
            }

            //- The immersed solids
            const PtrList<ibSolid>& solids() const noexcept
            {
                return solids_;
            }

            //- Cell status
            const labelList& status() const noexcept
            {
                return status_;
            }

            //- Nearest solid per cell
            const labelList& nearestSolid() const noexcept
            {
                return nearestSolid_;
            }

            //- IB cells
            const labelList& ibCells() const noexcept
            {
                return ibCells_;
            }

            //- Solid cells
            const labelList& solidCells() const noexcept
            {
                return solidCells_;
            }

            //- Geometric open-area fraction of the faces (empty unless cut
            //- faces are active)
            const scalarField& faceAlpha() const noexcept
            {
                return faceAlpha_;
            }

            //- True if the bed uses cut-face apertures
            bool cutFaces() const noexcept
            {
                return cutFaces_;
            }

            //- Fluid volume fraction of the cells
            const volScalarField& theta() const noexcept
            {
                return theta_;
            }

            //- Wall area inside each IB cell [m2]
            const scalarField& ibWallArea() const noexcept
            {
                return ibWallArea_;
            }

            //- Signed distance [m]
            const volScalarField& distance() const noexcept
            {
                return distance_;
            }

            //- Surface normal of the nearest solid (into the fluid)
            const volVectorField& normal() const noexcept
            {
                return normal_;
            }

            //- Face mask (1 open, 0 blocked)
            const surfaceScalarField& faceMask() const noexcept
            {
                return faceMask_;
            }

            //- Wall model
            const ibWallModel& wallModel() const;

            //- Wall model (non-const)
            ibWallModel& wallModel();

            //- True if any solid moves or the bed requires a new
            //- classification
            bool moving() const;

            //- True if an evolving bed is present
            bool hasBed() const noexcept
            {
                return bool(bed_);
            }

            //- The bed model (hasBed() must be true)
            bedModel& bed()
            {
                return *bed_;
            }

            //- The bed model (hasBed() must be true)
            const bedModel& bed() const
            {
                return *bed_;
            }

            //- Cells that changed from solid to fluid at the last update
            const labelList& freshCells() const noexcept
            {
                return freshCells_;
            }

            //- Cells that changed from fluid to solid at the last update
            const labelList& buriedCells() const noexcept
            {
                return buriedCells_;
            }


        // Geometry

            //- Cell extent along the unit direction n:
            //- h = 2 V / sum_f |n . S_f|
            scalar cellExtent(const label celli, const vector& n) const;

            //- Hydrodynamic force [N] of the fluid on solid solidi from the
            //- pressure pIB [Pa] and the density rhoIB [kg/m3] of the IB
            //- cells (indexed as ibCells()): sum of (-pIB n + rhoIB tau_w)
            //- A_w with tau_w the kinematic wall shear stress of the cell
            //- (first order: the cell pressure at the wall), summed over
            //- the processors
            vector solidForce
            (
                const label solidi,
                const scalarField& pIB,
                const scalarField& rhoIB
            ) const;

            //- Single-phase form: the kinematic pressure field p and the
            //- wall shear stress of the IB cells times rho
            vector solidForce(const label solidi, const scalar rho) const;

            //- Hydrodynamic moment [N m] of the fluid on solid solidi about
            //- the point xRef (the wall point of each IB cell is its
            //- centre moved to the wall along the normal)
            vector solidMoment
            (
                const label solidi,
                const scalarField& pIB,
                const scalarField& rhoIB,
                const point& xRef
            ) const;

            //- Single-phase form of the moment
            vector solidMoment
            (
                const label solidi,
                const scalar rho,
                const point& xRef
            ) const;

            //- Velocity of the nearest solid at the given points
            //- (one point per listed cell)
            tmp<vectorField> solidVelocity
            (
                const labelUList& cells,
                const pointField& points
            ) const;


        // Update

            //- Re-classify the cells if solids or the bed moved
            void update();

            //- Advance the bed (Exner) over the flow time step [s]
            void evolveBed(const scalar dt);

            //- Initialise the velocity of fresh cells (solid -> fluid)
            void initialiseFreshCells(volVectorField& U) const;

            //- Two-phase variant: the fresh cells also take the volume
            //- fraction and the piezometric pressure of the fluid cells
            //- next to them, and the water of the cells that became
            //- solid (alpha V) is moved to the fluid cell above them
            //- (what does not fit is reported as lost)
            void initialiseFreshCells
            (
                volVectorField& U,
                volScalarField& alpha,
                volScalarField& p_rgh
            ) const;

            //- Correct the wall-distance field of the turbulence models so
            //- Wall distance of the turbulence models: the mesh walls and
            //- the immersed solids and bed. Rebuilt from the mesh-wall
            //- distance every call, so it follows a moving bed; cells that
            //- no mesh wall reaches (a mesh without wall patches) take the
            //- immersed distance alone.
            void correctWallDistance() const;


        // Constraints

            //- Zero the flux on blocked faces
            void maskFlux(surfaceScalarField& phi) const;

            //- Faces with an IB or solid cell on either side (1), else 0
            tmp<surfaceScalarField> ibTouchMask() const;

            //- Whether the non-orthogonal correction is removed on those faces
            bool maskNonOrthCorrection() const noexcept
            {
                return maskNonOrthCorr_;
            }

            //- Remove the stored explicit non-orthogonal correction of a
            //- Laplacian (fluxRequired fields) on the masked faces: from
            //- the source and from the flux correction alike
            template<class Type>
            void maskNonOrthCorrection
            (
                fvMatrix<Type>& eqn,
                const surfaceScalarField& mask
            ) const;

            //- Return rAU interpolated to the faces with the faces between
            //- fluid and solid cells removed, for the pressure Laplacian
            tmp<surfaceScalarField> maskedLaplacianCoeff
            (
                const volScalarField& rAU
            ) const;

            //- True if the pressure is solved on the fluid-only subset
            bool pressureSubset() const noexcept
            {
                return pressureSubset_;
            }

            //- Solve the pressure on the fluid-only subset: same equation
            //- as on the full mesh with the solid cells removed, so the
            //- linear solvers see exactly a body-fitted problem. Maps the
            //- pressure back and sets the conservative flux phi.
            void solvePressure
            (
                volScalarField& p,
                const surfaceScalarField& rAUf,
                const surfaceScalarField& phiHbyA,
                surfaceScalarField& phi,
                pimpleControl& pimple,
                const label pRefCell,
                const scalar pRefValue
            ) const;

            //- Print the cell and the faces of an empty pressure row
            //- (diagnostics of solvePressure)
            void reportEmptyRow
            (
                const label celli,
                const surfaceScalarField& rAUf
            ) const;

            //- Hold the pressure of the solid cells, which are decoupled
            //- from the fluid, at its current value with a diagonal penalty
            //- (solid cells remain coupled among themselves, which keeps the
            //- matrix well conditioned for multigrid solvers)
            void fixSolidCells(fvScalarMatrix& pEqn) const;

            //- Fix the velocity in IB cells (wall model) and solid cells
            //- (solid velocity) in a momentum equation. Stress mode: the
            //- wall sink is added only with addWallSink (once per equation;
            //- false for a copy of an equation that already has it)
            void constrainU
            (
                fvVectorMatrix& UEqn,
                const bool addWallSink = true
            ) const;

            //- Reset the velocity in IB (velocity mode) and solid cells
            //- after the pressure correction
            void correctU(volVectorField& U) const;

            //- Wall-sink term of the momentum equation per unit volume in
            //- the IB cells, c (U - U_w)/V (tangential part only with
            //- wallSink tangential), as added by constrainU (stress mode;
            //- zero otherwise). Diagnostics.
            tmp<volVectorField> wallSinkTerm(const volVectorField& U) const;

            //- Remove the coupling of fluid cells to solid cells in an
            //- equation (zero gradient towards the solids). The blocked
            //- faces carry no flux, so only the diffusive part is removed.
            template<class Type>
            void decoupleSolid(fvMatrix<Type>& eqn) const;

            //- Cut faces: weight the diffusive coefficients of the faces
            //- between two fluid cells by their open fraction alpha_f
            //- (removes (1 - alpha_f) gamma_f |S_f| delta_f, the part of
            //- the implicit Laplacian of -laplacian(gamma, psi) through the
            //- closed part of the face). No-op without cut faces.
            template<class Type>
            void scaleCutDiffusion
            (
                fvMatrix<Type>& eqn,
                const volScalarField& gamma
            ) const;

            //- Pressure gradient of the momentum equation (pressureGradient:
            //- Gauss with zero normal gradient on the blocked faces of IB
            //- cells, or reconstructed from the open-face gradient fluxes)
            tmp<volVectorField> grad(const volScalarField& p) const;

            //- Reconstruct a cell vector field from face fluxes (normal
            //- component times |S_f|) over the open faces only: the
            //- least-squares form of fvc::reconstruct with the blocked
            //- faces removed, so that a flux balance on the open faces
            //- (gravity against the pressure gradient of a fluid at rest)
            //- is balanced in the IB cells too. Solid cells are zero.
            tmp<volVectorField> reconstruct
            (
                const surfaceScalarField& ssf
            ) const;

            //- Volume of the phase alpha in the non-solid cells [m3] (the
            //- solid cells carry the ghost extension of correctSolidAlpha)
            scalar fluidVolume(const volScalarField& alpha) const;

            //- alpha with the solid cells set to zero (a field for
            //- function objects that integrate alpha, e.g. interfaceHeight)
            void fluidOnly
            (
                const volScalarField& alpha,
                volScalarField& result
            ) const;

            //- Set the volume fraction of the solid cells from the fluid
            //- cells next to them (layer by layer, zero normal gradient),
            //- so that interface normals and compression fluxes next to the
            //- immersed boundary see a 90 degree contact
            void correctSolidAlpha(volScalarField& alpha) const;

            //- constrainU for the density-weighted momentum equation
            //- (ddt(rho, U) form): the wall sink is rho c (U - U_w)
            void constrainU
            (
                fvVectorMatrix& UEqn,
                const volScalarField& rho,
                const bool addWallSink = true
            ) const;

            //- Cut faces: the skew correction switch
            bool skewCorrection() const noexcept
            {
                return skewCorrection_;
            }

            //- Cut faces: the open-face interpolation switch
            bool faceInterpolation() const noexcept
            {
                return faceInterpolation_;
            }

            //- Least-squares gradient of a field at the fluid centroids over
            //- the fluid face neighbours, for the cells that touch a cut
            //- face (zero elsewhere)
            template<class Type>
            tmp
            <
                GeometricField
                <
                    typename outerProduct<vector, Type>::type,
                    fvPatchField,
                    volMesh
                >
            >
            cutLeastSquaresGrad
            (
                const GeometricField<Type, fvPatchField, volMesh>& vf
            ) const;

            //- Explicit non-orthogonal correction of the diffusion flux on
            //- the faces of cut cells, per unit volume, to be added to the
            //- right-hand side of the momentum equation:
            //- (1/V) sum_f +- gamma_f alpha_f |S_f| k_f . grad(U)_f
            tmp<volVectorField> cutSkewCorrection
            (
                const volVectorField& U,
                const volScalarField& gamma
            ) const;

            //- Explicit correction of the convection flux on the faces of cut
            //- cells, per unit volume, for the right-hand side of the
            //- momentum equation: -(1/V) sum_f +- phi_f (U_fo - U_f), with
            //- U_fo the value at the open-face centroid and U_f the linear
            //- (mesh-weight) face value the scheme uses (faceInterpolation)
            tmp<volVectorField> cutConvectionCorrection
            (
                const surfaceScalarField& phi,
                const volVectorField& U
            ) const;

            //- Interpolate a scalar field to the faces with the open-face
            //- centroid weights on the faces of cut cells (mesh weights
            //- elsewhere)
            tmp<surfaceScalarField> cutInterpolate
            (
                const volScalarField& vf
            ) const;

            //- Replace the flux of a cell vector field through the faces of
            //- cut cells by the flux of its value at the open-face centroid
            //- (centroid weights plus the tangential gradient term)
            void correctCutFlux
            (
                const volVectorField& vf,
                surfaceScalarField& phi
            ) const;


        // Diagnostics

            //- Report the flux through blocked faces and the boundary
            //- velocity flux into solid cells
            void checkFlux
            (
                const surfaceScalarField& phi,
                const volVectorField& U
            ) const;


        // MeshObject

            //- Topological changes are not supported
            virtual bool movePoints()
            {
                return false;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "immersedBoundaryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
