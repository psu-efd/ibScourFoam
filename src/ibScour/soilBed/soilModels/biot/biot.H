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
    Foam::biot

Description
    Level 2 soil model: quasi-static Biot consolidation (u-p) with a
    run-time-selectable constitutive law (soilConstitutiveModel).
    Increments with respect to the geostatic state:

        div(sigma'(D)) - grad(p) = 0
        (n/K') dp/dt + d(div D)/dt = div((k/gamma_w) grad p)

    solved in the segregated finite-volume form of Jasak & Weller (2000)
    and Tang, Hededal & Cardiff (2015): the displacement equation with the
    implicit stiffness 2 mu + lambda and the explicit correction
    div(sigma' - (2 mu + lambda) grad D), and the pressure equation with
    the skeleton compressibility 1/M added implicitly and removed
    explicitly (stabilised coupling, contractive for saturated soil),
    iterated within each time step until the relative change of p falls
    below 'tolerance' or 'nCorrectors' is reached. Pressures in Pa, excess
    over hydrostatic; sigma' tension positive.

    Boundary conditions: top (bed) p fixed from the flow; D traction with
    the effective traction given by topTraction() (the bed shear stress,
    zero normal increment since p equals the bed pressure there); bottom
    D = 0, p zero-flux (or 'bottomPressure'); sides and structure walls
    slip (zero normal displacement) unless listed in
    'fixedDisplacementPatches' (D = 0, e.g. a pile shaft).

    The effective stress sigma'(D) includes the pore-pressure effect
    through equilibrium, so the loss of vertical effective stress is
    sigma'_zz (see soilRegion).

    The hydraulic diffusivity of the coupled problem is
    c_v = k / (gamma_w (1/M + n/K')), M the constrained modulus.

Usage
    \verbatim
    soil
    {
        type            biot;
        porosity        0.4;
        conductivity    2e-4;       // k [m/s]
        saturation      0.99;       // S_r (or bulkModulus)
        constitutive    { type linearElastic; E 5e7; nu 0.3; }
        nCorrectors     30;
        tolerance       1e-6;       // relative change of p per corrector
        aitken          yes;        // Aitken on D (default no)
        volumetricSmoothing 1;      // B-bar passes of div(D) (default 0)
        impKScale       1;          // implicit stiffness scale (default 1)
        stabilisationScale 1;       // compact-vs-Gauss stabilisation (default 1)
        solver          segregated; // | newtonKrylov | blockCoupled
        blockInnerRelTol 1e-3;      // blockCoupled inner solves (default)
        fixedDisplacementPatches (cylinder);
        // bottomPressure constant 0;
        structureMotion             // free structure (default: prescribed)
        {
            density     1800;       // [kg/m3]
            stiffness   1e7;        // estimate of the soil stiffness [N/m]
            directions  (0 0 1);    // free components (default vertical)
            relaxation  1;          // (default 1)
            gate        1e-3;       // soil change before an update (1e-3)
        }
    }
    \endverbatim

    Free structure (quasi-static): the rigid-body
    displacement d of the structure cells is an unknown of the step, found
    within the segregated correctors from the static balance
    R(d) = F_eff + F_p + (rho_w - rho_s) g V e_z = 0 along the free
    directions, with the secant update d += relaxation R / K_s: K_s starts
    from 'stiffness' and follows -dR.dd/|dd|^2 of successive correctors
    (kept between steps), so that it adapts as the soil softens (liquefied)
    or stiffens. F_eff is the
    effective traction (geostatic part interpolated to the faces, the
    increment of the soil cell), F_p the excess pore pressure of the soil
    cells, V the structure volume; z is up. In a liquefied soil without
    viscosity the structure moves without bound (no inertia); see the
    liquefaction envelope of linearElastic for the viscous resistance.

SourceFiles
    biot.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_soilBiot_H
#define Foam_soilBiot_H

#include "soilModel.H"
#include "soilConstitutiveModel.H"
#include "rigidStructure.H"
#include "Function1.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class biot Declaration
\*---------------------------------------------------------------------------*/

class biot
:
    public soilModel
{
    // Private Data

        //- Porosity
        const scalar porosity_;

        //- Hydraulic conductivity k [m/s]
        const scalar conductivity_;

        //- Pore-fluid bulk modulus K' [Pa]
        const scalar bulkModulus_;

        //- Storage coefficient n/K' [1/Pa]
        const dimensionedScalar storage_;

        //- Mobility k/gamma_w [m3 s/kg]
        const dimensionedScalar mobility_;

        //- Outer correctors per time step and displacement tolerance
        const label nCorrectors_;
        const scalar tolerance_;

        //- Aitken acceleration of the displacement correctors
        const bool aitken_;

        //- Aitken relaxation of the pressure iterates (experimental: it
        //- can overshoot at the start of a step; default off)
        const bool aitkenPressure_;

        //- Solution of the coupled step: segregated (the stabilised
        //- fixed point), newtonKrylov (Jacobian-free Newton with
        //- restarted GMRES, the segregated sweep as preconditioner)
        //- or blockCoupled (the u-p system solved together each
        //- outer iteration)
        const word solverType_;

        //- Newton-Krylov settings: Krylov subspace dimension, restarts,
        //- Newton iterations, residual tolerance (relative to the first
        //- residual of the step), finite-difference perturbation
        const label krylovDim_;
        const label krylovRestarts_;
        const label maxNewton_;
        const scalar newtonTol_;
        const scalar fdEpsilon_;

        //- Block-coupled solve: relative tolerance of the inner solves of
        //- the fixed-stress preconditioner (flexible GMRES)
        const scalar blockInnerRelTol_;

        //- Volumetric smoothing (B-bar): the volumetric part of the strain
        //- given to the law is the neighbour-averaged divergence of D,
        //- applied nSmooth times (0: off)
        const label volumetricSmoothing_;

        //- Scale of the implicit stiffness (2 mu + lambda) of the
        //- displacement equation; the stabilisation of the segregated
        //- scheme scales with it (default 1)
        const scalar impKScale_;

        //- Scale of the stabilisation of the segregated scheme: the
        //- explicit part of the implicit Laplacian is taken with the
        //- compact scheme for the fraction (1 - s) and with the Gauss
        //- gradient for the fraction s, so that at convergence the
        //- compact-minus-Gauss stabilisation is scaled by s (default 1;
        //- the implicit operator is unchanged)
        const scalar stabilisationScale_;

        //- Optional fixed pressure on the bottom patch [Pa]
        autoPtr<Function1<scalar>> bottomPressure_;

        //- Sealed top: no pore flow through the top (zero-gradient p, as
        //- under an impermeable membrane); the total traction is kept,
        //- so the effective traction carries the excess pore pressure of
        //- the top cells (element tests; default false: the top pore
        //- pressure is the flow pressure)
        const bool sealedTop_;

        //- Structure of the 'structureMotion' dictionary (owned here and
        //- attached; nullptr otherwise) and whether its balance is solved
        //- within the correctors (generalised to six
        //- degrees of freedom by rigidStructure)
        autoPtr<rigidStructure> ownStructure_;

        //- True once the reference point of the attached structure has
        //- been set (the volume centroid of the masked cells when the
        //- dictionary gave none)
        mutable bool structureRefSet_;

        //- Restart: the displacement and pore pressure were read from
        //- the time directory
        bool restored_;

        //- The stress of the restored displacement is still to be formed
        //- (once the geostatic stress is set, at the first solve)
        bool restoredStressPending_;

        //- Top (flow) pressure of the last correct [Pa]
        scalarField pTop_;

        //- Constitutive law of the skeleton
        autoPtr<soilConstitutiveModel> law_;

        //- Pore pressure [Pa]
        volScalarField p_;

        //- Displacement increment [m]
        volVectorField D_;

        //- Displacement gradient
        volTensorField gradD_;

        //- Effective stress increment [Pa]
        volSymmTensorField sigma_;

        //- Geostatic effective stress [Pa] (tension positive)
        volSymmTensorField sigma0_;

        //- Volumetric strain div(D) (old time kept for the storage term)
        volScalarField divD_;

        //- Explicit part of div(sigma')
        volVectorField divSigmaExp_;

        //- Structure contact ('structureContact' sub-dictionary):
        //- bonded (default: displacement continuity on the structure-soil
        //- faces) or frictional: a face opens when its normal effective
        //- traction is tensile (beyond tensionTolerance [Pa]) and closes
        //- again when the gap would become negative; a closed face slips
        //- when the tangential traction exceeds friction x the normal
        //- pressure + adhesion [Pa] (Coulomb): the excess is converted
        //- into a tangential slip of the structure relative to the soil
        //- (return mapping), so that the transmitted tangential force is
        //- the bonded one minus the slip and stays at the Coulomb limit.
        //- Implemented on the segregated solve as an edit of the implicit
        //- coupling (open faces) and an explicit correction field (open
        //- faces: traction-free; closed faces: the slip offset).
        const word contactType_;
        const scalar contactFriction_;
        const scalar contactAdhesion_;
        const scalar contactTension_;

        //- Active set control: the states change only when the soil
        //- iterates have settled (relative change below 'gate', default
        //- 1e-3) and at most 'maxUpdates' times per step (default 30),
        //- against the chatter of the open/closed and stick/slip switches
        const scalar contactGate_;
        const label contactMaxUpdates_;

        //- Relaxation of the slip increment of the return mapping at a
        //- state update (1: the full increment that brings the tangential
        //- force of the face to the Coulomb limit; 'slipRelaxation')
        const scalar contactFrictionRelax_;

        //- Number of faces that changed state at the last state update
        label contactChanged_;

        //- True when any processor has contact faces on processor patches
        //- (the neighbour exchange of the processor-face branch is then
        //- made by every rank: an exchange inside a rank-dependent branch
        //- deadlocks)
        bool contactProcAny_;

        //- Contact faces: internal faces between a structure cell and a
        //- soil cell (face, soil cell, structure cell, unit normal out of
        //- the structure into the soil, state 0 closed/stick, 1 open,
        //- 2 slip) and the index of each internal face in that list (-1)
        labelList contactFaces_;
        labelList contactSoilCells_;
        labelList contactStructureCells_;
        vectorField contactNormals_;
        labelList contactStates_;
        labelList faceContactIndex_;

        //- Accumulated tangential slip of the structure relative to the
        //- soil on the closed faces [m]: the transmitted tangential force
        //- is that of the bonded coupling minus the slip times the
        //- coupling coefficient; updated by the return mapping at the
        //- state updates, so that the iterations in between are linear.
        //- The return mapping starts from the converged slip of the
        //- previous step (contactSlip0_) at every update, so that the slip
        //- of an overshooting iterate is undone when the iterate corrects
        //- (an accumulated slip made the result depend on the iteration
        //- path: +14 % rotation with the dynamic solver, 2026-10-28)
        vectorField contactSlip_;
        vectorField contactSlip0_;

        //- Contact faces on processor patches with the soil cell on this
        //- side (patch, patch face, soil cell, normal, state)
        labelList contactProcPatches_;
        labelList contactProcFaces_;
        labelList contactProcSoilCells_;
        vectorField contactProcNormals_;
        labelList contactProcStates_;
        vectorField contactProcSlip_;
        vectorField contactProcSlip0_;

        //- Explicit correction of the contact faces in the displacement
        //- equation [Pa/m]
        volVectorField contactCorr_;

        //- Effective traction on the top faces [Pa]
        vectorField topTraction_;

        //- Value fraction (1: prescribed displacement) and prescribed
        //- displacement [m] of the top faces; tangential: only the
        //- tangential components are prescribed (directionMixed)
        scalarField topFraction_;
        vectorField topDisplacement_;
        bool topTangential_;


    // Private Member Functions

        //- Pore-fluid bulk modulus from the dictionary
        static scalar readBulkModulus(const dictionary& dict);

        //- Patch field types of p and D
        wordList pTypes() const;
        wordList DTypes(const dictionary& dict) const;

        //- Update the traction gradient of the top patch
        void updateTraction();

        //- Update gradD, sigma, divSigmaExp and divD from D
        void updateStress();

        //- Build the contact face lists from the structure mask
        void buildContactFaces();

        //- Update the contact states and the explicit correction field
        //- from the current iterate (impK: the implicit stiffness)
        void updateContact
        (
            const volScalarField& impK,
            const bool updateStates
        );

        //- Remove the implicit coupling of the open contact faces
        void applyContact(fvVectorMatrix& DEqn) const;

        //- Open contact faces are free surfaces of the soil cell: remove
        //- the displacement jump across them from the Gauss gradient and
        //- divergence of the soil cell (and of the structure cell)
        void correctOpenFaceStrain();

        //- Contact state of the boundary faces of the structure side
        //- (soil across), from the soil side's lists (0 when bonded)
        tmp<labelField> procContactStates() const;

        //- Newton-Krylov solve of one step: the two matrices of
        //- the corrector supply the residual of the current iterate, the
        //- displacement of the top faces is unknown with the traction
        //- condition as its equation, the Jacobian acts by finite
        //- differences, one linear segregated sweep on the increment
        //- preconditions a restarted GMRES.
        //- Returns true when converged (the residual below newtonTol
        //- times its initial value); false when the Newton iterations
        //- ran out or the line search found no decrease, with the last
        //- accepted iterate left in D and p (the caller falls back to the
        //- segregated correctors). The number of Newton iterations and
        //- the relative changes of the last step are returned through the
        //- references
        bool solveNewtonKrylov
        (
            const volScalarField& impK,
            volScalarField& invM,
            const dimensionedScalar& invMsrc,
            const volScalarField::Internal& fGen,
            const label nSteps,
            const scalar pScale,
            label& nNewton,
            scalar& changeP,
            scalar& changeD
        );

        //- Block-coupled solve of one step: outer defect-
        //- correction iterations on the residual of the current state
        //- (law, traction condition and explicit stress refreshed), each
        //- solving the coupled linear u-p system with both coupling terms
        //- and the full linearised stress (law tangent moduli) implicit by
        //- flexible restarted GMRES on the block operator (homogeneous
        //- matrices), preconditioned by one fixed-stress sweep with the
        //- Laplacian stress operator. Returns true when the relative changes fall
        //- below the tolerance within nCorrectors outer iterations
        bool solveBlockCoupled
        (
            const volScalarField& impK,
            volScalarField& invM,
            const dimensionedScalar& invMsrc,
            const volScalarField::Internal& fGen,
            const label nSteps,
            const scalar pScale,
            label& nOuter,
            scalar& changeP,
            scalar& changeD
        );

        //- Residual of the coupled equations at the current D and p
        //- (internal fields) and the displacement bFace of the top faces,
        //- in the units of the matrices: the D block in rD (3 per cell),
        //- the p block in rp and the traction block in rb (the defect of
        //- the traction condition times the face area); assembles both
        //- equations as the corrector does. The face displacement is
        //- imposed through the gradient of the traction directions, so
        //- that the residual is a function of the state alone; with
        //- bFace null the traction condition is refreshed from the state
        //- instead (the incremental fixed-point update of the corrector)
        void coupledResidual
        (
            const volScalarField& impK,
            volScalarField& invM,
            const dimensionedScalar& invMsrc,
            const volScalarField::Internal& fGen,
            const label nSteps,
            const vectorField* bFace,
            vectorField& rD,
            scalarField& rp,
            vectorField& rb
        );

        //- Defect of the traction condition on the top faces, t - n.sigma'
        //- in the traction directions (zero in the prescribed ones) [Pa]
        tmp<vectorField> tractionDefect() const;

        //- Mismatch of the traction condition on the top faces: the
        //- normalised norm of t - n.sigma' over the traction faces
        scalar tractionMismatch(const scalar pScale) const;


public:

    //- Runtime type information
    TypeName("biot");


    // Constructors

        //- Construct from the soil mesh and dictionary
        biot(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~biot() = default;


    // Member Functions

        virtual const volScalarField& p() const { return p_; }

        virtual void initialise(const scalarField& pCells);

        virtual void remap
        (
            const labelList& c0,
            const labelList& c1,
            const scalarField& w
        );

        virtual void correct(const scalarField& pTop, const label nSteps);

        //- Restart: the state fields were read from the time directory
        virtual bool restored() const { return restored_; }

        //- Write the state of the owned structure (restart)
        virtual void writeState() const;

        //- c_v = k / (gamma_w (1/M + n/K'))
        virtual scalar diffusivity() const;

        virtual bool solvesDisplacement() const { return true; }

        virtual const volSymmTensorField& sigma() const { return sigma_; }

        virtual void topTraction(const vectorField& t) { topTraction_ = t; }

        virtual void topDisplacement
        (
            const scalarField& fraction,
            const vectorField& D,
            const bool tangential
        )
        {
            topFraction_ = fraction;
            topDisplacement_ = D;
            topTangential_ = tangential;
        }

        virtual tmp<vectorField> patchEffectiveTraction
        (
            const label patchi
        ) const;

        virtual void geostatic(const scalarField& sigmaV0, const scalar k0);

        //- Set the structure mask and rebuild the contact faces
        virtual void setStructureMask(const scalarField& mask);

        //- Pressure of the structure cells as the image of the soil
        //- pressure next to them (zero normal gradient at the impermeable
        //- structure face), one value per masked cell
        tmp<scalarField> maskedPressureImage() const;

        //- True for the frictional (unilateral) structure contact
        bool contactFrictional() const
        {
            return contactType_ == "frictional";
        }

        //- Number of open and slipping contact faces (all processors)
        void contactCounts(label& nOpen, label& nSlip) const;

        //- Geostatic effective stress [Pa]
        const volSymmTensorField& sigma0() const { return sigma0_; }

        //- Displacement increment [m]
        const volVectorField& D() const { return D_; }

        //- Effective force [N] of the soil on the structure cells: the
        //- cell effective stress of the soil side dotted with the face
        //- normal pointing out of the structure, summed over the
        //- soil-structure faces (internal and processor faces)
        virtual vector structureForce() const;

        //- Loads on the structure cells [N]: the effective traction Feff
        //- (geostatic stress interpolated to the soil-structure faces plus
        //- the stress increment of the soil cell) and the excess pore
        //- pressure Fp of the soil cells; the structure volume V [m3]
        void structureLoads(vector& Feff, vector& Fp, scalar& V) const;

        //- The loads with their moments Meff, Mp [N m] about xRef, and
        //- the volume centroid xV [m] of the structure cells
        void structureLoads
        (
            vector& Feff,
            vector& Fp,
            scalar& V,
            vector& Meff,
            vector& Mp,
            point& xV,
            const point& xRef
        ) const;

        //- Resultant [N] of Feff, Fp, the hydrostatic buoyancy and, for a
        //- free structure, its weight and external force (z up)
        vector structureResultant() const;

        //- Resultant moment [N m] about the structure reference point
        //- (as structureResultant)
        vector structureResultantMoment() const;

        //- Reference point of the attached structure, set from the
        //- volume centroid of the masked cells when the dictionary gave
        //- none
        void ensureStructureRefPoint() const;

        //- True for an attached free structure (its balance is solved
        //- within the correctors)
        bool structureFree() const
        {
            return structure_ && structure_->free();
        }

        //- The constitutive law
        const soilConstitutiveModel& law() const { return *law_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
