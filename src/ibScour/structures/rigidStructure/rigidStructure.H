/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Xiaofeng Liu, Penn State
-------------------------------------------------------------------------------
License
    This file is part of ibScourFoam, distributed under the GNU General
    Public License v3 (see the LICENSE file of OpenFOAM).

Class
    Foam::rigidStructure

Description
    A rigid structure (pile, footing, pile group, pipe, block) with one
    rigid-body state shared by every domain that sees it: the flow
    (immersed solid), the soil (masked cells) and the bed.

    State: translation d of the reference point x_ref and small rotation
    theta (rotation vector [rad]); a material point x displaces by
        D(x) = d + theta ^ (x - x_ref).
    Loads are accumulated by the owners as a force F and a moment M about
    x_ref (soil effective traction and pore pressure, buoyancy and weight,
    hydrodynamic force, external loads from the dictionary).

    Motion solvers:
    - prescribed (default when no 'free' entries): the owners set the
      state (Function1 'displacement' and 'rotation' of time, else zero);
    - quasiStatic: the static balance of the free degrees of freedom is
      solved by the owner with a secant stiffness matrix K (-dR/dq of the
      six degrees of freedom), K dq = relaxation R, started as the
      diagonal of the estimates 'stiffness' [N/m] and
      'rotationalStiffness' [N m/rad] and bounded to 1e-2..1e2 of them (a
      secant far above the estimate freezes the degree of freedom: the
      earlier component secant, bounded only relative to its previous
      value, ran away to 1e15 on the pier cases). 'update inner'
      (default): the structure is updated within the owner's correctors
      once the soil iterate has settled to 'gate', with the component
      secant of successive updates (kept when not positive, within a
      factor 100 per update); the reaction of a settled iterate lags the
      step, and the stiff bias of that secant is the under-relaxation the
      joint iteration needs. 'update outer': one update per solve from
      the converged reaction, the secant spanning the solves; K is then
      updated by the Broyden rank-one formula in scaled variables (step
      over load/stiffness, residual over the load), which learns the
      coupling of the degrees of freedom, from steps made by a residual
      above 'stiffnessResidualFloor' (default 1e-2 of the load scale;
      below it the converged reaction is within the soil accuracy), with
      the off-diagonals within the geometric mean of the diagonals and a
      diagonal fallback when the step does not decrease the residual
      form. The outer update converges a pile in a few solves at a
      fraction of the cost (no capped soil solve); the inner update is
      needed by the dynamic solver and the contact tests. This is the
      generalisation of the earlier vertical balance;
    - dynamic ('motionSolver dynamic'): the equations of motion of the free
      degrees of freedom, m a = F and I alpha = M with the mass (plus
      'addedMass'), the principal 'inertia' about the reference point (plus
      'addedInertia'), viscous 'damping' and 'rotationalDamping', are
      integrated by the Newmark scheme (beta 0.25, gamma 0.5: trapezoidal,
      unconditionally stable, no numerical damping) implicitly within the
      owner's correctors: each update solves the linearised dynamic balance
      with the same secant stiffness as the quasi-static solver as the
      tangent of the soil reaction. The soil itself stays quasi-static
      (massless); the structure carries the inertia, so the model is the
      lumped mass on the foundation springs of the soil solution.

    Mass: 'density' [kg/m3] with the volume reported by the owner (the
    masked soil cells and, later, the solid above the bed), or 'mass' [kg]
    and 'volume' [m3] given explicitly. The buoyancy of the submerged
    volume and the weight act at the volume centroid and at
    'centreOfMass' (default: the volume centroid).

Usage
    \verbatim
    structureMotion                     // in the biot dictionary, or
    {                                   // structures { <name> { ... } }
        density         1800;           // [kg/m3] (or mass + volume)
        freeTranslations (0 0 1);       // components solved (0/1)
        freeRotations   (0 0 0);        // (default none)
        referencePoint  (0 0 -1);       // x_ref [m] (default: volume
                                        // centroid of the masked cells)
        stiffness       1e7;            // estimate [N/m]
        rotationalStiffness 1e6;        // estimate [N m/rad]
        relaxation      1;
        gate            1e-3;           // soil change before an update
        externalForce   constant (0 0 -1e4);   // Function1 of time [N]
        externalForcePoint (0 0 0);     // its point of action (default
                                        // x_ref)
        externalMoment  constant (0 0 0);      // Function1 of time [N m]
        // backward compatible: 'directions' = freeTranslations
        motionSolver    dynamic;        // quasiStatic (default) | dynamic
        inertia         (1 1 1);        // [kg m2] about x_ref (dynamic,
                                        // free rotations)
        addedMass       (0 0 0);        // [kg]
        addedInertia    (0 0 0);        // [kg m2]
        damping         (0 0 0);        // [N s/m]
        rotationalDamping (0 0 0);      // [N m s/rad]
        newmarkBeta     0.25;
        newmarkGamma    0.5;
        dynamicStartTime 0.125;         // steps ending after it are dynamic;
                                        // before: quasi-static placement,
                                        // then dynamics from rest in balance
                                        // (default: dynamic from the start,
                                        // released from rest)
    }
    \endverbatim

SourceFiles
    rigidStructure.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_rigidStructure_H
#define Foam_rigidStructure_H

#include "dictionary.H"
#include "Function1.H"
#include "pointField.H"
#include "vectorField.H"
#include "scalarMatrices.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;

/*---------------------------------------------------------------------------*\
                       Class rigidStructure Declaration
\*---------------------------------------------------------------------------*/

class rigidStructure
{
    // Private Data

        //- Name
        word name_;

        //- Density [kg/m3] (0 with an explicit mass)
        scalar density_;

        //- Explicit mass [kg] and volume [m3] (negative: from the density
        //- and the volume reported by the owner)
        scalar mass_;
        scalar volume_;

        //- Free translations and rotations (components 0 or 1)
        vector freeTranslations_;
        vector freeRotations_;

        //- Reference point [m] and whether the dictionary gave it
        point refPoint_;
        bool refPointGiven_;

        //- Centre of mass [m] and whether the dictionary gave it
        point centreOfMass_;
        bool centreOfMassGiven_;

        //- External loads
        autoPtr<Function1<vector>> externalForce_;
        autoPtr<Function1<vector>> externalMoment_;
        point externalForcePoint_;
        bool externalForcePointGiven_;

        //- Loads applied by an owner (e.g. the hydrodynamic force and
        //- moment of the flow), added to the external loads
        vector appliedForce_;
        vector appliedMoment_;

        //- Prescribed motion (no free degrees of freedom)
        autoPtr<Function1<vector>> displacementFunction_;
        autoPtr<Function1<vector>> rotationFunction_;

        //- State: translation of the reference point [m] and rotation
        //- vector [rad]
        vector d_;
        vector theta_;

        //- Quasi-static solver: secant stiffnesses per component [N/m],
        //- [N m/rad], relaxation, gate, and the residual and step of the
        //- last update
        vector kT_;
        vector kR_;

        //- Secant stiffness matrix of the six degrees of freedom
        //- (translations 0-2, rotations 3-5), the estimates it starts
        //- from (its bounds) and the step floor of its update
        scalarSquareMatrix K_;
        vector kT0_;
        vector kR0_;
        scalar stiffnessResidualFloor_;
        scalar relaxation_;
        scalar gate_;
        vector RF_;
        vector RM_;
        vector dF_;
        vector dM_;
        label nUpdates_;

        //- Motion solver of a free structure: quasiStatic or dynamic
        word motionSolver_;

        //- Update of a free structure by the owner: inner (within the
        //- correctors) or outer (once per solve)
        word update_;

        //- Dynamic solver: principal inertia about the reference point
        //- [kg m2], added mass [kg] and inertia [kg m2], viscous damping
        //- [N s/m] and [N m s/rad], Newmark parameters
        vector inertia_;
        vector addedMass_;
        vector addedInertia_;
        vector damping_;
        vector rotationalDamping_;
        scalar newmarkBeta_;
        scalar newmarkGamma_;

        //- Dynamic solver: steps ending after dynamicStartTime are dynamic;
        //- before it the structure is placed quasi-statically and the
        //- dynamics then start from rest in balance (a structure resting
        //- on the bed) instead of from the release acceleration
        scalar dynamicStartTime_;
        bool dynamicActive_;
        label nStaticSteps_;

        //- Dynamic state: velocity [m/s], angular velocity [rad/s],
        //- acceleration [m/s2], angular acceleration [rad/s2]; the state
        //- at the start of the step and the step size [s]
        vector v_;
        vector omega_;
        vector a_;
        vector alpha_;
        vector d0_;
        vector theta0_;
        vector v0_;
        vector omega0_;
        vector a0_;
        vector alpha0_;
        scalar dt_;
        bool dynamicStarted_;


    // Private Member Functions

        //- Degree of freedom i (0-5) is free
        bool freeDof(const label i) const
        {
            return (i < 3 ? freeTranslations_[i] : freeRotations_[i - 3]) > 0.5;
        }

        //- Reset the stiffness matrix to the diagonal kT_, kR_
        void resetStiffness();

        //- Broyden update of the stiffness matrix from the residuals of
        //- this update and the step and residuals of the last one
        //- (quasi-static and dynamic)
        void updateStiffness
        (
            const vector& RF,
            const vector& RM,
            const scalar forceScale,
            const scalar momentScale
        );

        //- Step of the free degrees of freedom: K dq = R (dF, dM)
        void solveStep
        (
            const scalarSquareMatrix& K,
            const vector& RF,
            const vector& RM,
            vector& dF,
            vector& dM
        ) const;

        //- Newmark acceleration and velocity of the current state
        void newmarkRates();


public:

    // Constructors

        //- Construct from name and dictionary; defaultFree is the value
        //- of freeTranslations when the dictionary gives neither
        //- 'freeTranslations' nor 'directions'
        rigidStructure
        (
            const word& name,
            const dictionary& dict,
            const vector& defaultFree = Zero
        );


    //- Destructor
    ~rigidStructure() = default;


    // Member Functions

        // Access

            const word& name() const noexcept { return name_; }

            //- True if any degree of freedom is solved
            bool free() const;

            //- True for the dynamic motion solver
            bool dynamicSolver() const noexcept
            {
                return motionSolver_ == "dynamic";
            }

            //- True while the dynamic solver is active (after
            //- dynamicStartTime)
            bool dynamic() const noexcept
            {
                return dynamicActive_;
            }

            //- Velocity of the reference point [m/s] and angular velocity
            //- [rad/s] (dynamic)
            const vector& velocity() const noexcept { return v_; }
            const vector& angularVelocity() const noexcept { return omega_; }

            //- Free components (0/1)
            const vector& freeTranslations() const noexcept
            {
                return freeTranslations_;
            }
            const vector& freeRotations() const noexcept
            {
                return freeRotations_;
            }

            //- Density [kg/m3]
            scalar density() const noexcept { return density_; }

            //- Mass [kg] for the volume V reported by the owner
            scalar mass(const scalar V) const;

            //- Reference point [m]
            const point& refPoint() const noexcept { return refPoint_; }

            //- Set the reference point (owners that know the geometry;
            //- ignored when the dictionary gave one)
            void setRefPoint(const point& x);

            //- Centre of mass [m]; xV = the volume centroid of the owner
            point centreOfMass(const point& xV) const;

            //- Translation of the reference point [m]
            const vector& displacement() const noexcept { return d_; }

            //- Rotation vector [rad]
            const vector& rotation() const noexcept { return theta_; }

            //- Soil change below which a quasi-static update is made
            scalar gate() const noexcept { return gate_; }

            //- Secant stiffnesses of the last update
            const vector& secantTranslation() const noexcept { return kT_; }
            const vector& secantRotation() const noexcept { return kR_; }

            //- Number of quasi-static updates so far
            label nUpdates() const noexcept { return nUpdates_; }


        // Kinematics

            //- Displacement of the point x [m]
            vector displacementAt(const point& x) const
            {
                return d_ + (theta_ ^ (x - refPoint_));
            }

            //- Displacement of the points [m]
            tmp<vectorField> displacementAt(const pointField& x) const;

            //- Set the state
            void setState(const vector& d, const vector& theta);

            //- Prescribed state at time t (no free degrees of freedom)
            void setPrescribed(const scalar t);


        // Loads

            //- External force [N] and moment [N m] about the reference
            //- point at time t (the moment of the force and the applied
            //- loads included)
            vector externalForce(const scalar t) const;
            vector externalMoment(const scalar t) const;

            //- Set the loads applied by an owner [N], [N m] (about the
            //- reference point)
            void setAppliedLoads(const vector& F, const vector& M)
            {
                appliedForce_ = F;
                appliedMoment_ = M;
            }

            const vector& appliedForce() const noexcept
            {
                return appliedForce_;
            }
            const vector& appliedMoment() const noexcept
            {
                return appliedMoment_;
            }

            //- Weight and buoyancy [N] and their moment [N m] about the
            //- reference point for the submerged volume V with centroid
            //- xV, water density rhoW, gravity g (z up)
            void gravityLoads
            (
                const scalar V,
                const point& xV,
                const scalar rhoW,
                const scalar g,
                vector& F,
                vector& M
            ) const;


        // Balance of the free degrees of freedom

            //- Update the free components from the static residual force
            //- F [N] and moment M [N m] (both about the reference point)
            //- with the mass m [kg]: quasi-static or dynamic (Newmark)
            //- according to the motion solver; returns the residual
            //- relative to the scales
            scalar correct
            (
                const vector& F,
                const vector& M,
                const scalar m,
                const scalar forceScale,
                const scalar momentScale
            );

            //- Dynamic update (Newmark, implicit in the secant)
            scalar correctDynamic
            (
                const vector& F,
                const vector& M,
                const scalar m,
                const scalar forceScale,
                const scalar momentScale
            );

            //- Update the free components from the residual force F [N]
            //- and moment M [N m] (both about the reference point);
            //- returns the residual relative to the scales (the largest
            //- of |F_i|/forceScale and |M_i|/momentScale over the free
            //- components)
            scalar correctQuasiStatic
            (
                const vector& F,
                const vector& M,
                const scalar forceScale,
                const scalar momentScale
            );

            //- Start a new time step of size dt ending at time t (the
            //- secant is kept, the last residual forgotten; dynamic: the
            //- state of the step start is stored and the state predicted;
            //- the dynamics are activated once t > dynamicStartTime)
            void newStep(const scalar dt = 0, const scalar t = GREAT);

            //- The owner updates the structure once per solve, from the
            //- converged reaction (the secant then spans the solves)
            bool outerUpdate() const noexcept
            {
                return update_ == "outer";
            }


        // Output

            //- Report the state
            void report() const;

            //- Append the state and loads to
            //- postProcessing/structures/<name>.dat (master only; the
            //- header is written when the file is created)
            void writeHistory(const fileName& caseDir, const scalar t) const;

            //- Write the state to uniform/structures/<name> of the time
            //- directory (restart)
            void writeState(const Time& runTime) const;

            //- Read the state written at the current time when present
            //- (restart); returns whether it was read
            bool readState(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
