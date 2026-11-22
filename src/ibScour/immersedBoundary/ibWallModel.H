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
    Foam::ibWallModel

Description
    Wall model of the immersed boundary: velocity and turbulence values of
    the IB cells and the wall shear stress.

    For each IB cell c (wall distance y_c, unit normal n, wall velocity
    u_w) the wall law is applied at its live neighbours j in the
    cell-point-cell stencil (parallel-consistent through the stencil's
    mapDistribute), each at its own wall distance y_j > y_c:

        - u_tau,j from the tangential relative velocity |u_t,j| at y_j;
        - u_tau = sum w_j u_tau,j / sum w_j, with weights
          w_j = (n . r_j/|r_j|)^2/|r_j|^2, r_j = x_j - x_c, which favour
          neighbours in the wall-normal direction;
        - U_c = u_w + u_tau u+(y_c+) t + <u_n,j y_c/y_j> n, with t the
          weighted tangential direction.

    No interpolation to an artificial sampling point is needed, so u_tau
    does not depend on a sampling distance and the wall distance of every
    velocity sample is exact.

    Turbulence quantities of the IB cell:

        - omega = sqrt(omega_vis^2 + omega_log^2), k = nu_t,log omega,
          nu_t,log = kappa u_tau y_c (1 - exp(-y+/A+))^2,
          epsilon = Cmu k omega;
        - with 'consistentNut true' (default) the eddy viscosity is chosen
          so that the stress transmitted to the live neighbours equals
          u_tau^2: across [y_c, y_j] the wall law requires the effective
          viscosity nu_sec,j = u_tau (y_j - y_c)/(u+(y_j+) - u+(y_c+)); with
          linear interpolation to the face, nu_t,c = 2(<nu_sec> - nu)
          - <nu_t,j> (bounded below by 0) (cf. Roman et al. 2009, Comput.
          Fluids 38:1510-1527); otherwise nu_t,c = nu_t,log.

    The kinematic wall shear stress is tau_w = u_tau^2 t.

    Two modes ('mode' entry):

        - velocity: the IB-cell velocity is imposed (U_c above);
        - stress (default): the IB cells solve their momentum, including the
          pressure gradient, like wall-function cells of a body-fitted mesh.
          The wall shear stress acts on the blocked faces as c (U - u_w) A_w,
          with A_w the projected wall area in the cell and
          c = u_tau^2/max(|u_t|, u_tau u+(y_c+)): the stress is u_tau^2 when
          the cell velocity has reached the wall-law value and smaller while
          it is still slower (e.g. a cell just uncovered by the bed), so the
          wall force is bounded by u_tau^2 A_w.
          u_tau is taken from the IB cell's own tangential velocity at y_c
          when the cell centre is far enough from the wall,
          y_c >= sampleBelow h (default 0.25). Closer to the wall (a moving
          bed may bring y_c arbitrarily close to zero, and a cell that has
          just been uncovered has no developed velocity yet) u_tau and the
          tangential direction are blended linearly, with weight
          1 - y_c/(sampleBelow h), towards the values sampled at the live
          neighbours further from the wall. In this range the cell's own
          u_tau is limited to the viscous scaling
          u_tau(sampleBelow h) sqrt(sampleBelow h/y_c): a smooth or laminar
          law follows it, while a rough log law has no solution once
          y_c < Cs ks/E and would grow without bound. A cell without live
          neighbours uses its own wall law at sampleBelow h alone (as in
          the velocity mode;
          cf. the y+-adaptive wall function of Xu & Liu 2021, Int. J. Numer.
          Meth. Fluids 93:1929-1946). sampleBelow 0 uses the cell value
          only (then a solidThreshold of about 0.2 is needed).
          Cells uncovered by the bed start from the wall-law velocity
          u_w + u_tau u+(y_c+) t reconstructed from the sampled values.
          With the same weight the eddy viscosity of the IB cell is blended
          from the turbulence model's value towards the wall-consistent
          value of the velocity mode (consistentNut): at small y+ the model
          value is nearly zero and the faces to the live neighbours would
          transmit only part of u_tau^2, which slows the cell and its
          neighbours (a low-shear streak).
          The k and omega values use max(y_c, yTurbMin h), yTurbMin = 0.05
          by default, so that omega_vis stays bounded for y_c -> 0. k and omega are fixed to the wall-function
          values (kTreatment equilibrium, default; 'transported' solves k
          but over-predicts the production in IB cells).

Usage
    \verbatim
    wallModel
    {
        mode            stress;     // stress (default) | velocity
        sampleBelow     0.25;       // optional, stress mode
        directionHeight 0.02;       // optional: direction of the reported
                                    // wall shear from the velocity at this
                                    // height above the wall [m] (0: cell)
        roughnessScale  uTau;       // optional: uTau | k (ks+ of the k blend)
        wallStrain      model;      // optional: model | wall (IB eddy viscosity)
        ksPlusSmooth    2.25;       // optional: equilibrium k below
        ksPlusRough     5;          // optional: dk/dn = 0 above
        yTurbMin        0.05;       // optional, stress mode
        reportMaxShear  false;      // optional, diagnostics
        consistentNut   true;       // optional, velocity mode
        wallLaw
        {
            type        spalding;
        }
    }
    \endverbatim

SourceFiles
    ibWallModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ibWallModel_H
#define Foam_ibWallModel_H

#include "ibWallLaw.H"
#include "volFields.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class immersedBoundary;

/*---------------------------------------------------------------------------*\
                         Class ibWallModel Declaration
\*---------------------------------------------------------------------------*/

class ibWallModel
{
    // Private Data

        //- The immersed boundary
        const immersedBoundary& ib_;

        //- Wall law
        autoPtr<ibWallLaw> law_;

        //- Stress mode (true) or velocity mode (false)
        bool stressMode_;

        //- Make the IB-cell eddy viscosity consistent with u_tau^2
        bool consistentNut_;

        //- Stress mode: below this wall distance (relative to the cell
        //- extent) u_tau is blended towards the neighbour-sampled value
        scalar sampleBelow_;

        //- Report the IB cell with the largest wall shear stress and the
        //- cell with the largest velocity every step (diagnostics)
        bool reportMaxShear_;

        //- Stress mode: minimum distance for the turbulence values,
        //- relative to the cell extent
        scalar yTurbMin_;

        //- Model constants
        scalar Cmu_;
        scalar beta1_;
        scalar Aplus_;

        //- ks+ of the smooth/rough blend of k from u_k = Cmu^1/4 sqrt(k) of
        //- the neighbours (roughnessScale k) instead of u_tau (uTau)
        bool roughnessFromK_;

        //- ks+ bounds of the smooth/rough blend of k
        scalar ksPlusSmooth_;
        scalar ksPlusRough_;

        //- Height above the wall [m] at which the velocity giving the
        //- direction of the reported wall shear stress is sampled (0: the
        //- IB cell's own velocity); the magnitude and the momentum
        //- coupling are unchanged
        scalar directionHeight_;

        //- Stress mode: strain of the IB wall cells from the wall-consistent
        //- gradient (wallStrain wall) instead of the cell gradient (model)
        bool wallStrain_;


        // Geometry per IB cell

            //- Wall distance of the IB cell centre [m]
            scalarField yIB_;

            //- Cell extent in the wall-normal direction [m]
            scalarField hIB_;

            //- Unit wall normal (into the fluid)
            vectorField nIB_;

            //- Nearest wall point [m]
            pointField wallPoints_;

            //- Live neighbours (compact stencil addressing)
            labelListList nbrIdx_;

            //- Stress mode: IB cells with fewer live neighbours than this
            //- (minLiveNeighbours, default 3) do not solve their momentum;
            //- they take the wall-law velocity of the velocity mode. A cell
            //- boxed in by solid cells (a one-cell pit of a staircase bed)
            //- with one open face has no restoring force once its velocity
            //- grows: 10 -> 22 -> 90 000 m/s in the wave basin, 2026-10-31
            label minLiveNeighbours_;

            //- Indices (into ibCells) of the cells pinned that way
            DynamicList<label> pinned_;

            //- Weights of the live neighbours
            List<scalarList> nbrWeights_;

            //- Wall distance of the live neighbours [m]
            List<scalarList> nbrY_;

            //- Direction sampling: compact index of the cell at
            //- directionHeight above the wall (-1: none) and its map
            labelList dirIdx_;
            autoPtr<mapDistribute> dirMapPtr_;


        // Results per IB cell

            //- Velocity [m/s]
            vectorField Uib_;

            //- Friction velocity [m/s]
            scalarField uTau_;

            //- y+ of the IB cell centre
            scalarField yPlus_;

            //- Turbulence quantities
            scalarField kIB_;
            scalarField omegaIB_;
            scalarField nutIB_;

            //- Stress mode: weight of nutIB_ against the turbulence
            //- model's eddy viscosity in correctNut
            scalarField nutWeight_;

            //- Kinematic wall shear stress [m2/s2]
            vectorField tauW_;

            //- Stress mode: wall coefficient u_tau^2/|u_t| [m/s]
            scalarField wallCoeff_;

            //- Wall velocity at the wall points [m/s]
            vectorField Uw_;


        // Output fields

            //- Wall shear stress in IB cells (kinematic) [m2/s2]
            volVectorField wallShearStress_;

            //- y+ in IB cells
            volScalarField yPlusField_;


    // Private Member Functions

        //- Laminar kinematic viscosity per IB cell [m2/s]
        tmp<scalarField> nuIB() const;

        //- Write the reportMaxShear diagnostics
        void reportMaxShear
        (
            const volVectorField& U,
            const UList<vector>& flatU,
            const vectorField& Uw,
            const scalarField& nu
        ) const;

        //- Eddy viscosity of IB cell i for which the faces to its live
        //- neighbours transmit u_tau^2 according to the wall law
        scalar consistentNut
        (
            const label i,
            const scalar uTau,
            const scalar nu,
            const UList<scalar>& flatNut
        ) const;

        //- Specific dissipation omega and mixing-length eddy viscosity of
        //- an IB cell at wall distance y. omega blends the near-wall
        //- solution 6 nu/(beta1 (y + y_s)^2) with the log-layer value
        //- u_tau/(Cmu^0.5 kappa (y + dy)); y_s = sqrt(6 nu/(beta1 omega_w))
        //- gives the wall value omega_w = S_r u_tau^2/nu of the wall law
        //- (Wilcox rough-wall condition, y_s = 0 for a smooth wall) and dy
        //- is the wall law's roughness shift of the mixing length
        void wallTurbulence
        (
            const scalar uTau,
            const scalar y,
            const scalar nu,
            scalar& omega,
            scalar& nutLog
        ) const;

        //- Turbulent kinetic energy of IB cell i: the equilibrium value kEq
        //- for a smooth wall and the live neighbours' k (zero normal
        //- gradient, Roulund et al. 2005) for transitional and rough walls,
        //- blended linearly over 2.25 < ks+ < 5
        scalar wallK
        (
            const label i,
            const scalar uTau,
            const scalar nu,
            const scalar kEq,
            const UList<scalar>& flatK
        ) const;

        //- Wall law at the live neighbours of IB cell i: weighted friction
        //- velocity, tangential direction and normal-velocity ratio.
        //  Returns false if the cell has no live neighbours.
        bool sampleNeighbours
        (
            const label i,
            const UList<vector>& flatU,
            const vector& Uw,
            const scalar nu,
            scalar& uTau,
            vector& tHat,
            scalar& unRatio
        ) const;


public:

    // Constructors

        //- Construct from the immersed boundary and the wallModel dictionary
        ibWallModel(const immersedBoundary& ib, const dictionary& dict);


    //- Destructor
    ~ibWallModel() = default;


    // Member Functions

        // Update

            //- Build the map of the direction-sampling cells

            void updateDirectionSampling();


            //- Geometry and neighbour weights (after classification)
            void updateGeometry();

            //- Wall-model values of the IB cells from the velocity field
            void correct(const volVectorField& U);

            //- Set the turbulent viscosity in IB and solid cells
            void correctNut(volScalarField& nut) const;

            //- Stress mode, wallStrain wall: eddy viscosity of the IB cells
            //- from the SST closure with the strain of the wall-consistent
            //- gradient, i.e. the blocked bed faces carry the wall velocity
            //- instead of the average with the solid cell
            void wallStrainNut(volScalarField& nut) const;

            //- Initialise cells that changed from solid to fluid: IB cells
            //- with live neighbours get the wall-law velocity reconstructed
            //- from the neighbour-sampled u_tau, the others the wall velocity
            void initialiseFreshCells
            (
                const labelList& freshCells,
                volVectorField& U
            ) const;


        // Access (one value per IB cell, ordered as ib.ibCells())

            //- True in stress mode
            bool stressMode() const noexcept
            {
                return stressMode_;
            }

            //- Wall shear stress of the IB cells (kinematic) [m2/s2]
            const volVectorField& wallShearStress() const noexcept
            {
                return wallShearStress_;
            }

            //- Stress mode: wall coefficient u_tau^2/|u_t| [m/s]
            const scalarField& wallCoeff() const noexcept
            {
                return wallCoeff_;
            }

            //- Wall velocity [m/s]
            const vectorField& wallVelocity() const noexcept
            {
                return Uw_;
            }

            //- Velocity
            //- Stress mode: IB cells pinned to the wall-law velocity
            const UList<label>& pinned() const noexcept
            {
                return pinned_;
            }

            const vectorField& Uib() const noexcept
            {
                return Uib_;
            }

            //- Friction velocity
            const scalarField& uTau() const noexcept
            {
                return uTau_;
            }

            //- The wall law
            const ibWallLaw& law() const noexcept
            {
                return *law_;
            }

            //- Kinematic viscosity of the IB cells [m2/s]
            tmp<scalarField> nuOfIBCells() const
            {
                return nuIB();
            }

            //- Turbulent kinetic energy
            const scalarField& kIB() const noexcept
            {
                return kIB_;
            }

            //- Specific dissipation rate
            const scalarField& omegaIB() const noexcept
            {
                return omegaIB_;
            }

            //- Turbulent viscosity
            const scalarField& nutIB() const noexcept
            {
                return nutIB_;
            }

            //- Dissipation rate
            tmp<scalarField> epsilonIB() const;

            //- Kinematic wall shear stress
            const vectorField& tauW() const noexcept
            {
                return tauW_;
            }

            //- Wall distance of the IB cell centres
            const scalarField& yIB() const noexcept
            {
                return yIB_;
            }

            //- Wall points
            const pointField& wallPoints() const noexcept
            {
                return wallPoints_;
            }

            //- Unit wall normals
            const vectorField& normals() const noexcept
            {
                return nIB_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
