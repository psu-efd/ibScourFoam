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
    Foam::fv::wallLawTurbulence

Group
    grpFvOptionsConstraints

Description
    Turbulence of the wall cells of body-fitted walls from the wall law of
    the immersed boundary (ibWallLaw), as in the IB cells of the immersed
    boundary (ibWallModel), so that a body-fitted wall and an immersed
    wall differ only in their geometric representation.

    For every cell next to a face of the given patches, u_tau of each wall
    face follows from |U_t| at the normal distance y of the cell centre
    (wallLawFaces); with several wall faces the face values are
    area-weighted. Then

        - omega = sqrt(omega_vis^2 + omega_log^2) at y, with the roughness
          shift of the mixing length and the offset that gives
          omega_w = S_r u_tau^2/nu at the wall (ibWallLaw::wallTurbulence);
        - k: smooth walls (ks+ <= 2.25) the equilibrium value
          nu_t,log omega; transitional and rough walls (ks+ >= 5) dk/dn = 0
          (Roulund et al. 2005): the value of the wall-normal neighbour
          (the cell across the face opposite the wall face), as the live
          neighbours of an IB cell; in between (1 - f) k_eq + f k_n. A
          wall cell without such a neighbour (non-hexahedral, or across a
          processor boundary) takes k_eq.

    k is always fixed in the wall cells: solving it there with the
    production of the cell-centred velocity gradient and the fixed omega
    lets k (and nu_t) run away near rough walls.

    The wall distance is bounded below by 0.1 of the centre-to-face
    distance on distorted wall cells (wallLawFaces). Once per time step
    the number of such faces and the smallest n.d/|d| are reported when
    any face is bounded, and at write times always.

    Use it with nutWallLaw (same wall law) on the same patches.

Usage
    In constant/fvOptions:
    \verbatim
    wallTurbulence
    {
        type            wallLawTurbulence;
        patches         (bed pier);
        fields          (k omega);  // optional; default: those present
        wallLaw
        {
            type        vanDriest;
            ks          0.00065;
        }
        Cmu             0.09;       // optional
        beta1           0.075;      // optional
        Aplus           26;         // optional
        roughnessScale  uTau;       // optional: uTau | k, velocity scale
                                    // of ks+ in the smooth/rough k blend
        ksPlusSmooth    2.25;       // optional: equilibrium k below
        ksPlusRough     5;          // optional: dk/dn = 0 above
    }
    \endverbatim

See also
    Foam::nutWallLawFvPatchScalarField
    Foam::ibWallModel

SourceFiles
    wallLawTurbulence.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fv_wallLawTurbulence_H
#define Foam_fv_wallLawTurbulence_H

#include "fvOption.H"
#include "ibWallLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class wallLawTurbulence Declaration
\*---------------------------------------------------------------------------*/

class wallLawTurbulence
:
    public fv::option
{
    // Private Data

        //- Wall patches
        labelList patchIDs_;

        //- The wall law
        autoPtr<ibWallLaw> law_;

        //- Model constants
        scalar Cmu_;
        scalar beta1_;
        scalar Aplus_;

        //- ks+ of the smooth/rough blend of k from u_k = Cmu^1/4 sqrt(k)
        //- (roughnessScale k) instead of u_tau (uTau, default)
        bool roughnessFromK_;

        //- ks+ bounds of the smooth/rough blend of k
        scalar ksPlusSmooth_;
        scalar ksPlusRough_;

        //- Time index of the last distorted-cell report
        mutable label reportIndex_;


    // Private Member Functions

        //- Wall cells, their omega, equilibrium k, rough-wall fraction f
        //- and wall-normal neighbour cell (-1 if none)
        void wallValues
        (
            labelList& cells,
            scalarField& omega,
            scalarField& kEq,
            scalarField& fRough,
            labelList& nbrCells
        ) const;


public:

    //- Runtime type information
    TypeName("wallLawTurbulence");


    // Constructors

        //- Construct from components
        wallLawTurbulence
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        wallLawTurbulence(const wallLawTurbulence&) = delete;

        //- No copy assignment
        void operator=(const wallLawTurbulence&) = delete;


    //- Destructor
    virtual ~wallLawTurbulence() = default;


    // Member Functions

        //- Fix the wall-cell values of k or omega
        virtual void constrain(fvMatrix<scalar>& eqn, const label fieldi);

        //- Read the settings
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
