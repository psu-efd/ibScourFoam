/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Xiaofeng Liu, Penn State University
-------------------------------------------------------------------------------
License
    This file is part of ibScourFoam, built on OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::sedimentBedFluxFvPatchScalarField

Group
    grpWallBoundaryConditions

Description
    Erosion-deposition condition for the suspended-sediment volume
    concentration C at a bed patch: the Robin family of Cheng (1984), with
    the near-wall treatments of Liu (2014) (wall functions, HRN) and
    Liu (2018) (wall-resolved mesh, LRN) as members.

    The net upward sediment flux through the bed (into the fluid) is the
    entrainment minus the deposition,
    \f[
        \Gamma_f \, \mathrm{snGrad}(C) - w_n C_f = E - \beta C_f ,
    \f]
    where snGrad is the outward normal gradient at the face (so that
    \f$\Gamma_f\,\mathrm{snGrad}(C)\f$ is the diffusive flux into the
    fluid), \f$\Gamma_f\f$ the patch value of the diffusivity of the C
    equation (field DC of suspendedSedimentFoam), \f$w_n = w_s\cdot n\f$ the
    settling velocity through the patch (n outward) and \f$C_f\f$ the patch
    value. The left side is what the finite-volume equation carries
    through the face (diffusive plus settling flux), so the bed flux of the
    discrete equations is exactly \f$E - \beta C_f\f$. Written as an
    OpenFOAM mixed condition with
    \f$K = \beta - w_n \ge 0\f$:
    \f[
        f = \frac{K}{K + \Gamma_f \delta_c}, \quad
        \mathrm{refValue} = 0, \quad
        \mathrm{refGradient} = \frac{E}{\Gamma_f},
    \f]
    (\f$\delta_c\f$ = patch deltaCoeffs) which gives
    \f$C_f = (E + \Gamma_f\delta_c C_P)/(K + \Gamma_f\delta_c)\f$. For
    \f$\beta = w_n\f$ (deposition = settling flux) f = 0 and the condition is
    the fixed gradient \f$E/\Gamma_f\f$ of Liu (2014, 2018).

    The flux balance is Cheng (1984), J. Geophys. Res. 89(C5) 8209-8214,
    Eq. 3, Q_0 = -[k dc/dz + w c]_(z=0) = -beta (c - c_*) = E - D with
    E = beta c_* and D = beta c_0 (z upward, w the fall velocity): with the
    outward normal of the bed pointing down, Gamma_f snGrad(C) = -k dc/dz
    and -w_n C_f = -w c, i.e. the left side above is Cheng's Q_0.

    Modes (E = entrainment flux [m/s], C_b* = reference concentration at
    delta from the referenceConcentrationModel of sedimentProperties):
    \table
        mode           | E                                     | beta | where the patch is
        referenceLevel | w_n C_b*                              | w_n  | at y = delta (Liu 2014 Eq. 10); at the wall: the "naive" condition
        HRN            | w_n C_b* (1 + P2) or exp(P2), P2 = w_n delta/eps_s | w_n | wall, wall functions (Liu 2014 Eqs. 17-18)
        LRN            | w_n C_b* exp(int_0^delta w_n/eps_s dy) | w_n  | wall, y+ ~ 1 (Liu 2018 Eq. 9)
        Robin          | beta C_b*                             | beta | reference level (Cheng 1984; GITT benchmarks)
        nearestCell    | w_n C_b*                              | w_n  | wall; D = w_n C_P of the wall cell (Liu and Garcia 2008)
    \endtable
    Mode \c nearestCell is the exchange of Liu and Garcia (2008, J. Waterway
    Port Coastal Ocean Eng. 134(4) 203-217, Eq. 15) and of the immersed
    bed (ibSuspendedLoad): E = w_n C_b* at the reference level and
    D = w_n C_P with C_P the concentration of the wall cell ("the
    concentration at the nearest cell center"). The face carries no
    flux (zero gradient, no settling through the face: suspendedSedimentFoam
    sets phiC = phi there) and the solver adds E - D as a source of the
    wall cell, deposition implicit (CEqn.H, as ibSuspendedLoad does). This
    is the condition to use at a wall patch whose first-cell Peclet number
    w_n y1/Gamma_f is large (fine wall-function mesh, coarse sand): the
    other modes evaluate D from the face value, which the two-point flux
    then drives far above C_P and C_P below zero (wall-jet scour runs).

    Entry \c diffusivity (all modes) sets the diffusivity Gamma_f of the
    bed face, i.e. the eps_s of the theory between the wall and the first
    cell centre: \c firstCell (default for HRN) replaces the patch value of DC by the value of
    the wall-adjacent cell, the constant eps_s across the near-wall layer
    assumed by Liu (2014) (and implied by the Peclet numbers of his
    Table 3); \c patch (default for the other modes) keeps the patch value
    of DC, which contains the wall-function nut_w (HRN) or nu only (LRN).
    suspendedSedimentFoam applies it to DC before assembling the C equation
    (faceDiffusivity()). The HRN P2 = w_n delta/eps_s uses the same eps_s.
    With HRN and \c patch the first-cell Peclet number w_n y1/Gamma_f can
    exceed 1 (Liu 2014 case 2) and C in the wall cell then becomes negative.

    Findings of the verification runs: the HRN result is not mesh
    independent inside the guideline 0.2 k_s < y1 < 0.1 H of Liu (2014)
    when the Rouse number is of order 0.5 or more (C(y) scales with
    y1^Z); the LRN condition needs a near-wall mesh with w_s Delta/nu of
    about 0.1 or less (the equilibrium profile decays over nu/w_s).

    \c logLaw is an extension, not part of Liu (2014): it
    makes the two-point wall flux exact for the equilibrium profile in a
    log layer,
    eps_s = kappa u* y/sigma_c: Gamma_f = w_n y1/(1 - (delta/y1)^Z') with
    Z' = sigma_c w_n/(kappa u*), y1 = 1/deltaCoeffs, and E = w_n C_b*; the bed
    face value then plays the role of C at delta, and C1 = C_b* (delta/y1)^Z'
    at equilibrium independently of y1 (mesh study in a turbulent flume: the spread over
    y1 drops from 7.9 % to 1.4 % for Z = 0.06 and from 77 % to 19 % for
    Z = 0.6). Needs y1 > delta (else firstCell); entries kappa (0.41) and
    sigmaC (sedimentProperties).

    LRN, entry \c integral: \c computed (default) integrates w_n/DC from the
    wall to delta along the cell column above each face (hexahedral or
    prismatic cells, trapezoidal rule through the cell centres, the wall value
    being Gamma_f); \c table uses kappa Z I1(ks+, delta+) = (w_n/u*) I1,
    I1 = int_0^delta+ dy+/(1 + nut+) (Liu 2018 Eq. 10, nu + nut without
    sigma_c, as printed there), ks+ = ks u* / nu, from one of
    the tables selected by \c I1Model:
    - \c vanDriestCC (default): van Driest mixing length with the
      roughness shift of Cebeci & Chang (1978) Eqs. 8-9; ks+ <= 4.535 is
      the smooth van Driest table, which it reproduces exactly;
    - \c SST (option): computed with the k-omega SST model of OpenFOAM in
      1D columns, smooth (omegaWallFunction) and rough (omega_w = u*^2 S_R
      / nu), ks+ = 0-2000, delta+ = 0-1000 (consistent with the nut the C
      equation sees on an LRN mesh). Caveats: its smooth row misses the
      smooth van Driest table by 29 % (the omegaWallFunction near-wall
      state of the low-Re SST model is case dependent: nu_t+(15) = 0.78
      in the table column against 2.5 in flume columns with the same wall-unit grading); the Wilcox
      (2006) S_R constants of the rough rows were taken from the SU2
      source and still await a check against Wilcox (2006);
    - \c file: an interpolation2DTable given in the sub-dictionary I1Table.
    The built-in tables (I1Tables.H) are a recomputation by the authors, NOT the
    rough-wall table of Liu et al. (2013, IAHR Congress), which is not
    available; values outside the tables are clamped.

    With HRN and \c firstCell the wall-cell diffusivity of a log layer is
    Gamma_1 = kappa u* y1/sigma_c, so the first-cell Peclet number
    w_n y1/Gamma_1 equals Z' = sigma_c w_n/(kappa u*) for every y1; for
    Z' >= 1 the two-point wall flux drives C negative (van Rijn 1981 case, Z' = 1.19). The
    condition therefore stops with a FatalError when Z' >= 1 on any bed
    face (faces with u* = 0 are skipped) and names \c diffusivity logLaw
    as the remedy.

    The shear velocity is u* = sqrt(nuEff_f |snGrad(U)_t|) on the patch,
    i.e. the wall shear stress of the flow solution (wall function or
    resolved); it sets the HRN P2, the logLaw Z', the LRN wall units and
    the check above. The reference concentration uses the shear velocity
    of its Shields number: with \c shieldsNumber \c skinFriction (default
    of zysermanFredsoe and vanRijn1984, sedimentProperties) the grain
    shear velocity u*' from u*, the depth-mean tangential velocity V and
    the depth D of the cell column above the face (bed-normal walk through
    hexahedral or prismatic cells to the top boundary of the mesh; the
    column must not cross a processor boundary); with \c total the u*
    itself.

Usage
    \table
        Property     | Description                              | Required | Default
        mode         | referenceLevel, HRN, LRN, Robin or nearestCell | yes |
        diffusivity  | patch, firstCell or logLaw (bed-face eps_s) | no | firstCell (HRN), patch (others)
        form         | HRN: linear or exponential               | no  | linear
        integral     | LRN: computed or table                   | no  | computed
        I1Model      | LRN table: vanDriestCC, SST or file      | no  | vanDriestCC (file if I1Table given)
        ks           | LRN table: roughness height for ks+ [m]  | no  | 0
        I1Table      | LRN table: interpolation2DTable dict     | file |
        beta         | Robin: deposition velocity [m/s]         | Robin |
        DC           | name of the diffusivity field            | no  | DC
        ws           | name of the settling velocity (uniform)  | no  | ws
        U            | name of the velocity field               | no  | U
    \endtable
    \verbatim
    bottom
    {
        type        sedimentBedFlux;
        mode        HRN;
        diffusivity firstCell;      // default for HRN
        value       uniform 0;
    }
    \endverbatim

    Needs the registered fields DC and ws and the dictionary
    sedimentProperties of suspendedSedimentFoam.

SourceFiles
    sedimentBedFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef sedimentBedFluxFvPatchScalarField_H
#define sedimentBedFluxFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "referenceConcentrationModel.H"
#include "interpolation2DTable.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sedimentBedFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class sedimentBedFluxFvPatchScalarField
:
    public mixedFvPatchScalarField
{
public:

    //- Near-bed treatments
    enum modeType { referenceLevel, HRN, LRN, Robin, nearestCell };

    static const Enum<modeType> modeTypeNames;

    //- Bed-face diffusivity
    enum diffusivityType { patchValue, firstCell, logLaw };

    static const Enum<diffusivityType> diffusivityTypeNames;


private:

    // Private Data

        //- Near-bed treatment
        modeType mode_;

        //- eps_s at the bed face
        diffusivityType diffusivity_;

        //- von Karman constant (logLaw)
        scalar kappa_;

        //- HRN: exponential form (Eq. 17) instead of linear (Eq. 18)
        bool exponential_;

        //- LRN: integral from the table instead of the computed nut
        bool useTable_;

        //- LRN table: vanDriestCC (built in, default), SST (built in) or
        //  file (I1Table dictionary)
        word I1Model_;

        //- LRN table I1(ks+, delta+)
        autoPtr<interpolation2DTable<scalar>> I1Table_;

        //- LRN table: roughness height for ks+ [m]
        scalar ks_;

        //- Robin: deposition velocity [m/s]
        scalar beta_;

        //- Field names
        word DCName_, wsName_, UName_;

        //- Reference concentration model (from sedimentProperties)
        mutable autoPtr<suspendedSediment::referenceConcentrationModel> refModel_;

        //- Entrainment flux of the last update [m/s]
        scalarField E_;

        //- Deposition velocity beta of the last update [m/s]
        scalarField beta_f_;

        //- Reference concentration and shear velocity of the last update
        scalarField Cstar_, uStar_;

        //- Erodibility per face (1 erodible, 0 fixed: no pickup); empty:
        //- all erodible (set by deformingBed from its fixed strips)
        scalarField erodible_;

        //- Shear velocity of the Shields number of the last update (u*'
        //  with shieldsNumber skinFriction, else u*)
        scalarField uStarSkin_;

        //- LRN computed: cells above each face and their wall distances
        mutable List<labelList> columnCells_;
        mutable List<scalarList> columnY_;

        //- Skin friction: cells of the full column above each face (to the
        //  boundary at the top), their heights and the column depth
        mutable List<labelList> fullColumnCells_;
        mutable List<scalarList> fullColumnDy_;
        mutable scalarField columnDepth_;


    // Private Member Functions

        //- sedimentProperties (registered by the solver, else read)
        const dictionary& sedimentDict() const;

        //- Shear velocity from the wall shear stress of the flow solution
        tmp<scalarField> calcUStar() const;

        //- Model from sedimentProperties (read on first use)
        const suspendedSediment::referenceConcentrationModel& refModel() const;

        //- Depth-mean tangential velocity V and depth D of the cell column
        //  above each face (skin-friction Shields number)
        void depthMean(scalarField& V, scalarField& D) const;

        //- HRN firstCell: stop if the first-cell Peclet number
        //  sigma_c w_n/(kappa u*) >= 1 on any bed face
        void checkFirstCellPeclet(const scalarField& wn) const;

        //- Integral of w_n/DC from the wall to delta, per face (LRN computed)
        tmp<scalarField> wallIntegral
        (
            const scalarField& wn,
            const scalarField& Gamma,
            const scalar delta
        ) const;


public:

    //- Runtime type information
    TypeName("sedimentBedFlux");


    // Constructors

        //- Construct from patch and internal field
        sedimentBedFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sedimentBedFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        sedimentBedFluxFvPatchScalarField
        (
            const sedimentBedFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct setting internal field reference
        sedimentBedFluxFvPatchScalarField
        (
            const sedimentBedFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member Functions

        //- Bed-face diffusivity the solver must use for the C equation
        //  (diffusivity firstCell or logLaw); invalid tmp for patch
        tmp<scalarField> faceDiffusivity() const;

        //- Erodibility per face: C_b* and E are zero where it is 0
        void setErodible(const scalarField& erodible)
        {
            erodible_ = erodible;
        }

        //- Entrainment flux E [m/s] (volume per area and time)
        const scalarField& E() const
        {
            return E_;
        }

        //- Deposition flux D = beta C_f [m/s] (nearestCell: w_n C_P)
        tmp<scalarField> D() const;

        //- Exchange as a source of the wall cell (mode nearestCell): the
        //  face carries no flux; the solver adds E - w_n C_P
        bool cellExchange() const
        {
            return mode_ == modeType::nearestCell;
        }

        //- Deposition velocity beta of the last update [m/s] (w_n unless
        //  Robin)
        const scalarField& betaF() const
        {
            return beta_f_;
        }

        //- Reference concentration of the last update
        const scalarField& Cstar() const
        {
            return Cstar_;
        }

        //- Shear velocity of the last update [m/s]
        const scalarField& uStar() const
        {
            return uStar_;
        }

        //- Shear velocity of the Shields number of the last update [m/s]
        //  (skin friction u*' or the total u*)
        const scalarField& uStarSkin() const
        {
            return uStarSkin_;
        }

        //- Update the coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
