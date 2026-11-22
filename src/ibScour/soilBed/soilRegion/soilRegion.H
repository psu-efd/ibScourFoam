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
    Foam::soilRegion

Description
    Soil region under a height-field bed: a column mesh built in memory
    from the floor patch of the flow mesh (bedRegionMesh::NewFvMesh, one
    column of nLayers cells per floor face, top patch 'soilTop' at the bed
    elevation, bottom patch 'soilBottom' at a fixed elevation), a soil
    model selected at run time (soilModel) and the effective-stress
    bookkeeping of the momentary liquefaction criterion.

    Column geometry: point j of level k lies at
    z = zBed_j - eta_k (zBed_j - zBottom), with eta_k the normalised depth
    of the level below the bed, graded geometrically so that the cells are
    thinner near the bed ('grading' = thickness of the bottom cell over
    that of the top cell). The columns are repositioned (setBed) whenever
    the bed elevation at the floor-patch points changes.

    Top boundary of the soil: the dynamic pressure of the flow on the bed
    face [Pa] plus an optional progressive surface wave
    p0 cos(k x.e - omega t) ('surfaceWave': the flow solver has no free
    surface; the wave loading of a seabed under a current or an
    oscillatory flow is prescribed), plus an optional prescribed unloading
    p_u(t) (Function1 of
    time; p_u = -rho g (h(t) - h0) for a water level falling under a rigid
    lid). With 'pressureRelaxation' the bed pressure is filtered in time
    (first-order relaxation), removing the step-to-step jumps of the
    pressure of immersed-boundary cells at reclassifications.

    Diagnostics per column: the upward hydraulic gradient at the bed
    i = (1/gamma_w) dp/dz (positive for upward seepage), evaluated from the
    top cell or, with 'gradientDepth', as the pressure difference over that
    depth below the bed (the gradient relevant to the mobility of the top
    grain layers; a sudden change of the bed pressure gives a singular
    surface gradient), and the depth
    of the layer in which the excess pore pressure p - p_top reaches the
    initial vertical effective stress gamma' (zBed - z) (momentary
    liquefaction), with gamma' = (s - 1)(1 - n) gamma_w.

Usage
    \verbatim
    soil
    {
        nLayers         40;
        bottom          -2;             // elevation of the soil bottom [m]
        grading         20;             // bottom/top cell thickness
        specificGravity 2.65;           // of the grains (gamma')
        k0              0.5;            // lateral earth-pressure coefficient
        porosity        0.4;
        writeMesh       true;           // constant/soil/polyMesh for ParaView
        equilibrateInitial true;        // p = bed pressure at the first step
        gradientDepth   0.02;           // i over this depth [m] (0: top cell)
        pressureRelaxation 0.1;         // time filter of the bed pressure [s]
        pressureRampTime 5;             // flow pressure x min(t/T, 1) [s]
                                        // (start-up transients; default 0)
        interval        1;              // soil solve every N flow steps
        unloading       constant 0;     // p_u(t) [Pa], Function1 (optional)
        surfaceWave                     // optional progressive wave
        {
            amplitude   563;            // p0 [Pa]
            wavelength  2.89;           // [m]
            period      1.6;            // [s]
            direction   (1 0 0);        // of propagation (default x)
        }
        type            darcyStorage;   // soilModel and its coefficients
        ...
    }
    \endverbatim

SourceFiles
    soilRegion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_soilRegion_H
#define Foam_soilRegion_H

#include "fvMesh.H"
#include "volFields.H"
#include "soilModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class soilRegion Declaration
\*---------------------------------------------------------------------------*/

class soilRegion
{
    // Private Data

        //- The flow mesh
        const fvMesh& flowMesh_;

        //- The floor patch of the flow mesh
        const polyPatch& floorPatch_;

        //- Settings
        dictionary dict_;

        //- Number of cell layers
        const label nLayers_;

        //- Elevation of the soil bottom [m]
        const scalar zBottom_;

        //- Normalised depth of the point levels below the bed (size
        //- nLayers + 1, level 0 at the bottom)
        scalarList eta_;

        //- Submerged unit weight of the soil gamma' [N/m3]
        scalar gammaS_;

        //- Lateral earth-pressure coefficient k0 of the geostatic state
        const scalar k0_;

        //- Overburden [Pa] added to the geostatic vertical effective
        //- stress (soil above the modelled column; default 0)
        const scalar overburden_;

        //- The soil mesh
        autoPtr<fvMesh> meshPtr_;

        //- The soil model
        autoPtr<soilModel> modelPtr_;

        //- Prescribed unloading p_u(t) [Pa] added to the top pressure
        autoPtr<Function1<scalar>> unloading_;

        //- Ramp time [s] of the flow pressure passed to the soil (the
        //- start-up transient of a flow would otherwise be locked into a
        //- low-permeability soil; 0: none)
        scalar pressureRampTime_;

        //- Prescribed progressive surface wave on the top pressure:
        //- amplitude [Pa], wavenumber [1/m], angular frequency [1/s] and
        //- unit direction of propagation
        scalar waveAmp_;
        scalar waveK_;
        scalar waveOmega_;
        vector waveDir_;

        //- Bed elevation at the floor-patch points (patch-local numbering)
        scalarField zBedPoints_;

        //- Pore pressure on the top faces at the last step [Pa]
        scalarField pTop_;

        //- Loss of vertical effective stress [Pa]: p - p_top (rigid
        //- skeleton) or sigma'_zz of the Biot models
        volScalarField pExcess_;

        //- Initial vertical effective stress gamma' (zBed - z) [Pa]
        volScalarField sigmaV0_;

        //- Depth of the cell centres below the bed [m]
        scalarField depth_;

        //- Whether the mesh has been written
        mutable bool meshWritten_;

        //- Equilibrate the soil with the bed pressure at the first step
        //- (true unless the soil dictionary gives p0)
        bool equilibrate_;

        //- Depth below the bed over which the upward gradient is evaluated
        //- [m] (0: from the top cell)
        const scalar gradientDepth_;

        //- Relaxation time of the bed pressure seen by the soil [s]
        //- (0: none)
        const scalar pressureRelaxation_;

        //- The soil is advanced every 'interval' flow steps (the models
        //- scale their mobility; the Biot model is costly)
        const label interval_;

        //- Flow steps since the last soil solve
        label stepsSinceSolve_;

        //- The geostatic stress has been passed to the model
        bool geostaticSet_;

        //- Relaxed bed pressure [Pa]
        scalarField pBedRelaxed_;

        //- The relaxed bed pressure holds a value (set from the first flow
        //- pressure, or restored): a restart does not relax from zero
        bool pBedRelaxedSet_;

        //- Effective traction [Pa] on the bed faces passed to the model
        //- (the bed shear stress), empty if not set
        vectorField bedTraction_;

        //- Prescribed displacement of the bed faces: value fraction and
        //- displacement [m] (empty: none)
        scalarField bedFraction_;
        vectorField bedDisplacement_;
        bool bedTangential_ = false;


    // Private Member Functions

        //- Normalised depths of the point levels
        static scalarList levelDepths
        (
            const label nLayers,
            const scalar grading
        );

        //- Point levels of unit-height columns above zBottom
        static scalarList initialLevels
        (
            const scalarList& eta,
            const scalar zBottom
        );

        //- Elevation of the top face of column facei
        scalar zTop(const label facei) const;

        //- Update the diagnostic fields
        void updateDiagnostics();

        //- Move the column points under the bed elevation zBed at the
        //- floor-patch points (no remap of the fields)
        void positionColumns(const scalarField& zBedPoints);

        //- Restart: position the columns on the bed written at the
        //- time directory (uniform/<region>/zBedPoints) when present
        void readBedState();

        //- Write the bed elevation the columns stand under (restart)
        void writeBedState() const;


public:

    // Constructors

        //- Construct from the flow mesh, the floor patch and the soil
        //- dictionary
        soilRegion
        (
            const fvMesh& flowMesh,
            const word& floorPatch,
            const dictionary& dict
        );


    // Member Functions

        //- The soil mesh
        const fvMesh& mesh() const noexcept { return *meshPtr_; }

        //- The soil model
        const soilModel& model() const noexcept { return *modelPtr_; }

        //- The floor patch
        const polyPatch& floorPatch() const noexcept { return floorPatch_; }

        //- Number of columns (floor faces) on this processor
        label nColumns() const noexcept { return floorPatch_.size(); }

        //- Number of cell layers
        label nLayers() const noexcept { return nLayers_; }

        //- Submerged unit weight of the soil [N/m3]
        scalar gammaS() const noexcept { return gammaS_; }

        //- Lateral earth-pressure coefficient of the geostatic state
        scalar k0() const noexcept { return k0_; }

        //- Overburden added to the geostatic vertical effective stress [Pa]
        scalar overburden() const noexcept { return overburden_; }

        //- Elevation of the soil bottom [m]
        scalar bottom() const noexcept { return zBottom_; }

        //- Bed elevation (top face) of column facei [m]
        scalar bedElevation(const label facei) const { return zTop(facei); }

        //- Position the columns under the bed elevation zBed at the
        //- floor-patch points (patch-local numbering)
        void setBed(const scalarField& zBedPoints);

        //- Advance one flow time step with the dynamic bed pressure
        //- pBed [Pa] on the floor faces (the unloading is added)
        void evolve(const scalarField& pBed);

        //- Pore pressure on the top faces at the last step [Pa]
        const scalarField& pTop() const noexcept { return pTop_; }

        //- Set the structure mask of the soil cells (1 inside a rigid,
        //- impermeable structure); recomputed by the caller after setBed
        void setStructureMask(const scalarField& mask);

        //- Prescribe the rigid-body displacement [m] of the structure
        //- cells for the next step (models that solve the displacement)
        void setStructureDisplacement(const vector& d)
        {
            modelPtr_->setStructureDisplacement(d);
        }

        //- Effective force [N] of the soil on the structure cells
        vector structureForce() const
        {
            return modelPtr_->structureForce();
        }

        //- Attach the rigid structure the mask belongs to (not owned):
        //- its balance is then solved by the soil model
        void attachStructure(rigidStructure* s)
        {
            modelPtr_->attachStructure(s);
        }

        //- Structure mask per soil cell (empty when none)
        const scalarField& structureMask() const
        {
            return modelPtr_->cellMask();
        }

        //- Set the effective traction [Pa] on the bed faces for the next
        //- step (models that solve the displacement)
        void setBedTraction(const vectorField& t);

        //- Prescribe the displacement [m] of the bed faces with value
        //- fraction 1 (0: traction faces) for the next step (models that
        //- solve the displacement)
        void setBedDisplacement
        (
            const scalarField& fraction,
            const vectorField& D,
            const bool tangential = false
        );

        //- Effective traction [Pa] (tension positive) on the bed faces
        //- (models that solve the displacement; else empty)
        tmp<vectorField> bedTraction() const
        {
            return modelPtr_->topEffectiveTraction();
        }

        //- Upward hydraulic gradient at the bed, per floor face
        tmp<scalarField> topGradient() const;

        //- Pore pressure at the gradient-evaluation depth below the bed
        //- (the top cell if gradientDepth is 0), per floor face [Pa]
        tmp<scalarField> pressureAtDepth() const;

        //- Depth of the momentarily liquefied layer, per floor face [m]
        tmp<scalarField> liquefiedDepth() const;

        //- Capacity index per floor face [Pa m]: the integral over the
        //- column of the available vertical effective stress
        //- max(sigma'_v0 - loss, 0) (the Broms-type lateral resistance of
        //- a pile scales with it)
        tmp<scalarField> capacityIndex() const;

        //- Cell centre elevations of column facei, from the top
        tmp<scalarField> columnZ(const label facei) const;

        //- Loss of vertical effective stress of column facei, from the
        //- top [Pa] (the pExcess field)
        tmp<scalarField> columnExcess(const label facei) const;

        //- Pore pressure above the bed pressure, p - p_top, of column
        //- facei, from the top [Pa]
        tmp<scalarField> columnPorePressure(const label facei) const;

        //- Write the mesh (once, if writeMesh) and the diagnostic fields
        void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
