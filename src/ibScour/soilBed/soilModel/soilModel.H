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
    Foam::soilModel

Description
    Base class of the models of the soil region under a height-field bed:
    what is solved in the soil (pore pressure only, Biot u-p, ...), selected
    at run time from the 'type' entry of the soil dictionary, as a
    turbulence model is.

    The soil region is a column mesh (bedRegionMesh::NewFvMesh) whose top
    patch 'soilTop' is the bed (face i is bed face i) and whose bottom
    patch is 'soilBottom'. Pressures are in Pa, excess over the hydrostatic
    pressure of the still water level, with the same datum as the pressure
    given on the top patch.

Usage
    \verbatim
    soil
    {
        type            darcyStorage;   // soilModel
        waterDensity    1000;           // [kg/m3]
        gravity         9.81;           // [m/s2]
        ...                             // model coefficients
    }
    \endverbatim

SourceFiles
    soilModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_soilModel_H
#define Foam_soilModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class poreGenerationModel;
class rigidStructure;

/*---------------------------------------------------------------------------*\
                          Class soilModel Declaration
\*---------------------------------------------------------------------------*/

class soilModel
{
protected:

    // Protected Data

        //- The soil column mesh
        const fvMesh& mesh_;

        //- Settings
        dictionary dict_;

        //- Index of the top (bed) patch
        const label topPatchi_;

        //- Index of the bottom patch
        const label bottomPatchi_;

        //- Cells inside a structure (rigid, impermeable; empty: none)
        labelList maskedCells_;

        //- Per-cell structure mask (1 inside a structure)
        scalarField cellMask_;

        //- True if any processor has structure cells (the masked face
        //- mobility exchanges boundary values: all ranks must take part)
        bool anyMasked_;

        //- Rigid-body displacement of the structure cells: translation
        //- [m] of the reference point, rotation vector [rad] and the
        //- reference point [m]
        vector structureDisplacement_;
        vector structureRotation_;
        point structureRefPoint_;

        //- Displacement of the masked cells [m] (maskedCells_ order)
        vectorField maskedDisplacement_;

        //- The structure the mask belongs to (not owned; nullptr: none)
        rigidStructure* structure_;

        //- Residual pore-pressure generation model (nullptr
        //- without a poreGeneration sub-dictionary)
        autoPtr<poreGenerationModel> generation_;

        //- Vertical effective stress [Pa], depth below the bed [m] and
        //- excess pore pressure p - p_bed [Pa] of the cells at the last
        //- solve, set by the soil region
        scalarField sigmaV0_;
        scalarField depth_;
        scalarField pExcess_;

        //- Water density [kg/m3]
        const scalar rhoW_;

        //- Gravity [m/s2]
        const scalar g_;


public:

    //- Name of the top (bed) patch of the soil mesh
    static const word topPatchName;

    //- Name of the bottom patch of the soil mesh
    static const word bottomPatchName;

    //- Runtime type information
    TypeName("soilModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            soilModel,
            dictionary,
            (const fvMesh& mesh, const dictionary& dict),
            (mesh, dict)
        );


    // Constructors

        //- Construct from the soil mesh and dictionary
        soilModel(const fvMesh& mesh, const dictionary& dict);


    // Selectors

        //- Return the selected soil model
        static autoPtr<soilModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~soilModel() = default;


    // Member Functions

        //- The soil mesh
        const fvMesh& mesh() const noexcept { return mesh_; }

        //- Index of the top (bed) patch
        label topPatchi() const noexcept { return topPatchi_; }

        //- Index of the bottom patch
        label bottomPatchi() const noexcept { return bottomPatchi_; }

        //- Set the structure mask (1 for cells inside a rigid, impermeable
        //- structure): the displacement is fixed to zero there and the
        //- pore-fluid flux across their faces vanishes
        virtual void setStructureMask(const scalarField& mask);

        //- Cells inside a structure
        const labelList& maskedCells() const noexcept { return maskedCells_; }

        //- Structure mask per cell (empty when none)
        const scalarField& cellMask() const noexcept { return cellMask_; }

        //- Prescribe the rigid-body translation [m] of the structure
        //- cells (no rotation; models that solve the displacement)
        void setStructureDisplacement(const vector& d)
        {
            setStructureMotion(d, Zero, structureRefPoint_);
        }

        //- Prescribe the rigid-body motion of the structure cells:
        //- D(x) = d + theta ^ (x - xRef)
        void setStructureMotion
        (
            const vector& d,
            const vector& theta,
            const point& xRef
        );

        //- Translation of the structure reference point [m]
        const vector& structureDisplacement() const noexcept
        {
            return structureDisplacement_;
        }

        //- Rotation vector of the structure [rad]
        const vector& structureRotation() const noexcept
        {
            return structureRotation_;
        }

        //- Reference point of the structure motion [m]
        const point& structureRefPoint() const noexcept
        {
            return structureRefPoint_;
        }

        //- Displacement of the masked cells [m] (maskedCells() order)
        const vectorField& maskedDisplacement() const noexcept
        {
            return maskedDisplacement_;
        }

        //- Attach the structure the mask belongs to (not owned)
        void attachStructure(rigidStructure* s) noexcept
        {
            structure_ = s;
        }

        //- The attached structure (nullptr: none)
        rigidStructure* structure() const noexcept
        {
            return structure_;
        }

        //- Effective force [N] of the soil on the structure cells (sum of
        //- the effective traction over the soil-structure faces; models
        //- that solve the displacement, else zero)
        virtual vector structureForce() const
        {
            return Zero;
        }

        //- Set the vertical effective stress [Pa] and the depth below
        //- the bed [m] of the cells (for the generation model)
        void setColumnState
        (
            const scalarField& sigmaV0,
            const scalarField& depth,
            const scalarField& pExcess
        );

        //- The generation model (nullptr without one)
        const poreGenerationModel* generation() const noexcept
        {
            return generation_.get();
        }

        //- Undrained generation rate f [Pa/s] per cell over nSteps flow
        //- steps ending now, from the current excess pore pressure
        //- p - pTop (pTop per column, cells layer-major); empty field
        //- without a generation model; accumulates the generated pressure
        tmp<scalarField> generationRate
        (
            const label nSteps,
            const scalarField& pTop
        );

        //- Remap the generation state (called by the models' remap)
        void remapGeneration
        (
            const labelList& c0,
            const labelList& c1,
            const scalarField& w
        );

        //- Write the generation state fields (if any)
        void writeGeneration() const;

        //- Restart: true when the state fields of the model were read
        //- from the time directory
        virtual bool restored() const { return false; }

        //- Write the state that is not a registered field (restart)
        virtual void writeState() const {}

        //- Face mobility: the uniform value, zero on the faces of masked
        //- cells
        tmp<surfaceScalarField> mobilityFace
        (
            const dimensionedScalar& mobility
        ) const;

        //- Unit weight of water [N/m3]
        scalar gammaW() const noexcept { return rhoW_*g_; }

        //- Pore pressure [Pa]
        virtual const volScalarField& p() const = 0;

        //- Advance the soil over nSteps flow time steps (the old-time
        //- fields are those of the last call) with the pore pressure pTop
        //- [Pa] on the faces of the top patch
        virtual void correct(const scalarField& pTop, const label nSteps) = 0;

        //- Hydraulic diffusivity [m2/s] (uniform models)
        virtual scalar diffusivity() const = 0;

        //- Set the pore pressure of every cell (initial equilibrium with
        //- the bed pressure); default: nothing
        virtual void initialise(const scalarField& pCells)
        {}

        //- Remap the fields after the columns moved: new cell i takes
        //- (1 - w[i]) field[c0[i]] + w[i] field[c1[i]] of the old values
        virtual void remap
        (
            const labelList& c0,
            const labelList& c1,
            const scalarField& w
        )
        {}

        //- True if the model solves the skeleton displacement (effective
        //- stress available)
        virtual bool solvesDisplacement() const
        {
            return false;
        }

        //- Effective stress increment [Pa] (tension positive); only with
        //- solvesDisplacement()
        virtual const volSymmTensorField& sigma() const;

        //- Set the effective traction [Pa] on the faces of the top patch
        //- (e.g. the bed shear stress); default: ignored
        virtual void topTraction(const vectorField& t)
        {}

        //- Prescribe the displacement [m] on the top faces with value
        //- fraction 1 (0: traction faces); tangential: only the
        //- tangential components (the normal one keeps its traction);
        //- default: ignored
        virtual void topDisplacement
        (
            const scalarField& fraction,
            const vectorField& D,
            const bool tangential
        )
        {}

        //- Effective traction [Pa] (tension positive) exerted on the soil
        //- by the faces of a patch, from the total effective stress (the
        //- stress of the adjacent cells on faces with a prescribed
        //- displacement, whose one-sided face strain is first order);
        //- default: empty
        virtual tmp<vectorField> patchEffectiveTraction
        (
            const label patchi
        ) const
        {
            return tmp<vectorField>::New();
        }

        //- Effective traction on the top faces
        tmp<vectorField> topEffectiveTraction() const
        {
            return patchEffectiveTraction(topPatchi_);
        }

        //- Set the geostatic effective stress from the vertical effective
        //- stress sigmaV0 [Pa] (compression positive) per cell and the
        //- lateral earth-pressure coefficient k0; default: ignored
        virtual void geostatic(const scalarField& sigmaV0, const scalar k0)
        {}

        //- Helper: remap a field as in remap()
        template<class Type>
        static tmp<Field<Type>> remapField
        (
            const Field<Type>& old,
            const labelList& c0,
            const labelList& c1,
            const scalarField& w
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "soilModelTemplates.C"
#endif

#endif

// ************************************************************************* //
