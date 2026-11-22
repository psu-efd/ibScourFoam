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
    Foam::soilConstitutiveModel

Description
    Base class of the stress-strain laws of the soil skeleton used by the
    Biot soil models, selected at run time from the 'type' entry of the
    'constitutive' sub-dictionary of the soil dictionary (as a turbulence
    model is). A law provides the effective stress increment for a
    displacement gradient and the implicit stiffness 2 mu + lambda of the
    segregated displacement equation (the explicit correction is the
    difference between its stress and the implicit linear-elastic part).

Usage
    \verbatim
    constitutive
    {
        type            linearElastic;  // | mohrCoulomb (later)
        E               5e7;            // Young's modulus [Pa]
        nu              0.3;            // Poisson's ratio
    }
    \endverbatim

SourceFiles
    soilConstitutiveModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_soilConstitutiveModel_H
#define Foam_soilConstitutiveModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class soilConstitutiveModel Declaration
\*---------------------------------------------------------------------------*/

class soilConstitutiveModel
{
protected:

    // Protected Data

        //- The soil mesh
        const fvMesh& mesh_;

        //- Settings
        dictionary dict_;


public:

    //- Runtime type information
    TypeName("soilConstitutiveModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            soilConstitutiveModel,
            dictionary,
            (const fvMesh& mesh, const dictionary& dict),
            (mesh, dict)
        );


    // Constructors

        //- Construct from the soil mesh and dictionary
        soilConstitutiveModel(const fvMesh& mesh, const dictionary& dict);


    // Selectors

        //- Return the selected law from the 'constitutive' sub-dictionary
        //- of soilDict (linearElastic if absent)
        static autoPtr<soilConstitutiveModel> New
        (
            const fvMesh& mesh,
            const dictionary& soilDict
        );


    //- Destructor
    virtual ~soilConstitutiveModel() = default;


    // Member Functions

        //- Shear modulus mu [Pa]
        virtual const volScalarField& mu() const = 0;

        //- Lame's first parameter lambda [Pa]
        virtual const volScalarField& lambda() const = 0;

        //- Implicit stiffness 2 mu + lambda [Pa]
        virtual const volScalarField& impK() const = 0;

        //- Constrained (oedometer) modulus [Pa] (= impK for linear
        //- elasticity)
        virtual scalar constrainedModulus() const = 0;

        //- Current drained bulk modulus per cell [Pa] for the
        //- stabilisation of the segregated Biot scheme (the fixed-stress
        //- split is stable with the bulk modulus, Kim, Tchelepi & Juanes
        //- 2011; the constrained modulus under-estimates the compliance
        //- unless the kinematics are one-dimensional, and a law whose
        //- skeleton softens must report the softer modulus); default:
        //- lambda + 2 mu/3 of the law's fields
        virtual tmp<volScalarField> stabilisationModulusField() const;

        //- Effective stress increment sigma [Pa] (tension positive) for
        //- the displacement gradient gradD, given the geostatic effective
        //- stress sigma0 (total effective stress = sigma0 + sigma); a
        //- plastic law updates its working plastic strain here
        virtual void correct
        (
            const volTensorField& gradD,
            const volSymmTensorField& sigma0,
            volSymmTensorField& sigma
        ) = 0;

        //- Commit the state (plastic strain) at the end of a converged
        //- time step; default: nothing
        virtual void updateState()
        {}

        //- Start of a time step, before the first corrector (laws that
        //- freeze a decision over the correctors); default: nothing
        virtual void beginStep()
        {}

        //- Forget the decisions frozen within the step (the Newton-Krylov
        //- solve calls it before every residual evaluation, so that the
        //- stress is a function of the strain of that evaluation alone
        //- and not of a finite-difference perturbation seen first);
        //- default: nothing
        virtual void unfreeze()
        {}

        //- The geostatic effective stress has been set (laws with a
        //- stress-dependent state or moduli); default: nothing
        virtual void geostatic(const volSymmTensorField& sigma0)
        {}

        //- Remap the state after the soil columns moved (see
        //- soilModel::remap); default: nothing
        virtual void remap
        (
            const labelList& c0,
            const labelList& c1,
            const scalarField& w
        )
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
