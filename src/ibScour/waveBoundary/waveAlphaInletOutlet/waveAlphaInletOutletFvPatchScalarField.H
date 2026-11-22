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
    Foam::waveAlphaInletOutletFvPatchScalarField

Group
    grpWaveBoundaryConditions

Description
    Phase-fraction condition of a wave paddle that imposes the wave model's
    phase fraction on inflow faces only and is zero-gradient on outflow
    faces, in the manner of inletOutlet.

    The stock waveAlpha condition fixes the phase fraction on every face
    of the paddle. Where the flow leaves through the paddle (the trough
    of the wave) a face whose imposed value is water next to a cell that
    the interior has already left dry carries a water-density mass flux
    out of an air cell, and the velocity of that cell grows without bound
    (three-dimensional basin runs, 2026-10-31: 2 to 32 m/s in 0.1 s in the
    first cell above the water line at the paddle). With this condition
    the outflow faces carry the interior phase fraction, as any outlet
    does, and the inflow faces the wave model's.

Usage
    \table
        Property     | Description                  | Required | Default
        waveDictName | wave dictionary name         | no       | waveProperties
        phi          | flux field name              | no       | phi
    \endtable

    Example of the condition specification:
    \verbatim
    inlet
    {
        type            waveAlphaInletOutlet;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    waveAlphaInletOutletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_waveAlphaInletOutletFvPatchScalarField_H
#define Foam_waveAlphaInletOutletFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class waveAlphaInletOutletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class waveAlphaInletOutletFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    typedef waveAlphaInletOutletFvPatchScalarField this_bctype;
    typedef mixedFvPatchScalarField parent_bctype;

    // Private Data

        //- Wave dictionary name
        word waveDictName_;

        //- Flux field name
        word phiName_;


public:

    //- Runtime type information
    TypeName("waveAlphaInletOutlet");


    // Constructors

        //- Construct from patch and internal field
        waveAlphaInletOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveAlphaInletOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        waveAlphaInletOutletFvPatchScalarField
        (
            const this_bctype&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        waveAlphaInletOutletFvPatchScalarField
        (
            const this_bctype&,
            const DimensionedField<scalar, volMesh>&
        );

        //- No copy without an internal field
        waveAlphaInletOutletFvPatchScalarField(const this_bctype&) = delete;

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
