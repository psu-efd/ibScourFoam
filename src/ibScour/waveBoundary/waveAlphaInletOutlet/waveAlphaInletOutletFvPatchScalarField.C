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

\*---------------------------------------------------------------------------*/

#include "waveAlphaInletOutletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "waveModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveAlphaInletOutletFvPatchScalarField::
waveAlphaInletOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF),
    waveDictName_(waveModel::dictName),
    phiName_("phi")
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1;
}


Foam::waveAlphaInletOutletFvPatchScalarField::
waveAlphaInletOutletFvPatchScalarField
(
    const this_bctype& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper),
    waveDictName_(ptf.waveDictName_),
    phiName_(ptf.phiName_)
{}


Foam::waveAlphaInletOutletFvPatchScalarField::
waveAlphaInletOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF),
    waveDictName_(dict.getOrDefault<word>("waveDictName", waveModel::dictName)),
    phiName_(dict.getOrDefault<word>("phi", "phi"))
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1;

    if (!this->readValueEntry(dict))
    {
        fvPatchField<scalar>::operator=(refValue());
    }
}


Foam::waveAlphaInletOutletFvPatchScalarField::
waveAlphaInletOutletFvPatchScalarField
(
    const this_bctype& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(ptf, iF),
    waveDictName_(ptf.waveDictName_),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveAlphaInletOutletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    tmp<waveModel> tmodel
    (
        waveModel::lookupOrCreate
        (
            patch().patch(),
            internalField().mesh(),
            waveDictName_
        )
    );

    waveModel& model = tmodel.constCast();

    model.correct(db().time().value());

    refValue() = model.alpha();

    // Fixed value on inflow faces, zero gradient on outflow faces
    const auto& phip =
        patch().lookupPatchField<surfaceScalarField>(phiName_);

    valueFraction() = 1.0 - pos0(phip);

    parent_bctype::updateCoeffs();
}


void Foam::waveAlphaInletOutletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("waveDictName", waveDictName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        waveAlphaInletOutletFvPatchScalarField
    );
}

// ************************************************************************* //
