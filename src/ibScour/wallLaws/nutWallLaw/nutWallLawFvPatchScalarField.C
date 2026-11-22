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

#include "nutWallLawFvPatchScalarField.H"
#include "wallLawFaces.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutWallLawFvPatchScalarField::calcNut() const
{
    const label patchi = patch().index();

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const scalarField nuw(turbModel.nu(patchi));

    scalarField uTau;
    scalarField y;
    scalarField magUt;
    wallLawFaces::evaluate(*law_, Uw, nuw, uTau, y, magUt);

    auto tnutw = tmp<scalarField>::New(patch().size(), Zero);
    scalarField& nutw = tnutw.ref();

    forAll(nutw, facei)
    {
        if (magUt[facei] > VSMALL)
        {
            nutw[facei] =
                max(sqr(uTau[facei])*y[facei]/magUt[facei] - nuw[facei], 0);
        }
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutWallLawFvPatchScalarField::nutWallLawFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF),
    lawDict_(),
    law_()
{}


Foam::nutWallLawFvPatchScalarField::nutWallLawFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict),
    lawDict_(dict.subDict("wallLaw")),
    law_(ibWallLaw::New(lawDict_))
{}


Foam::nutWallLawFvPatchScalarField::nutWallLawFvPatchScalarField
(
    const nutWallLawFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper),
    lawDict_(ptf.lawDict_),
    law_(ibWallLaw::New(lawDict_))
{}


Foam::nutWallLawFvPatchScalarField::nutWallLawFvPatchScalarField
(
    const nutWallLawFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(wfpsf, iF),
    lawDict_(wfpsf.lawDict_),
    law_(ibWallLaw::New(lawDict_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutWallLawFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const scalarField nuw(turbModel.nu(patchi));

    scalarField uTau;
    scalarField y;
    scalarField magUt;
    wallLawFaces::evaluate(*law_, Uw, nuw, uTau, y, magUt);

    return y*uTau/nuw;
}


void Foam::nutWallLawFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.beginBlock("wallLaw");
    lawDict_.write(os, false);
    os.endBlock();
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nutWallLawFvPatchScalarField
    );
}


// ************************************************************************* //
