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

#include "darcyStorage.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(darcyStorage, 0);
    addToRunTimeSelectionTable(soilModel, darcyStorage, dictionary);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::darcyStorage::readBulkModulus(const dictionary& dict)
{
    if (dict.found("bulkModulus"))
    {
        return dict.get<scalar>("bulkModulus");
    }

    const scalar Sr = dict.getOrDefault<scalar>("saturation", 1);
    const scalar pAbs = dict.getOrDefault<scalar>("absolutePressure", 1.2e5);
    const scalar Kw = dict.getOrDefault<scalar>("waterBulkModulus", 2.2e9);

    return 1/(1/Kw + (1 - Sr)/pAbs);
}


Foam::wordList Foam::darcyStorage::patchTypes() const
{
    const fvBoundaryMesh& bm = mesh_.boundary();

    wordList types(bm.size(), zeroGradientFvPatchScalarField::typeName);

    forAll(bm, patchi)
    {
        if (bm[patchi].coupled())
        {
            types[patchi] = bm[patchi].type();
        }
    }

    types[topPatchi_] = fixedValueFvPatchScalarField::typeName;

    if (bottomPressure_)
    {
        types[bottomPatchi_] = fixedValueFvPatchScalarField::typeName;
    }

    return types;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::darcyStorage::darcyStorage(const fvMesh& mesh, const dictionary& dict)
:
    soilModel(mesh, dict),
    porosity_(dict.getOrDefault<scalar>("porosity", 0.4)),
    conductivity_(dict.get<scalar>("conductivity")),
    bulkModulus_(readBulkModulus(dict)),
    storage_
    (
        "storage",
        dimensionSet(-1, 1, 2, 0, 0),
        porosity_/bulkModulus_
    ),
    mobility_
    (
        "mobility",
        dimensionSet(-1, 3, 1, 0, 0),
        conductivity_/gammaW()
    ),
    bottomPressure_
    (
        Function1<scalar>::NewIfPresent
        (
            "bottomPressure", dict, word::null, &mesh
        )
    ),
    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPressure, dict.getOrDefault<scalar>("p0", 0)),
        patchTypes()
    ),
    restored_(false)
{
    Info<< "    Darcy storage soil: porosity " << porosity_
        << ", conductivity " << conductivity_ << " m/s, bulk modulus "
        << bulkModulus_ << " Pa, diffusivity " << diffusivity() << " m2/s"
        << (bottomPressure_ ? ", fixed bottom pressure" : "") << endl;

    // Restart: the pore pressure written at the restart time
    restored_ =
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ).typeHeaderOk<volScalarField>(true);
    if (restored_)
    {
        p_.oldTime() = p_;
        Info<< "    Darcy storage soil: pore pressure restored from time "
            << mesh.time().timeName() << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::darcyStorage::initialise(const scalarField& pCells)
{
    p_.primitiveFieldRef() = pCells;
    p_.correctBoundaryConditions();
    p_.oldTime() = p_;
}


void Foam::darcyStorage::remap
(
    const labelList& c0,
    const labelList& c1,
    const scalarField& w
)
{
    initialise(remapField(p_.primitiveField(), c0, c1, w));
    remapGeneration(c0, c1, w);
}


void Foam::darcyStorage::correct(const scalarField& pTop, const label nSteps)
{
    p_.boundaryFieldRef()[topPatchi_] == pTop;

    if (bottomPressure_)
    {
        p_.boundaryFieldRef()[bottomPatchi_] ==
            bottomPressure_->value(mesh_.time().value());
    }

    // Over nSteps flow time steps: the time derivative spans nSteps dt,
    // i.e. the mobility is scaled by nSteps in the dt-based equation
    fvScalarMatrix pEqn
    (
        fvm::ddt(storage_, p_)
      - fvm::laplacian(nSteps*mobilityFace(mobility_), p_)
    );

    // Residual pore-pressure generation: the undrained rate f
    // over the same nSteps interval, with the storage of the equation
    const tmp<scalarField> tf(generationRate(nSteps, pTop));
    if (tf().size())
    {
        volScalarField::Internal f
        (
            IOobject("poreGenerationRate", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar(dimPressure/dimTime, Zero)
        );
        f.field() = tf();
        pEqn -= nSteps*storage_*f;
    }

    pEqn.solve();
}


// ************************************************************************* //
