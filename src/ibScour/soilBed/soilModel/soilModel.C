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

#include "soilModel.H"
#include "poreGenerationModel.H"
#include "fvc.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(soilModel, 0);
    defineRunTimeSelectionTable(soilModel, dictionary);
}

const Foam::word Foam::soilModel::topPatchName("soilTop");
const Foam::word Foam::soilModel::bottomPatchName("soilBottom");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::soilModel::soilModel(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    dict_(dict),
    topPatchi_(mesh.boundaryMesh().findPatchID(topPatchName)),
    bottomPatchi_(mesh.boundaryMesh().findPatchID(bottomPatchName)),
    anyMasked_(false),
    structureDisplacement_(Zero),
    structureRotation_(Zero),
    structureRefPoint_(Zero),
    maskedDisplacement_(),
    structure_(nullptr),
    generation_(poreGenerationModel::New(mesh, dict)),
    sigmaV0_(),
    depth_(),
    pExcess_(),
    rhoW_(dict.getOrDefault<scalar>("waterDensity", 1000)),
    g_(dict.getOrDefault<scalar>("gravity", 9.81))
{
    if (topPatchi_ < 0 || bottomPatchi_ < 0)
    {
        FatalErrorInFunction
            << "Soil mesh " << mesh.name() << " needs the patches "
            << topPatchName << " and " << bottomPatchName
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::soilModel> Foam::soilModel::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "    Selecting soil model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "soilModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<soilModel>(ctorPtr(mesh, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volSymmTensorField& Foam::soilModel::sigma() const
{
    FatalErrorInFunction
        << "Soil model " << type() << " does not solve the displacement;"
        << " no effective stress" << exit(FatalError);
    return NullObjectRef<volSymmTensorField>();
}


void Foam::soilModel::setStructureMask(const scalarField& mask)
{
    cellMask_ = mask;

    DynamicList<label> cells(mask.size()/8);
    forAll(mask, celli)
    {
        if (mask[celli] > 0.5)
        {
            cells.append(celli);
        }
    }
    maskedCells_ = labelList(std::move(cells));
    maskedDisplacement_.setSize(maskedCells_.size());
    setStructureMotion
    (
        structureDisplacement_,
        structureRotation_,
        structureRefPoint_
    );

    const label nMasked = returnReduce(maskedCells_.size(), sumOp<label>());
    anyMasked_ = (nMasked > 0);

    Info<< "    Soil model: " << nMasked << " cells inside structures" << endl;
}


void Foam::soilModel::setStructureMotion
(
    const vector& d,
    const vector& theta,
    const point& xRef
)
{
    structureDisplacement_ = d;
    structureRotation_ = theta;
    structureRefPoint_ = xRef;

    const pointField& C = mesh_.cellCentres();
    maskedDisplacement_.setSize(maskedCells_.size());
    forAll(maskedCells_, i)
    {
        maskedDisplacement_[i] = d + (theta ^ (C[maskedCells_[i]] - xRef));
    }
}


void Foam::soilModel::setColumnState
(
    const scalarField& sigmaV0,
    const scalarField& depth,
    const scalarField& pExcess
)
{
    sigmaV0_ = sigmaV0;
    depth_ = depth;
    pExcess_ = pExcess;
}


Foam::tmp<Foam::scalarField> Foam::soilModel::generationRate
(
    const label nSteps,
    const scalarField& pTop
)
{
    if (!generation_ || sigmaV0_.size() != mesh_.nCells())
    {
        return tmp<scalarField>::New();
    }

    const scalar t = mesh_.time().value();
    const scalar dt = nSteps*mesh_.time().deltaTValue();

    // Current excess pore pressure (not the one of the last diagnostics)
    const scalarField& pCells = p().primitiveField();
    const label nCol = pTop.size();
    scalarField pExcess(pCells.size());
    forAll(pExcess, celli)
    {
        pExcess[celli] = pCells[celli] - pTop[celli % nCol];
    }

    return generation_->rate(t, dt, sigmaV0_, depth_, pExcess);
}


void Foam::soilModel::remapGeneration
(
    const labelList& c0,
    const labelList& c1,
    const scalarField& w
)
{
    if (generation_)
    {
        generation_->remap(c0, c1, w);
    }
}


void Foam::soilModel::writeGeneration() const
{
    if (generation_)
    {
        generation_->write();
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::soilModel::mobilityFace
(
    const dimensionedScalar& mobility
) const
{
    auto tmob = tmp<surfaceScalarField>::New
    (
        IOobject
        (
            "mobilityFace",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        mobility
    );

    if (cellMask_.size() == mesh_.nCells() && anyMasked_)
    {
        volScalarField mask
        (
            IOobject("structureMask", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar(dimless, Zero),
            zeroGradientFvPatchField<scalar>::typeName
        );
        mask.primitiveFieldRef() = cellMask_;
        mask.correctBoundaryConditions();

        // The interpolate is 0.5 or 1 when one or both cells are masked
        tmob.ref() *= pos0(0.25 - fvc::interpolate(mask));
    }

    return tmob;
}


// ************************************************************************* //
