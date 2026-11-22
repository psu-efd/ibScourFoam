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

#include "deformingBed.H"
#include "bedRegionMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "fixedValuePointPatchFields.H"
#include "sedimentBedFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(deformingBed, 0);
    addToRunTimeSelectionTable(functionObject, deformingBed, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::deformingBed::calcShearStress()
{
    const volVectorField& U = lookupObject<volVectorField>(UName_);
    const fvPatchVectorField& Up = U.boundaryField()[patchID_];

    // Effective viscosity of the wall faces (nu_t of the wall condition)
    scalarField nuEff(Up.size(), nu_);
    const auto* nutPtr = findObject<volScalarField>("nut");
    if (nutPtr)
    {
        nuEff += nutPtr->boundaryField()[patchID_];
    }

    const vectorField n(mesh_.boundary()[patchID_].nf());
    const vectorField snGradU(Up.snGrad());

    // Bed faces of the bed mesh are the bed faces of the flow mesh, in order
    const polyPatch& bp =
        bedMeshPtr_->boundaryMesh()[sedimentPtr_->patchID()];
    const labelList& faceLabels = sedimentPtr_->aMesh().faceLabels();

    areaVectorField& tauField = sedimentPtr_->tau();
    vectorField& tau = tauField.primitiveFieldRef();

    forAll(faceLabels, i)
    {
        const label facei = faceLabels[i] - bp.start();
        tau[i] =
           -nuEff[facei]*(snGradU[facei] - n[facei]*(n[facei] & snGradU[facei]));
        tau[i].z() = 0;
    }
    tauField.correctBoundaryConditions();
}


void Foam::functionObjects::deformingBed::setBedDisplacement()
{
    // Bed elevation at the bed points. Point j of the bed patch of the bed
    // mesh is local point meshPoints[j] of the bed patch of the flow mesh.
    const polyPatch& bp =
        bedMeshPtr_->boundaryMesh()[sedimentPtr_->patchID()];
    const scalarField zBed(sedimentPtr_->pointElevation());

    const polyPatch& pp = mesh_.boundaryMesh()[patchID_];
    const labelList& meshPoints = pp.meshPoints();
    const pointField& points = mesh_.points();

    // Change of the bed points from their current position (vertical)
    scalarField dz(pp.nPoints(), Zero);
    forAll(bp.meshPoints(), j)
    {
        const label pointi = bp.meshPoints()[j];
        dz[pointi] = zBed[j] - points[meshPoints[pointi]].z();
    }

    // The displacement is relative to the initial points: add the change.
    // Vertical component solver (displacementComponentLaplacian, component
    // z: scalar pointDisplacementz) or vector solver (pointDisplacement).
    if (foundObject<pointScalarField>("pointDisplacementz"))
    {
        pointScalarField& pD =
            lookupObjectRef<pointScalarField>("pointDisplacementz");
        pointPatchScalarField& bcBase = pD.boundaryFieldRef()[patchID_];

        if (!isA<fixedValuePointPatchScalarField>(bcBase))
        {
            FatalErrorInFunction
                << "pointDisplacementz of patch " << patchName_
                << " must be fixedValue" << exit(FatalError);
        }
        auto& bc = refCast<fixedValuePointPatchScalarField>(bcBase);
        bc == scalarField(bc) + dz;
        return;
    }

    if (!foundObject<pointVectorField>("pointDisplacement"))
    {
        FatalErrorInFunction
            << "No pointDisplacement(z) field: the flow mesh must be a"
            << " dynamicMotionSolverFvMesh with a displacement motion solver"
            << " (displacementComponentLaplacian with component z, or"
            << " displacementLaplacian)" << exit(FatalError);
    }

    pointVectorField& pD =
        lookupObjectRef<pointVectorField>("pointDisplacement");
    pointPatchVectorField& bcBase = pD.boundaryFieldRef()[patchID_];

    if (!isA<fixedValuePointPatchVectorField>(bcBase))
    {
        FatalErrorInFunction
            << "pointDisplacement of patch " << patchName_
            << " must be fixedValue" << exit(FatalError);
    }
    auto& bc = refCast<fixedValuePointPatchVectorField>(bcBase);

    vectorField disp(bc);
    forAll(disp, pointi)
    {
        disp[pointi].z() += dz[pointi];
    }
    bc == disp;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::deformingBed::deformingBed
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    patchName_(dict.get<word>("bedPatch")),
    patchID_(mesh_.boundaryMesh().findPatchID(patchName_)),
    UName_("U"),
    CName_("C"),
    moveMesh_(true),
    nu_(0)
{
    if (patchID_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Bed patch " << patchName_ << " not found"
            << exit(FatalIOError);
    }

    read(dict);

    Info<< type() << " " << name << ":" << nl;

    bedMeshPtr_ = bedRegionMesh::New(mesh_, patchName_, "bedRegion");

    sedimentPtr_.reset
    (
        new bedSediment(mesh_, *bedMeshPtr_, dict, patchName_)
    );

    // Without a z_b field the bathymetry is that of the mesh: z_b = the
    // elevation of the bed-face centres
    if (!sedimentPtr_->zbRead())
    {
        const polyPatch& bp =
            bedMeshPtr_->boundaryMesh()[sedimentPtr_->patchID()];
        const labelList& faceLabels = sedimentPtr_->aMesh().faceLabels();
        const vectorField Cf(mesh_.boundaryMesh()[patchID_].faceCentres());

        scalarField zb(faceLabels.size());
        forAll(faceLabels, i)
        {
            zb[i] = Cf[faceLabels[i] - bp.start()].z();
        }
        sedimentPtr_->setElevation(zb);

        Info<< "    Initial z_b from the bed-face centres of the mesh: "
            << gMin(zb) << " to " << gMax(zb) << " m" << endl;
    }

    Info<< endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::deformingBed::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    UName_ = dict.getOrDefault<word>("U", "U");
    CName_ = dict.getOrDefault<word>("C", "C");
    moveMesh_ = dict.getOrDefault<bool>("moveMesh", true);

    const IOdictionary& transportProperties =
        mesh_.lookupObject<IOdictionary>("transportProperties");
    nu_ = dimensionedScalar("nu", dimViscosity, transportProperties).value();

    return true;
}


bool Foam::functionObjects::deformingBed::execute()
{
    calcShearStress();

    // Fixed strips entrain no suspended sediment: erodibility per bed
    // face to the sedimentBedFlux condition (for the next C solve)
    if (auto* CRef = getObjectPtr<volScalarField>(CName_))
    {
        auto* bedRef = dynamic_cast<sedimentBedFluxFvPatchScalarField*>
        (
            &CRef->boundaryFieldRef()[patchID_]
        );

        if (bedRef)
        {
            const polyPatch& bp =
                bedMeshPtr_->boundaryMesh()[sedimentPtr_->patchID()];
            const labelList& faceLabels =
                sedimentPtr_->aMesh().faceLabels();
            const scalarField& fixed =
                sedimentPtr_->fixed().primitiveField();

            scalarField erodible(bedRef->size(), scalar(1));
            forAll(faceLabels, i)
            {
                erodible[faceLabels[i] - bp.start()] = 1 - fixed[i];
            }
            bedRef->setErodible(erodible);
        }
    }

    // Suspended load: net deposition D - E of the bed patch condition
    // (sedimentBedFlux of libsuspendedSediment) to the Exner equation
    const auto* CPtr = findObject<volScalarField>(CName_);
    if (CPtr)
    {
        const auto* bedPtr =
            isA<sedimentBedFluxFvPatchScalarField>
            (
                CPtr->boundaryField()[patchID_]
            );

        if (bedPtr)
        {
            const scalarField netPatch(bedPtr->D() - bedPtr->E());

            // Bed faces of the bed mesh are the bed faces of the flow
            // mesh, in order (the face labels are those of the bed mesh)
            const polyPatch& bp =
                bedMeshPtr_->boundaryMesh()[sedimentPtr_->patchID()];
            const labelList& faceLabels =
                sedimentPtr_->aMesh().faceLabels();

            scalarField net(faceLabels.size(), Zero);
            forAll(faceLabels, i)
            {
                net[i] = netPatch[faceLabels[i] - bp.start()];
            }

            sedimentPtr_->setSuspendedExchange(net);
        }
    }

    // Flow spin-up: the bed is frozen before morphStartTime
    if (sedimentPtr_->morphActive())
    {
        sedimentPtr_->advance(time_.deltaTValue());
    }

    // After the advance: the surface then holds the shear stress and the
    // bed load of this step and the new z_b (as the z_b field files)
    if (time_.writeTime())
    {
        sedimentPtr_->writeSurface();
        sedimentPtr_->writeMorphologyState();
    }

    // Also during the spin-up, so that an initial z_b field (bathymetry)
    // is taken over by the mesh at the first step
    if (moveMesh_)
    {
        setBedDisplacement();
    }

    return true;
}


bool Foam::functionObjects::deformingBed::write()
{
    return true;
}


// ************************************************************************* //
