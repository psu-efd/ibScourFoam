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

#include "wallLawTurbulence.H"
#include "wallLawFaces.H"
#include "fvMatrices.H"
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(wallLawTurbulence, 0);
    addToRunTimeSelectionTable(option, wallLawTurbulence, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::wallLawTurbulence::wallValues
(
    labelList& cells,
    scalarField& omega,
    scalarField& kEq,
    scalarField& fRough,
    labelList& nbrCells
) const
{
    const auto& turbModel = mesh_.lookupObject<turbulenceModel>
    (
        turbulenceModel::propertiesName
    );
    const volVectorField& U = turbModel.U();

    // Area-weighted face sums per wall cell
    Map<label> cellIndex;
    DynamicList<label> wallCells;
    DynamicList<scalar> sumA;
    DynamicList<scalar> sumOmegaA;
    DynamicList<scalar> sumKA;
    DynamicList<scalar> sumFA;
    DynamicList<label> nbr;

    label nClamped = 0;
    scalar minRatio = GREAT;

    const cellList& meshCells = mesh_.cells();
    const faceList& meshFaces = mesh_.faces();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    const volScalarField* kPtr = mesh_.cfindObject<volScalarField>("k");
    const scalar Cmu25 = Foam::pow(Cmu_, 0.25);

    for (const label patchi : patchIDs_)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const labelUList& faceCells = Up.patch().faceCells();
        const scalarField& magSf = Up.patch().magSf();
        const scalarField nuw(turbModel.nu(patchi));

        scalarField uTau;
        scalarField y;
        scalarField magUt;
        wallLawFaces::evaluate
        (
            *law_, Up, nuw, uTau, y, magUt, &nClamped, &minRatio
        );

        forAll(faceCells, facei)
        {
            scalar om = 0;
            scalar nutLog = 0;
            law_->wallTurbulence
            (
                uTau[facei], y[facei], nuw[facei], Cmu_, beta1_, Aplus_,
                om, nutLog
            );

            const label celli = faceCells[facei];
            auto iter = cellIndex.find(celli);
            label i = -1;
            if (iter.good())
            {
                i = iter.val();
            }
            else
            {
                i = wallCells.size();
                cellIndex.insert(celli, i);
                wallCells.append(celli);
                sumA.append(0);
                sumOmegaA.append(0);
                sumKA.append(0);
                sumFA.append(0);

                // Wall-normal neighbour: across the face opposite the wall
                // face (hexahedra; internal faces only)
                const label wallFace = Up.patch().start() + facei;
                const label oppFace =
                    meshCells[celli].opposingFaceLabel(wallFace, meshFaces);
                label n = -1;
                if (oppFace >= 0 && mesh_.isInternalFace(oppFace))
                {
                    n = (own[oppFace] == celli ? nei[oppFace] : own[oppFace]);
                }
                nbr.append(n);
            }

            // Smooth (ks+ <= 2.25): equilibrium k; rough (ks+ >= 5): dk/dn = 0.
            // ks+ from u_tau, or (roughnessScale k) from u_k = Cmu^1/4 sqrt(k)
            // of the wall-normal neighbour, which stays finite where u_tau
            // vanishes (separation, reattachment)
            scalar uScale = uTau[facei];
            if (roughnessFromK_ && kPtr)
            {
                const label nb = (nbr[i] >= 0 ? nbr[i] : celli);
                uScale = Cmu25*Foam::sqrt(max((*kPtr)[nb], scalar(0)));
            }
            const scalar ksPlus = law_->ks()*uScale/nuw[facei];
            const scalar f =
                min
                (
                    max
                    (
                        (ksPlus - ksPlusSmooth_)
                       /max(ksPlusRough_ - ksPlusSmooth_, SMALL),
                        scalar(0)
                    ),
                    scalar(1)
                );

            const scalar a = magSf[facei];
            sumA[i] += a;
            sumOmegaA[i] += a*om;
            sumKA[i] += a*max(nutLog*om, SMALL);
            sumFA[i] += a*f;
        }
    }

    // Distorted wall cells: report once per time step if any, and at
    // write times always
    const Time& runTime = mesh_.time();
    if (reportIndex_ != runTime.timeIndex())
    {
        reportIndex_ = runTime.timeIndex();
        reduce(nClamped, sumOp<label>());
        reduce(minRatio, minOp<scalar>());

        if (nClamped > 0 || runTime.writeTime())
        {
            Info<< "wallLawTurbulence: " << nClamped
                << " wall faces with distorted cells (n.d < "
                << wallLawFaces::yMinFraction << " |d|, wall distance "
                << "bounded); min n.d/|d| = " << minRatio << endl;
        }
    }

    cells.transfer(wallCells);
    omega = scalarField(sumOmegaA)/scalarField(sumA);
    kEq = scalarField(sumKA)/scalarField(sumA);
    fRough = scalarField(sumFA)/scalarField(sumA);
    nbrCells.transfer(nbr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::wallLawTurbulence::wallLawTurbulence
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(name, modelType, dict, mesh),
    patchIDs_(),
    law_(),
    Cmu_(0.09),
    beta1_(0.075),
    Aplus_(26),
    roughnessFromK_(false),
    ksPlusSmooth_(2.25),
    ksPlusRough_(5),
    reportIndex_(-1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::wallLawTurbulence::constrain
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const word& fieldName = eqn.psi().name();

    labelList cells;
    scalarField omega;
    scalarField kEq;
    scalarField f;
    labelList nbrCells;
    wallValues(cells, omega, kEq, f, nbrCells);

    if (fieldName == "omega")
    {
        eqn.setValues(cells, omega);
    }
    else if (fieldName == "k")
    {
        // Smooth: equilibrium; rough: the wall-normal neighbour (dk/dn = 0)
        const scalarField& k = eqn.psi().primitiveField();
        scalarField values(cells.size());

        forAll(cells, i)
        {
            const scalar kN =
                (nbrCells[i] >= 0 ? max(k[nbrCells[i]], SMALL) : kEq[i]);
            values[i] = (1 - f[i])*kEq[i] + f[i]*kN;
        }

        eqn.setValues(cells, values);
    }
    else
    {
        FatalErrorInFunction
            << "Field " << fieldName << " is not supported; use k or omega"
            << exit(FatalError);
    }
}


bool Foam::fv::wallLawTurbulence::read(const dictionary& dict)
{
    if (!fv::option::read(dict))
    {
        return false;
    }

    if (!coeffs_.readIfPresent("fields", fieldNames_))
    {
        fieldNames_.clear();

        for (const word f : {"k", "omega"})
        {
            if (mesh_.foundObject<volScalarField>(f))
            {
                fieldNames_.push_back(f);
            }
        }
    }

    patchIDs_ = mesh_.boundaryMesh().indices
    (
        coeffs_.get<wordRes>("patches")
    );

    law_ = ibWallLaw::New(coeffs_.subDict("wallLaw"));
    Cmu_ = coeffs_.getOrDefault<scalar>("Cmu", 0.09);
    beta1_ = coeffs_.getOrDefault<scalar>("beta1", 0.075);
    Aplus_ = coeffs_.getOrDefault<scalar>("Aplus", 26);

    const word scale(coeffs_.getOrDefault<word>("roughnessScale", "uTau"));
    if (scale != "uTau" && scale != "k")
    {
        FatalIOErrorInFunction(coeffs_)
            << "Unknown roughnessScale " << scale << "; use uTau or k"
            << exit(FatalIOError);
    }
    roughnessFromK_ = (scale == "k");
    ksPlusSmooth_ = coeffs_.getOrDefault<scalar>("ksPlusSmooth", 2.25);
    ksPlusRough_ = coeffs_.getOrDefault<scalar>("ksPlusRough", 5);

    Info<< "    wallLawTurbulence on patches "
        << UIndirectList<word>(mesh_.boundaryMesh().names(), patchIDs_)
        << " applied to " << fieldNames_ << endl;

    fv::option::resetApplied();

    return true;
}


// ************************************************************************* //
