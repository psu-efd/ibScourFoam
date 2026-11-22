/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Xiaofeng Liu, Penn State University
-------------------------------------------------------------------------------
License
    This file is part of ibScourFoam, built on OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "sedimentBedExchange.H"
#include "sedimentBedFluxFvPatchScalarField.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcVolumeIntegrate.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sedimentBedExchange, 0);
    addToRunTimeSelectionTable(functionObject, sedimentBedExchange, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::sedimentBedExchange::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Bed exchange of suspended sediment (volumes in m3)");
    writeCommented(os, "Time");
    writeTabbed(os, "intE");
    writeTabbed(os, "intD");
    writeTabbed(os, "intEminusD");
    writeTabbed(os, "intNetFaceFlux");
    writeTabbed(os, "otherOutflow");
    writeTabbed(os, "VC");
    writeTabbed(os, "deposit");
    writeTabbed(os, "budgetError");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sedimentBedExchange::sedimentBedExchange
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    CName_("C"),
    DCName_("DC"),
    phiName_("phi"),
    writeFields_(true),
    bedExchange_
    (
        IOobject
        (
            "bedExchange",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimVelocity, Zero)
    ),
    bedDeposit_
    (
        IOobject
        (
            "bedDeposit",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimLength, Zero)
    ),
    deposit_(getProperty<scalar>("deposit", 0)),
    otherOut_(getProperty<scalar>("otherOutflow", 0)),
    V0_(getProperty<scalar>("V0", 0)),
    initialised_(getProperty<bool>("initialised", false)),
    storage_(getProperty<scalar>("storage", 0)),
    startTimeIndex_(runTime.timeIndex()),
    ddtSchemeName_("Euler")
{
    read(dict);
    writeFileHeader(file());

    // Constructed at the first time-loop call, after the solver's fields:
    // the sediment volume at the start of the run for the budget error
    if (!initialised_ && foundObject<volScalarField>(CName_))
    {
        V0_ = fvc::domainIntegrate(lookupObject<volScalarField>(CName_)).value();
        initialised_ = true;
        setProperty("V0", V0_);
        setProperty("initialised", true);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sedimentBedExchange::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        CName_ = dict.getOrDefault<word>("field", "C");
        DCName_ = dict.getOrDefault<word>("DC", "DC");
        phiName_ = dict.getOrDefault<word>("phi", "phi");
        writeFields_ = dict.getOrDefault<bool>("writeFields", true);

        // Time scheme of the C equation, for the storage term of the budget
        const ITstream& is = mesh_.ddtScheme("ddt(" + CName_ + ')');
        ddtSchemeName_ =
        (
            is.size() && is[0].isWord() ? is[0].wordToken() : word("Euler")
        );
        if
        (
            ddtSchemeName_ != "Euler"
         && ddtSchemeName_ != "backward"
         && ddtSchemeName_ != "steadyState"
        )
        {
            WarningInFunction
                << "ddt scheme " << ddtSchemeName_ << " of " << CName_
                << ": the budget uses the Euler rule and closes only to the"
                << " time-discretisation error" << endl;
        }
        return true;
    }

    return false;
}


bool Foam::functionObjects::sedimentBedExchange::execute()
{
    const volScalarField& C = lookupObject<volScalarField>(CName_);
    const scalar VC = fvc::domainIntegrate(C).value();

    if (!initialised_)
    {
        FatalErrorInFunction
            << "Field " << CName_ << " not found at construction" << exit(FatalError);
    }

    const volScalarField& DC = lookupObject<volScalarField>(DCName_);
    const surfaceScalarField& phi = lookupObject<surfaceScalarField>(phiName_);
    const scalar dt = time().deltaTValue();

    scalar intE = 0, intD = 0, intNet = 0, other = 0;

    forAll(C.boundaryField(), patchi)
    {
        const fvPatchScalarField& Cp = C.boundaryField()[patchi];

        if (Cp.coupled())
        {
            continue;
        }

        const scalarField& magSf = Cp.patch().magSf();
        const scalarField diffIn(DC.boundaryField()[patchi]*Cp.snGrad());

        if (isA<sedimentBedFluxFvPatchScalarField>(Cp))
        {
            const auto& bed = refCast<const sedimentBedFluxFvPatchScalarField>(Cp);

            const scalarField E(bed.E());
            const scalarField D(bed.D());

            // Net upward flux from the face fluxes of the C equation:
            // diffusion in, settling (phi + ws.Sf) out
            const uniformDimensionedVectorField& ws =
                lookupObject<uniformDimensionedVectorField>("ws");
            const scalarField out
            (
                (phi.boundaryField()[patchi]/magSf + (ws.value() & Cp.patch().nf()))*Cp
            );

            intE += sum(E*magSf);
            intD += sum(D*magSf);

            // Mode nearestCell: the exchange is a source of the wall cell,
            // not a face flux
            if (bed.cellExchange())
            {
                intNet += sum((E - D)*magSf);
            }
            else
            {
                intNet += sum((diffIn - out)*magSf);
            }

            bedExchange_.boundaryFieldRef()[patchi] == (E - D);
            bedDeposit_.boundaryFieldRef()[patchi] ==
                bedDeposit_.boundaryField()[patchi] - (E - D)*dt;
        }
        else
        {
            other += sum(phi.boundaryField()[patchi]*Cp - diffIn*magSf);
        }
    }

    reduce(intE, sumOp<scalar>());
    reduce(intD, sumOp<scalar>());
    reduce(intNet, sumOp<scalar>());
    reduce(other, sumOp<scalar>());

    // Storage term with the time scheme of the C equation. Steps since the
    // start of this run: 0 = the call at the start (no step taken), 1 = the
    // first step (backward starts with Euler: no old-old field yet)
    const label step = time().timeIndex() - startTimeIndex_;

    if (step >= 1)
    {
        const scalar V1 = fvc::domainIntegrate(C.oldTime()).value();
        scalar dV = VC - V1;                                // Euler

        if (ddtSchemeName_ == "backward" && step >= 2 && C.nOldTimes() >= 2)
        {
            const scalar V2 =
                fvc::domainIntegrate(C.oldTime().oldTime()).value();
            const scalar dt0 = time().deltaT0Value();
            const scalar c = 1 + dt/(dt + dt0);
            const scalar c00 = dt*dt/(dt0*(dt + dt0));
            const scalar c0 = c + c00;
            dV = c*VC - c0*V1 + c00*V2;
        }

        storage_ += dV;
        deposit_ -= (intE - intD)*dt;
        otherOut_ += other*dt;
    }

    setProperty("storage", storage_);
    setProperty("deposit", deposit_);
    setProperty("otherOutflow", otherOut_);

    const scalar budgetError = storage_ + deposit_ + otherOut_;

    Log << type() << " " << name() << ": E = " << intE << ", D = " << intD
        << ", E - D = " << intE - intD << " (face fluxes " << intNet << ") m3/s"
        << ", V_C = " << VC << ", deposit = " << deposit_
        << ", budget error = " << budgetError << " m3" << endl;

    if (Pstream::master())
    {
        writeCurrentTime(file());
        file()
            << tab << intE << tab << intD << tab << intE - intD << tab << intNet
            << tab << other << tab << VC << tab << deposit_ << tab << budgetError
            << endl;
    }

    return true;
}


bool Foam::functionObjects::sedimentBedExchange::write()
{
    // Face fields at the solver's write times only (the time series is
    // written by execute at every step)
    if (writeFields_ && time().writeTime())
    {
        bedExchange_.instance() = time().timeName();
        bedDeposit_.instance() = time().timeName();
        bedExchange_.write();
        bedDeposit_.write();
    }

    return true;
}


// ************************************************************************* //
