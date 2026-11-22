/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Xiaofeng Liu, Penn State University
-------------------------------------------------------------------------------
    Part of suspendedSedimentFoam (no code from OpenFOAM in this file).
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

Description
    Sediment fields and settings of suspendedSedimentFoam, read from
    constant/sedimentProperties:

    \verbatim
    solveFlow       yes;        // no: frozen U, phi, nut (C only)
    solveSediment   yes;        // no: flow only (e.g. a precursor run)
    diffusivity     turbulence; // DC = nu/ScC + nut/sigmaC; or "field":
                                // DC read from 0/DC and kept fixed
    ScC             1;          // molecular Schmidt number nu/D_m [-]
    sigmaC          1;          // turbulent Schmidt number nut/eps_t [-]
    d50             2e-4;       // grain diameter [m]
    s               2.65;       // relative density rho_s/rho [-]
    settlingVelocityModel  vanRijn1984;    // or constant (ws), stokes, soulsby
    referenceConcentrationModel  zysermanFredsoe; // used by sedimentBedFlux
    \endverbatim

    The settling velocity is the vector ws = |ws| g/|g|, registered as the
    uniform field "ws" for the bed condition sedimentBedFlux; the
    diffusivity is the registered field "DC" (the patch condition and the
    bedExchange function object use the same field as the C equation).

\*---------------------------------------------------------------------------*/

#include "readGravitationalAcceleration.H"

Info<< "Reading sedimentProperties\n" << endl;
IOdictionary sedimentProperties
(
    IOobject
    (
        "sedimentProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const bool solveFlow(sedimentProperties.getOrDefault<bool>("solveFlow", true));
const bool solveSediment
(
    sedimentProperties.getOrDefault<bool>("solveSediment", true)
);
const bool fixedDiffusivity
(
    sedimentProperties.getOrDefault<word>("diffusivity", "turbulence") == "field"
);

const dimensionedScalar ScC("ScC", dimless, sedimentProperties);
const dimensionedScalar sigmaC("sigmaC", dimless, sedimentProperties);

// Settling velocity, along g
const scalar nu0 = gAverage(laminarTransport.nu()().primitiveField());
const scalar wsMag =
    suspendedSediment::settlingVelocityModel::New(sedimentProperties)
   ->ws(nu0, mag(g).value());

uniformDimensionedVectorField ws
(
    IOobject
    (
        "ws",
        runTime.constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    dimensionedVector(dimVelocity, wsMag*g.value()/mag(g.value()))
);

Info<< "Settling velocity ws = " << ws.value() << " m/s (|ws| = " << wsMag
    << ", nu = " << nu0 << " m2/s)" << nl
    << "solveFlow " << solveFlow << ", solveSediment " << solveSediment
    << ", diffusivity " << (fixedDiffusivity ? "field (0/DC)" : "nu/ScC + nut/sigmaC")
    << nl << endl;

Info<< "Reading field C (sediment volume concentration)\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Diffusivity of the C equation (registered: used by sedimentBedFlux).
// diffusivity field: read from 0/DC (READ_IF_PRESENT reads the file)
volScalarField DC
(
    IOobject
    (
        "DC",
        runTime.timeName(),
        mesh,
        fixedDiffusivity ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimViscosity, Zero),
    fvPatchFieldBase::calculatedType()
);

if (fixedDiffusivity && !DC.headerOk())
{
    FatalErrorInFunction
        << "diffusivity field: 0/DC not found" << exit(FatalError);
}

if (!fixedDiffusivity)
{
    DC = turbulence->nu()/ScC + turbulence->nut()/sigmaC;
}

// ************************************************************************* //
