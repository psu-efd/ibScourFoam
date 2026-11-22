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


Description
    Construct the immersed boundary, register its turbulence constraint and
    include the immersed solids in the wall distance.

\*---------------------------------------------------------------------------*/

immersedBoundary& ib =
    const_cast<immersedBoundary&>(immersedBoundary::New(mesh));

// Add the IB turbulence constraint unless the case already defines one or
// the flow is laminar (no k, omega or epsilon fields)
{
    PtrList<fv::option>& optionList = fvOptions;

    bool found = false;

    forAll(optionList, i)
    {
        if (isA<fv::ibTurbulenceConstraint>(optionList[i]))
        {
            found = true;
        }
    }

    const bool turbulent =
    (
        mesh.foundObject<volScalarField>("k")
     || mesh.foundObject<volScalarField>("omega")
     || mesh.foundObject<volScalarField>("epsilon")
    );

    if (!found && turbulent)
    {
        dictionary constraintDict;
        constraintDict.add("type", fv::ibTurbulenceConstraint::typeName);

        optionList.push_back
        (
            fv::option::New("ibTurbulence", constraintDict, mesh)
        );
    }
}

// Turbulence models use the wall distance; include the immersed solids
ib.correctWallDistance();

// Registered turbulent viscosity field (absent for laminar flow)
volScalarField* nutPtr = mesh.getObjectPtr<volScalarField>("nut");

// Pressure gradient of the momentum equation and the velocity correction
// (blocked faces at the cell value). Written for post-processing, where
// fvc::grad(p) would use the decoupled solid-cell pressure in the IB cells.
volVectorField ibGradP
(
    IOobject
    (
        "ibGradP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    ib.grad(p)
);

// ************************************************************************* //
