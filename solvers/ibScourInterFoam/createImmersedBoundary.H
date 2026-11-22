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
    include the immersed solids in the wall distance (as ibScourPimpleFoam).

\*---------------------------------------------------------------------------*/

immersedBoundary& ib =
    const_cast<immersedBoundary&>(immersedBoundary::New(mesh));

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

ib.correctWallDistance();

volScalarField* nutPtr = mesh.getObjectPtr<volScalarField>("nut");

// The volume fraction of the solid cells follows the fluid next to them
ib.correctSolidAlpha(alpha1);
alpha2 = 1.0 - alpha1;
rho == alpha1*rho1 + alpha2*rho2;
rho.oldTime() = rho;
mixture.correct();

// Optional ghost-free volume fraction (solid cells zero) for function objects
// that integrate alpha, e.g. interfaceHeight gauges: the solid cells carry the
// ghost extension of correctSolidAlpha, which is not water
autoPtr<volScalarField> alphaFluidPtr;
if (ib.dict().getOrDefault<bool>("fluidAlphaField", false))
{
    alphaFluidPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                alpha1.name() + "Fluid",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            alpha1
        )
    );
    ib.fluidOnly(alpha1, *alphaFluidPtr);
    Info<< "Writing the ghost-free volume fraction " << alphaFluidPtr->name()
        << " for function objects" << endl;
}
