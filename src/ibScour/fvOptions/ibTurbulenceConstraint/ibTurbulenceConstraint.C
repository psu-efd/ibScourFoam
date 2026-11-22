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

#include "ibTurbulenceConstraint.H"
#include "immersedBoundary.H"
#include "ibWallModel.H"
#include "fvMatrices.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(ibTurbulenceConstraint, 0);
    addToRunTimeSelectionTable(option, ibTurbulenceConstraint, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::ibTurbulenceConstraint::ibTurbulenceConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(name, modelType, dict, mesh),
    transportedK_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::ibTurbulenceConstraint::constrain
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const immersedBoundary& ib = immersedBoundary::New(mesh_);
    const ibWallModel& wm = ib.wallModel();

    const word& fieldName = eqn.psi().name();

    // Wall-model values in the IB cells
    scalarField ibValues;

    if (fieldName == "k")
    {
        ibValues = wm.kIB();
    }
    else if (fieldName == "omega")
    {
        ibValues = wm.omegaIB();
    }
    else if (fieldName == "epsilon")
    {
        ibValues = wm.epsilonIB();
    }
    else
    {
        FatalErrorInFunction
            << "Field " << fieldName << " is not supported; use k, omega"
            << " or epsilon" << exit(FatalError);
    }

    const labelList& ibCells = ib.ibCells();
    const labelList& solidCells = ib.solidCells();

    // Transported k: the IB cells solve their own k with zero gradient
    // towards the solid; only the solid cells are fixed
    const bool fixIB = !(transportedK_ && fieldName == "k");

    if (!fixIB)
    {
        ib.decoupleSolid(eqn);
    }

    labelList cells((fixIB ? ibCells.size() : 0) + solidCells.size());
    scalarField values(cells.size());

    label n = 0;

    if (fixIB)
    {
        forAll(ibCells, i)
        {
            cells[n] = ibCells[i];
            values[n++] = ibValues[i];
        }
    }

    // Small positive values in the solid cells
    const scalar solidValue = (fieldName == "omega" ? 1 : SMALL);

    forAll(solidCells, i)
    {
        cells[n] = solidCells[i];
        values[n++] = solidValue;
    }

    eqn.setValues(cells, values);
}


bool Foam::fv::ibTurbulenceConstraint::read(const dictionary& dict)
{
    if (!fv::option::read(dict))
    {
        return false;
    }

    if (!coeffs_.readIfPresent("fields", fieldNames_))
    {
        fieldNames_.clear();

        for (const word f : {"k", "omega", "epsilon"})
        {
            if (mesh_.foundObject<volScalarField>(f))
            {
                fieldNames_.push_back(f);
            }
        }
    }

    // Settings from the option itself or, when added by the solver, from
    // constant/immersedBoundaryDict (turbulence sub-dictionary)
    const immersedBoundary& ib = immersedBoundary::New(mesh_);
    const dictionary& ibTurbDict =
        ib.dict().subOrEmptyDict("turbulence");

    const word kTreatment
    (
        coeffs_.getOrDefault<word>
        (
            "kTreatment",
            ibTurbDict.getOrDefault<word>("kTreatment", "equilibrium")
        )
    );

    if (kTreatment != "equilibrium" && kTreatment != "transported")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown kTreatment " << kTreatment
            << "; use equilibrium or transported" << exit(FatalIOError);
    }

    transportedK_ = (kTreatment == "transported");

    Info<< "    ibTurbulenceConstraint applied to " << fieldNames_
        << ", k treatment " << kTreatment << endl;

    fv::option::resetApplied();

    return true;
}

// ************************************************************************* //
