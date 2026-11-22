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

Application
    kBudget

Description
    Terms of the k-omega SST turbulent kinetic energy equation per unit
    volume, from the written fields U, k, omega, nut and phi of an
    incompressible run (pimpleFoam or ibScourPimpleFoam), with the
    discretisation of kOmegaSSTBase::correct():

        kConv  = div(phi, k)
        kProd  = min(G, c1 betaStar k omega),  G = nut (grad U && dev2 ...)
        kDiss  = betaStar omega k
        kDiff  = laplacian(DkEff, k),  DkEff = alphaK(F1) nut + nu
        kDivU  = (2/3) (div U) k
        kGbyNu = G/nut (strain-rate part of the production)
        kResidual = kConv - kDiff - kProd + kDiss + kDivU

    The residual is minus the time derivative in a free cell and the source
    imposed by a constraint in a pinned cell (wall-law or immersed-boundary
    cells). The wall distance is wallDist of the mesh, limited by the
    immersed-boundary distance 'ibDistance' when that field is present, as
    in immersedBoundary::correctWallDistance(). The SST coefficients are the
    OpenFOAM defaults.

Usage
    kBudget -time 15

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

// k-omega SST coefficients (OpenFOAM defaults)
const Foam::scalar betaStar = 0.09;
const Foam::scalar c1 = 10.0;
const Foam::scalar alphaK1 = 0.85;
const Foam::scalar alphaK2 = 1.0;
const Foam::scalar alphaOmega2 = 0.856;

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Terms of the k-omega SST k equation from the written fields"
    );

    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createMesh.H"

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dimensionedScalar nu("nu", dimViscosity, transportProperties);

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        auto read = [&](const word& name)
        {
            return IOobject
            (
                name,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            );
        };

        volVectorField U(read("U"), mesh);
        volScalarField k(read("k"), mesh);
        volScalarField omega(read("omega"), mesh);
        volScalarField nut(read("nut"), mesh);
        surfaceScalarField phi(read("phi"), mesh);

        // Wall distance: mesh walls, limited by the immersed boundary
        volScalarField y
        (
            IOobject("kBudget:y", runTime.timeName(), mesh),
            wallDist::New(mesh).y()
        );

        Info<< "    mesh wall distance (wall patches only): min "
            << gMin(y.primitiveField()) << " max "
            << gMax(y.primitiveField()) << endl;

        {
            IOobject ibHeader
            (
                "ibDistance",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            );

            if (ibHeader.typeHeaderOk<volScalarField>(true))
            {
                Info<< "    wall distance limited by ibDistance" << endl;

                const volScalarField ibDistance(ibHeader, mesh);
                scalarField& yi = y.primitiveFieldRef();
                const scalarField& d = ibDistance.primitiveField();

                forAll(yi, celli)
                {
                    yi[celli] = max(min(yi[celli], d[celli]), SMALL);
                }

                y.correctBoundaryConditions();
            }
        }

        y.max(dimensionedScalar(dimLength, SMALL));

        Info<< "    wall distance used for F1: min "
            << gMin(y.primitiveField()) << " max "
            << gMax(y.primitiveField()) << endl;

        // Production (kOmegaSSTBase::GbyNu0, Pk)
        tmp<volTensorField> tgradU = fvc::grad(U);
        volScalarField kGbyNu
        (
            IOobject("kGbyNu", runTime.timeName(), mesh),
            tgradU() && devTwoSymm(tgradU())
        );
        tgradU.clear();

        volScalarField kProd
        (
            IOobject("kProd", runTime.timeName(), mesh),
            min(nut*kGbyNu, (c1*betaStar)*k*omega)
        );

        volScalarField kDiss
        (
            IOobject("kDiss", runTime.timeName(), mesh),
            betaStar*omega*k
        );

        // Diffusion with the SST blending (kOmegaSSTBase::F1, DkEff)
        const volScalarField CDkOmega
        (
            (2*alphaOmega2)*(fvc::grad(k) & fvc::grad(omega))/omega
        );

        const volScalarField CDkOmegaPlus
        (
            max(CDkOmega, dimensionedScalar(dimless/sqr(dimTime), 1.0e-10))
        );

        const volScalarField arg1
        (
            min
            (
                min
                (
                    max
                    (
                        (scalar(1)/betaStar)*sqrt(k)/(omega*y),
                        scalar(500)*nu/(sqr(y)*omega)
                    ),
                    (4*alphaOmega2)*k/(CDkOmegaPlus*sqr(y))
                ),
                scalar(10)
            )
        );

        const volScalarField F1(tanh(pow4(arg1)));

        const volScalarField DkEff
        (
            "DkEff",
            (F1*(alphaK1 - alphaK2) + alphaK2)*nut + nu
        );

        volScalarField kDiff
        (
            IOobject("kDiff", runTime.timeName(), mesh),
            fvc::laplacian(DkEff, k)
        );

        volScalarField kConv
        (
            IOobject("kConv", runTime.timeName(), mesh),
            fvc::div(phi, k)
        );

        volScalarField kDivU
        (
            IOobject("kDivU", runTime.timeName(), mesh),
            (2.0/3.0)*fvc::div(phi)*k
        );

        volScalarField kResidual
        (
            IOobject("kResidual", runTime.timeName(), mesh),
            kConv - kDiff - kProd + kDiss + kDivU
        );

        volScalarField kF1(IOobject("kF1", runTime.timeName(), mesh), F1);

        kGbyNu.write();
        kProd.write();
        kDiss.write();
        kDiff.write();
        kConv.write();
        kDivU.write();
        kResidual.write();
        kF1.write();

        Info<< "    wrote kProd kDiss kDiff kConv kDivU kGbyNu kResidual kF1"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
