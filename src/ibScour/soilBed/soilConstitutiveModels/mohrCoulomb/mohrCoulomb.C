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

#include "mohrCoulomb.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "soilModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace soilConstitutiveModels
{
    defineTypeNameAndDebug(mohrCoulomb, 0);
    addToRunTimeSelectionTable(soilConstitutiveModel, mohrCoulomb, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::symmTensor Foam::soilConstitutiveModels::mohrCoulomb::elasticStress
(
    const symmTensor& eps
) const
{
    return 2*mu0_*eps + lambda0_*tr(eps)*symmTensor::I;
}


Foam::symmTensor Foam::soilConstitutiveModels::mohrCoulomb::elasticStrain
(
    const symmTensor& sigma
) const
{
    return
        sigma/(2*mu0_)
      - lambda0_*tr(sigma)/(2*mu0_*(3*lambda0_ + 2*mu0_))*symmTensor::I;
}


Foam::label Foam::soilConstitutiveModels::mohrCoulomb::returnMapping
(
    vector& s
) const
{
    // s = (s1, s2, s3) ordered s1 >= s2 >= s3, tension positive
    const scalar k = k_;
    const scalar m = m_;
    const scalar sqrtk = Foam::sqrt(k);
    const scalar f = k*s.x() - s.z() - 2*c_*sqrtk;

    if (f <= 0)
    {
        return 0;
    }

    // Elastic stiffness in principal space
    const scalar a = lambda0_ + 2*mu0_;
    const scalar b = lambda0_;

    auto Dmul = [&](const vector& v) -> vector
    {
        return vector
        (
            a*v.x() + b*v.y() + b*v.z(),
            b*v.x() + a*v.y() + b*v.z(),
            b*v.x() + b*v.y() + a*v.z()
        );
    };

    // Apex (all planes meet): for phi = 0 (k = 1) there is no apex
    const bool hasApex = (k > 1 + SMALL);
    const scalar sApex = (hasApex ? 2*c_*sqrtk/(k - 1) : GREAT);

    // 1. Return to the plane f1 = k s1 - s3 - 2 c sqrt(k): sigma_p =
    //    sigma - lambda D r_g with r_f = (k, 0, -1), r_g = (m, 0, -1)
    const vector rf(k, 0, -1);
    const vector rg(m, 0, -1);
    const vector Drg(Dmul(rg));
    {
        const scalar lam = f/(rf & Drg);
        const vector sp(s - lam*Drg);

        if (sp.x() >= sp.y() - SMALL && sp.y() >= sp.z() - SMALL)
        {
            s = sp;
            return 1;
        }
    }

    // 2. Return to an edge: two planes active, sigma_p = sigma - D(l1 rg1
    //    + l2 rg2), with f1(sigma_p) = f2(sigma_p) = 0 (linear in l1, l2)
    auto edgeReturn =
        [&](const vector& rf2, const vector& rg2, vector& sp) -> bool
    {
        const vector Drg2(Dmul(rg2));
        const scalar f2 = (rf2 & s) - 2*c_*sqrtk;
        // [rf.Drg  rf.Drg2 ] [l1]   [f ]
        // [rf2.Drg rf2.Drg2] [l2] = [f2]
        const scalar a11 = rf & Drg, a12 = rf & Drg2;
        const scalar a21 = rf2 & Drg, a22 = rf2 & Drg2;
        const scalar det = a11*a22 - a12*a21;
        if (mag(det) < SMALL)
        {
            return false;
        }
        const scalar l1 = (f*a22 - f2*a12)/det;
        const scalar l2 = (a11*f2 - a21*f)/det;
        if (l1 < 0 || l2 < 0)
        {
            return false;
        }
        sp = s - l1*Drg - l2*Drg2;
        return true;
    };

    {
        vector sp;
        // Edge s1 = s2 >= s3: planes f1 (k s1 - s3) and f(k s2 - s3)
        if (edgeReturn(vector(0, k, -1), vector(0, m, -1), sp))
        {
            const scalar s12 = 0.5*(sp.x() + sp.y());
            if (s12 >= sp.z() - SMALL && (!hasApex || sp.z() <= sApex + SMALL))
            {
                s = vector(s12, s12, sp.z());
                return 2;
            }
        }
        // Edge s1 >= s2 = s3: planes f1 (k s1 - s3) and (k s1 - s2)
        if (edgeReturn(vector(k, -1, 0), vector(m, -1, 0), sp))
        {
            const scalar s23 = 0.5*(sp.y() + sp.z());
            if (sp.x() >= s23 - SMALL && (!hasApex || sp.x() <= sApex + SMALL))
            {
                s = vector(sp.x(), s23, s23);
                return 2;
            }
        }
    }

    // 3. Apex
    if (hasApex)
    {
        s = vector(sApex, sApex, sApex);
        return 3;
    }

    // phi = 0 and no admissible edge: project on the plane (Tresca)
    const scalar lam = f/(rf & Drg);
    s -= lam*Drg;
    return 1;
}


Foam::label Foam::soilConstitutiveModels::mohrCoulomb::returnPoint
(
    const symmTensor& eps,
    const symmTensor& s0,
    const symmTensor& epsPOld,
    symmTensor& sig,
    symmTensor& epsP
) const
{
    // Trial total effective stress
    const symmTensor sigTrial(s0 + elasticStress(eps - epsPOld));

    // Principal stresses, ordered s1 >= s2 >= s3 (eigenValues returns
    // them in ascending order, eigenVectors as rows)
    const vector ev(eigenValues(sigTrial));
    const tensor evec(eigenVectors(sigTrial, ev));
    vector s(ev.z(), ev.y(), ev.x());

    const label region = returnMapping(s);

    if (region == 0)
    {
        sig = sigTrial - s0;
        epsP = epsPOld;
        return 0;
    }

    const vector n1(evec.z()), n2(evec.y()), n3(evec.x());
    const symmTensor sigRet
    (
        s.x()*symm(n1*n1) + s.y()*symm(n2*n2) + s.z()*symm(n3*n3)
    );

    sig = sigRet - s0;
    epsP = eps - elasticStrain(sig);

    return region;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::soilConstitutiveModels::mohrCoulomb::mohrCoulomb
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    linearElastic(mesh, dict),
    c_(dict.get<scalar>("cohesion")),
    phi_(degToRad(dict.get<scalar>("frictionAngle"))),
    psi_(degToRad(dict.getOrDefault<scalar>("dilationAngle", 0))),
    k_((1 + Foam::sin(phi_))/(1 - Foam::sin(phi_))),
    m_((1 + Foam::sin(psi_))/(1 - Foam::sin(psi_))),
    mu0_(E_/(2*(1 + nu_))),
    lambda0_(E_*nu_/((1 + nu_)*(1 - 2*nu_))),
    epsPOld_
    (
        IOobject
        (
            "epsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor(dimless, Zero)
    ),
    epsP_
    (
        IOobject("epsilonPWork", mesh.time().timeName(), mesh),
        mesh,
        dimensionedSymmTensor(dimless, Zero)
    ),
    nPlane_(0),
    nEdge_(0),
    nApex_(0)
{
    Info<< "    Mohr-Coulomb soil: c " << c_ << " Pa, phi "
        << radToDeg(phi_) << " deg, psi " << radToDeg(psi_) << " deg"
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::soilConstitutiveModels::mohrCoulomb::correct
(
    const volTensorField& gradD,
    const volSymmTensorField& sigma0,
    volSymmTensorField& sigma
)
{
    const symmTensorField eps(symm(gradD.primitiveField()));
    const symmTensorField& s0 = sigma0.primitiveField();
    const symmTensorField& epsPOld = epsPOld_.primitiveField();
    symmTensorField& epsP = epsP_.primitiveFieldRef();
    symmTensorField& sig = sigma.primitiveFieldRef();

    nPlane_ = nEdge_ = nApex_ = 0;

    forAll(sig, celli)
    {
        const label region =
            returnPoint
            (
                eps[celli],
                s0[celli],
                epsPOld[celli],
                sig[celli],
                epsP[celli]
            );

        if (region == 1) ++nPlane_;
        else if (region == 2) ++nEdge_;
        else if (region == 3) ++nApex_;
    }

    // Boundary values from the patch displacement gradient and the
    // plastic strain of the adjacent cells
    auto& sigBf = sigma.boundaryFieldRef();

    forAll(sigBf, patchi)
    {
        const symmTensorField epsB(symm(gradD.boundaryField()[patchi]));
        const symmTensorField& s0B = sigma0.boundaryField()[patchi];
        const symmTensorField epsPB
        (
            epsPOld_.boundaryField()[patchi].patchInternalField()
        );
        symmTensorField& sigB = sigBf[patchi];
        symmTensor epsPDummy;

        forAll(sigB, facei)
        {
            returnPoint
            (
                epsB[facei],
                s0B[facei],
                epsPB[facei],
                sigB[facei],
                epsPDummy
            );
        }
    }

    if (debug && (nPlane_ + nEdge_ + nApex_) > 0)
    {
        Info<< "    Mohr-Coulomb: " << nPlane_ << " plane, " << nEdge_
            << " edge, " << nApex_ << " apex returns" << endl;
    }
}


void Foam::soilConstitutiveModels::mohrCoulomb::updateState()
{
    epsPOld_ = epsP_;

    const label n =
        returnReduce(nPlane_ + nEdge_ + nApex_, sumOp<label>());

    if (n > 0)
    {
        Info<< "    Mohr-Coulomb: " << n << " plastic cells ("
            << returnReduce(nPlane_, sumOp<label>()) << " plane, "
            << returnReduce(nEdge_, sumOp<label>()) << " edge, "
            << returnReduce(nApex_, sumOp<label>()) << " apex)" << endl;
    }
}


void Foam::soilConstitutiveModels::mohrCoulomb::remap
(
    const labelList& c0,
    const labelList& c1,
    const scalarField& w
)
{
    epsPOld_.primitiveFieldRef() =
        soilModel::remapField(epsPOld_.primitiveField(), c0, c1, w);
    epsPOld_.correctBoundaryConditions();
    epsP_ = epsPOld_;
}


// ************************************************************************* //
