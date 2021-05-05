/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "kOmegaSSTSASIB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSSTSASIB<BasicTurbulenceModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    volScalarField::Internal L
    (
        sqrt(this->k_())/(pow025(this->betaStar_)*this->omega_())
    );

    volScalarField::Internal Lvk
    (
        max
        (
            kappa_*sqrt(S2)
           /(
                mag(fvc::laplacian(this->U_))()()
              + dimensionedScalar
                (
                    "ROOTVSMALL",
                    dimensionSet(0, -1, -1, 0, 0),
                    ROOTVSMALL
                )
            ),
            Cs_*sqrt(kappa_*zeta2_/(beta/this->betaStar_ - gamma))*delta()()
        )
    );

    return fvm::Su
    (
        this->alpha_()*this->rho_()
       *min
        (
            max
            (
                zeta2_*kappa_*S2*sqr(L/Lvk)
              - (2*C_/sigmaPhi_)*this->k_()
               *max
                (
                    magSqr(fvc::grad(this->omega_)()())/sqr(this->omega_()),
                    magSqr(fvc::grad(this->k_)()())/sqr(this->k_())
                ),
                dimensionedScalar("0", dimensionSet(0, 0, -2, 0, 0), 0)
            ),
            // Limit SAS production of omega for numerical stability,
            // particularly during start-up
            this->omega_()/(0.1*this->omega_.time().deltaT())
        ),
        this->omega_
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTSASIB<BasicTurbulenceModel>::kOmegaSSTSASIB
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSST<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.11
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    zeta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "zeta2",
            this->coeffDict_,
            3.51
        )
    ),
    sigmaPhi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaPhi",
            this->coeffDict_,
            2.0/3.0
        )
    ),
    C_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C",
            this->coeffDict_,
            2
        )
    ),

    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", alphaRhoPhi.group()),
            *this,
            this->coeffDict_
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaSSTSASIB<BasicTurbulenceModel>::read()
{
    if (kOmegaSST<BasicTurbulenceModel>::read())
    {
        Cs_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
        sigmaPhi_.readIfPresent(this->coeffDict());
        zeta2_.readIfPresent(this->coeffDict());
        C_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kOmegaSSTSASIB<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    surfaceScalarField alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    BasicTurbulenceModel::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu(dev(twoSymm(tgradU()())) && tgradU()());
    volScalarField::Internal G(this->GName(), nut()*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    this->omega_.boundaryFieldRef().updateCoeffs();

    // create ibMesh
    const immersedBoundaryFvMesh& ibMesh = const_cast<immersedBoundaryFvMesh&>
            (this->db().parent().objectRegistry::lookupObject<immersedBoundaryFvMesh>(this->db().name()));
            
    

    // main correction of IB method
    ibMesh.kOmegaCorrection(*this); // move it before momentum predictor
    
     
    //ibMesh.ibCorrectPhi(alphaRhoPhi);

    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, this->omega_)
            + fvm::div(alphaRhoPhi, this->omega_)
            - fvm::laplacian(alpha*rho*this->DomegaEff(F1), this->omega_)
            ==
            alpha()*rho()*gamma
            *min
            (
                GbyNu,
                (this->c1_/this->a1_)*this->betaStar_*this->omega_()
                *max(this->a1_*this->omega_(), this->b1_*F23()*sqrt(S2()))
            )
            - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, this->omega_)
            - fvm::Sp(alpha()*rho()*beta*this->omega_(), this->omega_)
            - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/this->omega_(),
                this->omega_
            )
            + Qsas(S2(), gamma, beta)
            + this->omegaSource()
            + fvOptions(alpha, rho, this->omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
    //    correctDiag IB method
   //    ibMesh.correctDiag(omegaEqn.ref());
        
        omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
        ibMesh.manipulateMatrix(omegaEqn.ref());
        solve(omegaEqn);
        fvOptions.correct(this->omega_);
        bound(this->omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, this->k_)
        + fvm::div(alphaRhoPhi, this->k_)
        - fvm::laplacian(alpha*rho*this->DkEff(F1), this->k_)
        ==
        alpha()*rho()*this->Pk(G)
        - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, this->k_)
        - fvm::Sp(alpha()*rho()*this->epsilonByk(F1, F23),this->k_)
        + this->kSource()
        + fvOptions(alpha, rho, this->k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    ibMesh.manipulateMatrix(kEqn.ref());
    
    // correctOffDiag IB method
  //  ibMesh.evaluatek();
  //  ibMesh.correctDiag(kEqn.ref());
  //  ibMesh.correctOffDiag(kEqn.ref());
   
   
    solve(kEqn);
  
  
   // ibMesh.evaluatek();
    fvOptions.correct(this->k_);
    bound(this->k_, this->kMin_);
   
    this->correctNut(S2, F23);
    
    
    ibMesh.nutCorrection();
   
    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
