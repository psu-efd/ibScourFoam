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

Class
    Foam::soilConstitutiveModels::mohrCoulomb

Description
    Elastic-perfectly plastic Mohr-Coulomb law with a non-associated flow
    rule, returned in principal stress space (Clausen, Damkilde &
    Andersen 2006, 2007). With the principal effective stresses ordered
    sigma_1 >= sigma_2 >= sigma_3 (tension positive) the yield function
    and the plastic potential are

        f = k sigma_1 - sigma_3 - 2 c sqrt(k),  k = (1 + sin phi)/(1 - sin phi)
        g = m sigma_1 - sigma_3,                m = (1 + sin psi)/(1 - sin psi)

    with cohesion c, friction angle phi and dilation angle psi. The trial
    stress sigma0 + C:(eps - eps_p) is returned to the yield plane, to one
    of the two edges (sigma_1 = sigma_2 or sigma_2 = sigma_3) or to the
    apex, whichever keeps the principal stresses ordered, and the plastic
    strain is updated from the returned stress. The working plastic
    strain is committed by updateState() at the end of a converged time
    step.

    References:
        Clausen, J., Damkilde, L., & Andersen, L. (2006). Efficient return
        algorithms for associated plasticity with multiple yield planes.
        Int. J. Numer. Methods Eng., 66(6), 1036-1059.
        Clausen, J., Damkilde, L., & Andersen, L. (2007). An efficient
        return algorithm for non-associated plasticity with linear yield
        criteria in principal stress space. Comput. Struct., 85(23-24),
        1795-1807.

Usage
    \verbatim
    constitutive
    {
        type            mohrCoulomb;
        E               5e6;            // [Pa]
        nu              0.3;
        cohesion        5e3;            // c [Pa]
        frictionAngle   0;              // phi [deg]
        dilationAngle   0;              // psi [deg] (default 0)
    }
    \endverbatim

SourceFiles
    mohrCoulomb.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_soilMohrCoulomb_H
#define Foam_soilMohrCoulomb_H

#include "linearElastic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace soilConstitutiveModels
{

/*---------------------------------------------------------------------------*\
                         Class mohrCoulomb Declaration
\*---------------------------------------------------------------------------*/

class mohrCoulomb
:
    public linearElastic
{
    // Private Data

        //- Cohesion [Pa]
        const scalar c_;

        //- Friction and dilation angles [rad]
        const scalar phi_;
        const scalar psi_;

        //- k = (1 + sin phi)/(1 - sin phi), m likewise with psi
        const scalar k_;
        const scalar m_;

        //- Lame constants [Pa] (uniform)
        const scalar mu0_;
        const scalar lambda0_;

        //- Committed and working plastic strain
        volSymmTensorField epsPOld_;
        volSymmTensorField epsP_;

        //- Number of cells returned (plane, edge, apex) at the last call
        label nPlane_;
        label nEdge_;
        label nApex_;


    // Private Member Functions

        //- Elastic stress of a strain (uniform moduli)
        symmTensor elasticStress(const symmTensor& eps) const;

        //- Elastic strain of a stress
        symmTensor elasticStrain(const symmTensor& sigma) const;

        //- Return the trial principal stresses (ordered s1 >= s2 >= s3)
        //- to the yield surface; returns the region (0 elastic, 1 plane,
        //- 2 edge, 3 apex)
        label returnMapping(vector& s) const;

        //- Stress increment sig and plastic strain epsP of a point with
        //- the strain eps, geostatic stress s0 and committed plastic
        //- strain epsPOld; returns the region
        label returnPoint
        (
            const symmTensor& eps,
            const symmTensor& s0,
            const symmTensor& epsPOld,
            symmTensor& sig,
            symmTensor& epsP
        ) const;


public:

    //- Runtime type information
    TypeName("mohrCoulomb");


    // Constructors

        //- Construct from the soil mesh and dictionary
        mohrCoulomb(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~mohrCoulomb() = default;


    // Member Functions

        virtual void correct
        (
            const volTensorField& gradD,
            const volSymmTensorField& sigma0,
            volSymmTensorField& sigma
        );

        virtual void updateState();

        virtual void remap
        (
            const labelList& c0,
            const labelList& c1,
            const scalarField& w
        );

        //- Committed plastic strain
        const volSymmTensorField& plasticStrain() const { return epsPOld_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace soilConstitutiveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
