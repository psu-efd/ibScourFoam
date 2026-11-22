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
    Foam::nutWallLawFvPatchScalarField

Group
    grpWallFunctions

Description
    Wall-function condition for the turbulent viscosity of a body-fitted
    wall that uses the wall law of the immersed boundary (ibWallLaw), so
    that a body-fitted wall and an immersed wall differ only in their
    geometric representation.

    u_tau follows from the tangential velocity of the first cell relative
    to the wall, |U_t|, at the normal distance y of its centre (see
    wallLawFaces); the wall value

        nu_t,w = max(u_tau^2 y/|U_t| - nu, 0)

    makes the face transmit the wall-law stress u_tau^2. Use it with the
    wallLawTurbulence constraint (omega and k of the wall cells) and
    kqRWallFunction for k.

Usage
    \verbatim
    bed
    {
        type            nutWallLaw;
        wallLaw
        {
            type        vanDriest;  // spalding, roughLogLaw, vanDriest, ...
            ks          0.00065;    // equivalent sand roughness [m]
        }
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::ibWallLaw
    Foam::fv::wallLawTurbulence

SourceFiles
    nutWallLawFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_nutWallLawFvPatchScalarField_H
#define Foam_nutWallLawFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"
#include "ibWallLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class nutWallLawFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutWallLawFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
    typedef nutWallLawFvPatchScalarField this_bctype;
    typedef nutWallFunctionFvPatchScalarField parent_bctype;

    // Private Data

        //- Wall-law settings
        dictionary lawDict_;

        //- The wall law
        autoPtr<ibWallLaw> law_;


protected:

    // Protected Member Functions

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> calcNut() const;


public:

    //- Runtime type information
    TypeName("nutWallLaw");


    // Constructors

        //- Construct from patch and internal field
        nutWallLawFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutWallLawFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        nutWallLawFvPatchScalarField
        (
            const this_bctype&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        nutWallLawFvPatchScalarField
        (
            const this_bctype&,
            const DimensionedField<scalar, volMesh>&
        );

        //- No copy without an internal field
        nutWallLawFvPatchScalarField(const this_bctype&) = delete;

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member Functions

        //- The wall law
        const ibWallLaw& law() const
        {
            return *law_;
        }

        //- Calculate and return the yPlus at the boundary
        virtual tmp<scalarField> yPlus() const;

        //- Write
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Inherit assignment
        using parent_bctype::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
