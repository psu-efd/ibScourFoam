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

Class
    Foam::suspendedSediment::referenceConcentrationModel

Description
    Run-time selectable equilibrium (reference) volume concentration
    C_b* at the reference level delta above the bed, as a function of the
    bed shear velocity, read from constant/sedimentProperties:
    \verbatim
    d50                         7.9e-5;     // [m]
    s                           2.65;       // rho_s/rho [-]
    referenceConcentrationModel zysermanFredsoe;
    delta                       1.58e-4;    // reference level [m]
                                            // (default 2 d50 for ZF)
    shieldsNumber               skinFriction;   // or total (default
                                            // skinFriction for ZF and vR)
    kN                          5e-4;       // skin roughness [m] (default
                                            // 2.5 d50 for ZF, 3 d90 for vR)
    d90                         3.2e-4;     // vanRijn1984 skin friction
    Cb                          0.032;      // model "constant" only
    \endverbatim

    Models (same laws as python/reference_concentration.py):
    \table
        Model            | C_b*                                        | delta
        zysermanFredsoe  | 0.331 x/(1 + 0.72 x), x = (theta' - 0.045)^1.75 | 2 d50 (default)
        vanRijn1984      | 0.015 (d50/a) T^1.5/D*^0.3 (van Rijn 1984b Eq. 38) | a (required)
        constant         | Cb given                                    | required
    \endtable
    with theta' = u*'^2/((s-1) g d50) and T = (u*'^2 - u*cr^2)/u*cr^2,
    u*cr from the Shields curve of van Rijn (1984a).

    Shear velocity in the Shields number (entry shieldsNumber):
    - skinFriction (default for zysermanFredsoe and vanRijn1984): the
      grain (skin-friction) shear velocity u*' of the original papers,
      from the depth-mean velocity V and the depth D of the flow above the
      bed face (supplied by the caller, sedimentBedFlux) and the total u*:
      - zysermanFredsoe: Zyserman & Fredsoe (1994) Eqs. 9-10,
        V/u*' = 6 + 2.5 ln(D'/k_N) and u*' = (g D' I)^(1/2), k_N = 2.5 d50
        (p. 1023 and procedure step 1, p. 1038); with g D I = u*^2 (uniform
        flow, wide channel) D' = D (u*'/u*)^2, and u*' solves
        V = u*' [6 + 2.5 ln(D u*'^2/(u*^2 k_N))];
      - vanRijn1984: u*' = g^(1/2) V/C', C' = 18 log10(12 D/k_N),
        k_N = 3 d90 (van Rijn 1984a,b; R_b = D for a wide channel);
      in both u*' <= u* (skin friction is part of the total; a bed
      smoother than k_N gives u*' = u*);
    - total: the total shear velocity of the flow solution (the wall shear
      of the bed face), as used by Liu (2014, 2018).

SourceFiles
SourceFiles
    referenceConcentrationModel.C

\*---------------------------------------------------------------------------*/

#ifndef suspendedSediment_referenceConcentrationModel_H
#define suspendedSediment_referenceConcentrationModel_H

#include "scalarField.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace suspendedSediment
{

/*---------------------------------------------------------------------------*\
                  Class referenceConcentrationModel Declaration
\*---------------------------------------------------------------------------*/

class referenceConcentrationModel
{
protected:

    // Protected Data

        //- Grain diameter [m]
        const scalar d50_;

        //- Relative density rho_s/rho [-]
        const scalar s_;

        //- Reference level above the bed [m]
        scalar delta_;

        //- Shields number from the skin-friction shear velocity
        //  (shieldsNumber skinFriction) instead of the total one
        bool skinFriction_;

        //- Skin (grain) roughness k_N [m] (set by the models)
        scalar kN_;


public:

    //- Runtime type information
    TypeName("referenceConcentrationModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        referenceConcentrationModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );


    // Constructors

        //- Construct from sedimentProperties
        explicit referenceConcentrationModel(const dictionary& dict);

    //- Select from the entry "referenceConcentrationModel"
    static autoPtr<referenceConcentrationModel> New(const dictionary& dict);

    //- Destructor
    virtual ~referenceConcentrationModel() = default;


    // Member Functions

        //- Reference level [m]
        scalar delta() const
        {
            return delta_;
        }

        //- True if the model needs the bed shear velocity
        virtual bool needsShear() const
        {
            return true;
        }

        //- True if C_b* uses the skin-friction shear velocity, i.e. the
        //  caller must supply the depth-mean velocity and the depth
        bool skinFriction() const
        {
            return skinFriction_;
        }

        //- Skin roughness k_N [m]
        scalar kN() const
        {
            return kN_;
        }

        //- Skin-friction shear velocity u*' [m/s] from the total shear
        //  velocity uStar, the depth-mean velocity V [m/s] and the depth
        //  D [m] of the flow above each bed face (base class: uStar)
        virtual tmp<scalarField> uStarSkin
        (
            const scalarField& uStar,
            const scalarField& V,
            const scalarField& D,
            const scalar g
        ) const;

        //- Shields parameter u*^2/((s-1) g d50)
        tmp<scalarField> theta(const scalarField& uStar, const scalar g) const;

        //- Reference volume concentration for the shear velocity uStar [m/s]
        //  (u*' with shieldsNumber skinFriction, else the total u*)
        virtual tmp<scalarField> Cb
        (
            const scalarField& uStar,
            const scalar nu,
            const scalar g
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace suspendedSediment
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
