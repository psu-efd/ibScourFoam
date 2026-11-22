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
    Foam::seepageModel

Description
    Effect of the seepage through the bed on sediment mobility, selected at
    run time from the 'type' entry of the 'seepage' sub-dictionary of the
    bed dictionary (default: none). Given the upward hydraulic gradient i
    at a bed face (positive for upward seepage), the model returns the
    factors applied to the flat-bed critical Shields number theta_c0 and
    to the friction coefficient tan(phi_r) of the slope correction and of
    the sand slide. The critical gradient i_c = gamma'/gamma_w =
    (s - 1)(1 - n) is that of the bed (overridden by 'criticalGradient').

Usage
    \verbatim
    seepage
    {
        type            chengChiew;     // none | chengChiew
        reposeEffect    true;           // also reduce tan(phi_r)
        // criticalGradient 1.0;        // i_c (default (s - 1)(1 - n))
    }
    \endverbatim

SourceFiles
    seepageModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_seepageModel_H
#define Foam_seepageModel_H

#include "dictionary.H"
#include "scalarField.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class seepageModel Declaration
\*---------------------------------------------------------------------------*/

class seepageModel
{
protected:

    // Protected Data

        //- Settings
        dictionary dict_;

        //- Critical hydraulic gradient i_c
        const scalar iCrit_;


public:

    //- Runtime type information
    TypeName("seepageModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            seepageModel,
            dictionary,
            (const dictionary& dict, const scalar iCrit),
            (dict, iCrit)
        );


    // Constructors

        //- Construct from dictionary and the critical gradient of the bed
        seepageModel(const dictionary& dict, const scalar iCrit);


    // Selectors

        //- Return the selected model; 'none' if the dictionary is absent
        static autoPtr<seepageModel> New
        (
            const dictionary& bedDict,
            const scalar iCrit
        );


    //- Destructor
    virtual ~seepageModel() = default;


    // Member Functions

        //- Critical hydraulic gradient
        scalar iCrit() const noexcept { return iCrit_; }

        //- True if the model changes the mobility (false for 'none')
        virtual bool active() const { return true; }

        //- Factor on the flat-bed critical Shields number for the upward
        //- gradient i
        virtual scalar thetaCFactor(const scalar i) const = 0;

        //- Factor on the friction coefficient tan(phi_r)
        virtual scalar reposeFactor(const scalar i) const = 0;

        //- Factor on the Shields number of the transport closure (the
        //- generalised Shields number of Francalanci et al. 2008: the
        //- bed shear is scaled by the seepage-reduced effective weight);
        //- default 1
        virtual scalar thetaFactor(const scalar i) const
        {
            return 1;
        }

        //- True if the horizontal seepage force is applied as an apparent
        //- bed slope (pressureGradient model)
        virtual bool horizontalForce() const
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
