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

Namespace
    Foam::suspendedSediment

Description
    Grain properties, settling velocity and reference concentration of
    libsuspendedSediment.

Class
    Foam::suspendedSediment::settlingVelocityModel

Description
    Run-time selectable settling velocity of a single grain size, read
    from constant/sedimentProperties:
    \verbatim
    d50                     2.0e-4;     // grain diameter [m]
    s                       2.65;       // rho_s/rho [-]
    settlingVelocityModel   vanRijn1984;
    ws                      0.022;      // model "constant" only [m/s]
    \endverbatim

    Models (same laws as python/settling_velocity.py of the sediment cases):
    \table
        Model        | Law                                     | Source
        constant     | ws given                                | -
        stokes       | (s-1) g d^2/(18 nu)                     | van Rijn (1984b) Eq. 11
        vanRijn1984  | Stokes / Zanke / 1.1 sqrt((s-1) g d),   | van Rijn (1984b) Eqs. 11-13
                     | limits 100 and 1000 micron              |
        soulsby      | nu/d [(10.36^2 + 1.049 D*^3)^0.5 - 10.36] | Soulsby (1997)
    \endtable

SourceFiles
    settlingVelocityModel.C

\*---------------------------------------------------------------------------*/

#ifndef suspendedSediment_settlingVelocityModel_H
#define suspendedSediment_settlingVelocityModel_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace suspendedSediment
{

/*---------------------------------------------------------------------------*\
                    Class settlingVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class settlingVelocityModel
{
protected:

    // Protected Data

        //- Grain diameter [m]
        const scalar d50_;

        //- Relative density rho_s/rho [-]
        const scalar s_;


public:

    //- Runtime type information
    TypeName("settlingVelocityModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        settlingVelocityModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );


    // Constructors

        //- Construct from sedimentProperties
        explicit settlingVelocityModel(const dictionary& dict);

    //- Select from the entry "settlingVelocityModel" of sedimentProperties
    static autoPtr<settlingVelocityModel> New(const dictionary& dict);

    //- Destructor
    virtual ~settlingVelocityModel() = default;


    // Member Functions

        //- Settling velocity [m/s] for the kinematic viscosity nu [m2/s]
        //  and the acceleration of gravity g [m/s2]
        virtual scalar ws(const scalar nu, const scalar g) const = 0;

        //- Dimensionless grain size D* = d [(s-1) g/nu^2]^(1/3)
        scalar Dstar(const scalar nu, const scalar g) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace suspendedSediment
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
