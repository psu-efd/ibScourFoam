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
    Foam::fv::ibTurbulenceConstraint

Group
    grpFvOptionsConstraints

Description
    Fixes the turbulence quantities of the immersed boundary: in IB cells
    k, omega and epsilon take the wall-model values (ibWallModel); in solid
    cells they are held at small positive values.

    With 'kTreatment transported' the turbulent kinetic energy of the IB
    cells is not fixed but solved, as in a body-fitted wall-function cell:
    the coupling of the k equation across faces between IB and solid cells
    is removed (zero gradient towards the wall), while omega (epsilon)
    keeps the wall-model value.

    The solver ibScourPimpleFoam adds this option automatically when it is
    not present in constant/fvOptions or system/fvOptions; its settings are
    then read from the 'turbulence' sub-dictionary of
    constant/immersedBoundaryDict.

Usage
    \verbatim
    ibTurbulence
    {
        type            ibTurbulenceConstraint;
        fields          (k omega);  // optional; default: those present
        kTreatment      equilibrium;// equilibrium (default) | transported
    }
    \endverbatim

SourceFiles
    ibTurbulenceConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fv_ibTurbulenceConstraint_H
#define Foam_fv_ibTurbulenceConstraint_H

#include "fvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class ibTurbulenceConstraint Declaration
\*---------------------------------------------------------------------------*/

class ibTurbulenceConstraint
:
    public fv::option
{
    // Private Data

        //- Solve k in IB cells (true) or fix it to the wall-model value
        bool transportedK_;



public:

    //- Runtime type information
    TypeName("ibTurbulenceConstraint");


    // Constructors

        //- Construct from components
        ibTurbulenceConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        ibTurbulenceConstraint(const ibTurbulenceConstraint&) = delete;

        //- No copy assignment
        void operator=(const ibTurbulenceConstraint&) = delete;


    //- Destructor
    virtual ~ibTurbulenceConstraint() = default;


    // Member Functions

        //- Fix the IB and solid cell values of a scalar equation
        virtual void constrain(fvMatrix<scalar>& eqn, const label fieldi);

        //- Read the settings
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
