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
    Foam::functionObjects::sedimentBedExchange

Group
    grpFieldFunctionObjects

Description
    Bed exchange of suspended sediment through the sedimentBedFlux patches:
    the net upward flux E - D per bed face (the source term of the Exner
    equation) and the sediment budget of the domain.

    Per bed face (fields written at write times, values on the bed patches,
    zero elsewhere):
    - bedExchange   E - D [m/s] of the last time step;
    - bedDeposit    cumulative deposit int (D - E) dt [m] (volume of grains
                    per bed area; divide by 1 - porosity for a bed level).

    Per time step (file postProcessing/<name>/<time>/bedExchange.dat):
    time, int E dA, int D dA, int (E - D) dA [m3/s], the same net flux from
    the face fluxes of the C equation (check), the outflow through all other
    non-coupled patches [m3/s], the sediment volume in the fluid
    V_C = int C dV [m3], the cumulative deposit int int (D - E) dA dt [m3]
    and the budget error

        sum_n (dV_C/dt)^n dt^n + deposit + cumulative other outflow,

    which is zero to solver tolerance for a conservative C equation. The
    storage term (dV_C/dt)^n is formed with the time scheme of the C
    equation (ddtSchemes entry for ddt(C)), so that the budget is the
    conservation statement of the discrete equations:
    - Euler (and steadyState): (V^n - V^n-1)/dt, i.e. the budget error is
      V_C - V_C(start) + deposit + other outflow;
    - backward: (c V^n - c0 V^n-1 + c00 V^n-2)/dt with the coefficients of
      OpenFOAM's backwardDdtScheme for variable dt (c = 1 + dt/(dt + dt0),
      c00 = dt^2/(dt0 (dt + dt0)), c0 = c + c00); the first step of a run
      is Euler, as in backwardDdtScheme (no old-old field yet);
    - other schemes (CrankNicolson, localEuler): Euler with a warning (the
      budget then closes only to the time-discretisation error).

    Execute after the C equation (the default execute control timeStep
    does this in suspendedSedimentFoam). The cumulative per-face deposit is
    read back from the time directory at a restart; the integrated totals
    restart from the stored function-object properties.

Usage
    \verbatim
    bedExchange
    {
        type        sedimentBedExchange;
        libs        (suspendedSediment);
        field       C;              // default C
        DC          DC;             // default DC
        phi         phi;            // default phi
        writeFields yes;            // default yes
    }
    \endverbatim

SourceFiles
    sedimentBedExchange.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_sedimentBedExchange_H
#define functionObjects_sedimentBedExchange_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class sedimentBedExchange Declaration
\*---------------------------------------------------------------------------*/

class sedimentBedExchange
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- Field names
        word CName_, DCName_, phiName_;

        //- Write the face fields
        bool writeFields_;

        //- E - D of the last step and cumulative deposit, per face
        volScalarField bedExchange_;
        volScalarField bedDeposit_;

        //- Integrated cumulative deposit and other outflow [m3]
        scalar deposit_;
        scalar otherOut_;

        //- Sediment volume at the start [m3] (kept for reference)
        scalar V0_;
        bool initialised_;

        //- Cumulative storage change sum (dV_C/dt) dt with the ddt scheme
        //  of the C equation [m3]
        scalar storage_;

        //- Time index at construction (first step of this run: Euler)
        label startTimeIndex_;

        //- Time scheme of the C equation (Euler, backward, ...)
        word ddtSchemeName_;


    // Private Member Functions

        //- Write the file header
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("sedimentBedExchange");


    // Constructors

        //- Construct from Time and dictionary
        sedimentBedExchange
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        sedimentBedExchange(const sedimentBedExchange&) = delete;

        //- No copy assignment
        void operator=(const sedimentBedExchange&) = delete;


    //- Destructor
    virtual ~sedimentBedExchange() = default;


    // Member Functions

        //- Read the settings
        virtual bool read(const dictionary&);

        //- Compute the exchange and the budget
        virtual bool execute();

        //- Write the face fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
