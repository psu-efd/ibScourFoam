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
    Foam::darcyStorage

Description
    Level 1 soil model: pore-pressure diffusion in a rigid skeleton with a
    compressible pore fluid (storage equation),

        (n/K') dp/dt = div((k/gamma_w) grad p),

    with porosity n, pore-fluid bulk modulus K' (set by the degree of
    saturation S_r: 1/K' = 1/K_w + (1 - S_r)/p_abs), hydraulic conductivity
    k and unit weight of water gamma_w. The hydraulic diffusivity is
    c_v = k K'/(n gamma_w). The pressure is fixed on the top (bed) patch
    from the flow, zero-flux on the sides and, unless 'bottomPressure' is
    given, on the bottom; a fixed bottom pressure drives a steady seepage
    (upward for a positive value) through the column.

Usage
    \verbatim
    soil
    {
        type            darcyStorage;
        porosity        0.4;
        conductivity    1e-4;       // k [m/s]
        saturation      0.99;       // S_r; or bulkModulus K' [Pa] directly
        absolutePressure 1.2e5;     // p_abs [Pa] in 1/K' (default)
        waterBulkModulus 2.2e9;     // K_w [Pa] (default)
        waterDensity    1000;
        gravity         9.81;
        p0              0;          // initial pore pressure [Pa]
        // bottomPressure constant 0; // Function1 of time [Pa], optional
    }
    \endverbatim

SourceFiles
    darcyStorage.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_darcyStorage_H
#define Foam_darcyStorage_H

#include "soilModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class darcyStorage Declaration
\*---------------------------------------------------------------------------*/

class darcyStorage
:
    public soilModel
{
    // Private Data

        //- Porosity
        const scalar porosity_;

        //- Hydraulic conductivity k [m/s]
        const scalar conductivity_;

        //- Pore-fluid bulk modulus K' [Pa]
        const scalar bulkModulus_;

        //- Storage coefficient n/K' [1/Pa]
        const dimensionedScalar storage_;

        //- Mobility k/gamma_w [m3 s/kg]
        const dimensionedScalar mobility_;

        //- Optional fixed pressure on the bottom patch [Pa]
        autoPtr<Function1<scalar>> bottomPressure_;

        //- Pore pressure [Pa]
        volScalarField p_;

        //- Restart: the pore pressure was read from the time directory
        bool restored_;


    // Private Member Functions

        //- Pore-fluid bulk modulus from the dictionary
        static scalar readBulkModulus(const dictionary& dict);

        //- Patch field types of p
        wordList patchTypes() const;


public:

    //- Runtime type information
    TypeName("darcyStorage");


    // Constructors

        //- Construct from the soil mesh and dictionary
        darcyStorage(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~darcyStorage() = default;


    // Member Functions

        //- Pore pressure [Pa]
        virtual const volScalarField& p() const
        {
            return p_;
        }

        //- Set the pore pressure of every cell
        virtual void initialise(const scalarField& pCells);

        //- Restart: the state field was read from the time directory
        virtual bool restored() const { return restored_; }

        //- Remap the pore pressure after the columns moved
        virtual void remap
        (
            const labelList& c0,
            const labelList& c1,
            const scalarField& w
        );

        //- Advance one flow time step with the pore pressure pTop [Pa]
        //- on the faces of the top patch
        virtual void correct(const scalarField& pTop, const label nSteps);

        //- Hydraulic diffusivity c_v = k K'/(n gamma_w) [m2/s]
        virtual scalar diffusivity() const
        {
            return conductivity_*bulkModulus_/(porosity_*gammaW());
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
