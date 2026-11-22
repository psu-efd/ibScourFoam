# suspendedSedimentFoam

Solver for suspended sediment (one grain size, one-way coupled) over a
body-fitted bed on a static mesh, OpenFOAM v2606. Derived from `pimpleFoam`
(OpenFOAM Foundation, OpenCFD Ltd; notices kept in the derived files),
GPL-3.0-or-later. The bed exchange is the `sedimentBedFlux` condition of
`libsuspendedSediment`; the same formulation serves the immersed bed of
`ibScourPimpleFoam` and `ibScourInterFoam` (`ibSuspendedLoad`).

## Equations

Flow: `pimpleFoam` (incompressible, RANS/LES/laminar, fvOptions), static
mesh. Sediment (`CEqn.H`), after the PIMPLE loop of each time step:

    ddt(C) + div(phi + ws·Sf, C) − laplacian(DC, C) = fvOptions(C),
    DC = nu/ScC + nut/sigmaC

- `ws = |ws| g/|g|` from a `settlingVelocityModel`, registered as the
  uniform field `ws`; `DC` is a registered field (written).
- The settling part of the drift flux is kept only on patches where C has
  the condition `sedimentBedFlux` (bed: deposition and entrainment). On
  every other non-coupled patch (lid, walls, inlet, outlet) `phiC = phi`.
- For `sedimentBedFlux` with `diffusivity firstCell|logLaw` the bed-face
  value of DC is replaced by the condition's `faceDiffusivity()`.
- No density coupling (one-way).

## Input: constant/sedimentProperties

    solveFlow       yes;        // no: frozen U, phi, nut (C only)
    solveSediment   yes;        // no: flow only (precursor)
    diffusivity     turbulence; // or field: DC read from 0/DC, kept fixed
    ScC             1;          // nu/D_m
    sigmaC          1;          // turbulent Schmidt number
    d50             2e-4;       // m
    s               2.65;       // rho_s/rho
    settlingVelocityModel  vanRijn1984;   // constant (ws), stokes, soulsby
    referenceConcentrationModel  zysermanFredsoe;   // vanRijn1984 (delta), constant (Cb, delta)

plus `constant/g`, `transportProperties`, `turbulenceProperties`, fields
`U`, `p`, `C` (and the turbulence fields). The function object
`sedimentBedExchange` (libsuspendedSediment) writes the bed exchange and the
sediment budget.

## Build and tests

    ./Allwmake       # at the repository root (library and solver)

Tests: `tests/suspendedLoad/cases/settlingColumn` (deposition, budget to
7e-15), `settlingColumnBed` and `rouseColumn` (see `tests/README.md`).
