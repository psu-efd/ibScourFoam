# Free-surface tests (`ibScourInterFoam`)

Small 2-D cases, serial, seconds to minutes each. `fsChannelIB` needs the inlet and outlet conditions of
`libopenChannelBCs` from [openChannelVoF](https://github.com/psu-efd/openChannelVoF), built
separately; its `Allrun` reports SKIP without it. Every case has an `Allrun [MULES|isoAdvector]`
that meshes, initialises, runs and prints PASS/FAIL from `check.py`.

| Case | What it checks | Pass criterion |
|---|---|---|
| `fsDamBreakIB` | dam break over an immersed floor (plane at z = 0 on a cell face, two solid rows below) against the same case with the floor as a wall patch, run with stock `interFoam` (stage 1a step 2) | front position and column height within 2 % of the body-fitted run at 0.1--0.5 s, alpha L1 difference under 5 % of the water area (0.05--0.3 % before the splash at 0.3 s, 3--5 % after) |
| `fsChannelIB` | uniform open-channel flow over an immersed flat bed with the openChannelBCs discharge inlet and Manning normal-depth outlet (library built from the PI's singlePhaseVoF repository, loaded by `libs`), slope as tilted gravity, kOmegaSST with `density variable` (stage 1a step 3) | depth within 3 % of the Manning normal depth and flat to 2 %, bed shear of the IB row within 10 % of g h S, discharge conserved to 2 % (found: +0.1 %, 0.1 %, -7 %, 0.01 %) |
| `fsWaveIB` | Stokes II waves (H 0.02 m, T 1.2 s, depth 0.2 m) from the v2606 wave conditions over an immersed flat bed against the bed as a wall patch run by stock `interFoam`, surface elevation at four gauges (stage 1a step 4) | wave height at each gauge within 2 % of the body-fitted run and mean elevation difference under 5 % of the height over the last five periods (found: 0.00 % and 0.01 %) |
| `fsDamBreakMobile` | dam break (0.5 x 0.4 m column) over a dry mobile sand bed (height field at z = 0.05, d50 0.1 mm, bed load and suspended load, laminar) in a 2 x 0.9 m box (stage 1b): wet/dry gate on the bed shear, water transfer under a moving bed, sediment confined to the water | no sediment in air cells (alpha < 1e-3), C in cells with alpha < 0.5 bounded by the water value, suspended sediment equal to the integrated net entrainment within 10 %, water volume constant to 1e-3 up to 0.7 s (the splash leaves through the top later), no water lost under the bed, the dry bed inactive ahead of the bore (found: 0, ratio <= 1, 0.7 %, 2e-6, 0, 150 dry IB cells) |
| `fsLakeCut` | `fsLakeAtRest` with the sloping bed as a height-field bed with cut-face apertures (`cutFaces true`, z_b = 0.05 + 0.1 x as an `initialElevation` expression, 97 cut cells) and the immersed cylinder: the balance of the gravity and pressure fluxes on the open parts of the cut faces | same criteria as `fsLakeAtRest` (found: max 2.2e-8 m/s, last 1.2e-8, volume exact, hydrostatic to 1e-3) |
| `fsLakeAtRest` | water at rest over a sloping immersed bed with an immersed cylinder (stage 1a step 1): balanced gravity and pressure fluxes at the blocked faces, alpha of the solid cells, hydrostatic pressure | max spurious velocity in the water < 1e-5 m/s over 2 s and < 1e-6 m/s at the end, water volume constant to 1e-8, pressure hydrostatic to 2e-3 |

Lessons recorded by the cases:

- `maxDeltaT` must respect the gravity-wave Courant limit (sqrt(g h) dt/dx < 1): at 0.01 s on 1 cm cells
  the surface develops 3e-4 m/s oscillations, with stock interFoam as well; 0.002 s gives 1e-8 m/s.
- MULES turns the continuity error of a loosely converged pressure into volume-fraction dips next to the
  solids (alpha 0.99999 at tolerance 1e-8 / relTol 0.01), which drive a buoyancy current of 6e-5 m/s;
  the pressure is solved to 1e-10 (final 1e-12). isoAdvector is insensitive to this.
- Stage 1b: the bed shear of an IB cell counts only if its water fraction is at or above `wetThreshold` (0.5,
  bed dictionary); a column without a wet IB cell is inactive. The suspended load is solved for alpha C with
  the water flux of the alpha equation (MULES alphaPhi or the isoAdvector flux), the settling and diffusive
  fluxes weighted by the face water fraction, the bed exchange in wet IB cells only, and C set to zero where
  alpha < 1e-3; the cells buried by the bed hand their water to the cell above and fresh cells take the
  volume fraction and pressure of the cell above (logged every step).
- The fluid-subset pressure solve replaces derived pressure conditions (totalPressure,
  fixedFluxPressure) by mapped fixedValue and fixedGradient conditions on the subset copy.
