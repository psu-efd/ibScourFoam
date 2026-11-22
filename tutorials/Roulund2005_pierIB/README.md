# Roulund2005_pierIB: live-bed scour at a vertical pier (immersed bed)

Clear-water to live-bed scour around a vertical circular pier, after the
experiment and model of Roulund et al. (2005, J. Fluid Mech. 534, 351-401):
pier diameter D = 0.1 m, flow depth 0.2 m under a rigid lid, depth-mean
velocity 0.46 m/s, sand d50 = 0.26 mm.

| Item | Setting |
|---|---|
| Solver | `ibScourPimpleFoam` |
| Mesh | blockMesh + snappyHexMesh (pier as a wall patch), 8.0 M cells; the mesh floor lies 0.1 m below the bed so that every bed column has cells below the bed |
| Bed | immersed height-field bed on the floor patch `bottom`, cut-face apertures with cut-cell corrections, van Driest wall law (k_s = 0.7 mm) |
| Turbulence | k-omega SST |
| Sediment | Engelund-Fredsoe bed load with the slope-corrected critical Shields number (theta_c0 = 0.05), Koch-Flokstra gravitational transport, sand slide at 32 deg |
| Inflow | `constant/boundaryData/inlet`: developed profile from a precursor channel run |
| Time | 20 s of flow on the fixed bed, then morphological factor 2 to t = 170 s (300 s of morphological time) |

Run with `./Allrun` (`NP` ranks, default 128). Meshing takes about 15 min and
13 GB on one core; the flow run about 30 h on 128 ranks. The bed level is
written to `postProcessing/bed/<time>/bed.vtp` and the morphological time to
`postProcessing/bed/morphTime.dat`. `./Allclean` restores the initial state.

The same configuration with a body-fitted, deforming bed is
`../Roulund2005_pierBF`.
