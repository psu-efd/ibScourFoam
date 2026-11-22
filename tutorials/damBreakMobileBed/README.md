# damBreakMobileBed: dam break over a mobile bed

The IAHR benchmark of Soares-Frazao et al. (2012, J. Hydraul. Res. 50(4),
364-375): a 1 m wide breach between two impervious blocks releases a
reservoir into a 3.6 m wide flume over an 85 mm layer of uniform sand
(d50 = 1.61 mm) on a rigid floor.

* **Case 1:** reservoir at 0.47 m, dry sand downstream.
* **Case 2:** reservoir at 0.51 m, sand under 0.15 m of water downstream.

| Item | Setting |
|---|---|
| Solver | `ibScourInterFoam` (VOF free surface, MULES) |
| Domain | half flume (symmetry plane through the breach), reservoir and upstream tank included; the dam blocks are body-fitted walls |
| Bed | immersed height-field bed over the whole floor with a non-erodible floor (`zRigid`: the rigid flume floor and the downstream sill) |
| Turbulence | k-omega SST, rough wall law k_s = 2.5 d50 |
| Sediment | Engelund-Fredsoe bed load, critical Shields number 0.047 with slope correction, Koch-Flokstra transport, sand slide at 32 deg |
| Time | 20 s of flow, gate removed instantaneously |

`make_case.py` writes the case; `./Allrun` runs the whole sequence (set-up,
blockMesh, setFields, bed fields, decomposition, run). Options: `CASE=1|2`,
`MESH=coarse|fine` (5 or 2.5 cm cells over the sand; 0.49 M or 1.97 M cells),
`NP` (default 32). Run times on Intel Ice Lake: about 15 min (coarse, 32 ranks)
and 50 min (fine, 64 ranks). Water levels at the gauges of the paper are
written to `postProcessing/gauges`, the bed to `postProcessing/bed`.

The measured data are not part of this repository; they are in the paper.
