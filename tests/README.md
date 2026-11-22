# ibScourFoam test suite

Small cases that run in minutes (serial or up to 8 ranks). Each one
isolates one mechanism of the immersed-boundary (IB) scour solvers
`ibScourPimpleFoam` and `ibScourInterFoam` (library `libibScour`) or of the
body-fitted bed (`deformingBed`). They serve as regression and verification
tests.

## Layout

```
tests/
├── run_tests.sh            run selected cases (locally or as a batch job), collect metrics
├── common/                 shared Python modules and shell functions
│   ├── runFunctions.sh     environment, ibs_solve (serial/mpirun + reconstruct), cleaning
│   ├── foam_file.py        ASCII OpenFOAM I/O (headers, polyMesh, fields, dictionaries)
│   ├── bed_field.py        height-field bed level (<time>/finite-area/zb) and floor geometry
│   ├── decompose_zb.py     decomposition of the finite-area bed fields
│   ├── reconstruct_zb.py   reconstruction of the bed level of a parallel run
│   ├── stl_tools.py        height-field/cylinder STL writers, STL reader, slopes, volumes
│   ├── channel_profiles.py laminar / turbulent open-channel inflow profiles
│   ├── case_setup.py       template copy, boundaryData, initial fields
│   ├── log_tools.py        solver-log parser (Courant, continuity, timings, crashes)
│   ├── anomalies.py        first-anomaly times (Courant, continuity, p iterations, dz, tau, NaN)
│   ├── metrics.py          metrics.json writer (PASS / FAIL / REPORT)
│   ├── make_view.py        ParaView view/ of a run (time dirs with fields + view.foam)
│   └── collect.py          summary tables over all metrics.json
└── <case>/                 Allrun, Allclean, check.py and the case folders:
    └── <variant>/          a complete OpenFOAM case (0.orig, constant, system,
                            Allrun with blockMesh, setExprFields, makeFaMesh and
                            the solver; case_setup.json holds the parameters
                            that check.py uses)
```

The suite cases (`sandSlideRelax`, `ibChannelFlat`, `humpMigration`,
`pipelineScour`) consist of ready-made case folders, one per variant. Their
`Allrun` (or `run_tests.sh`) copies each folder to
`$IBS_RUN_ROOT/<case>/<variant>` (default `<repo>/run/tests`, ignored by git),
runs the folder's own `Allrun` there, prepares a ParaView view and writes
`metrics.json` with `check.py`. A variant can also be run by hand: copy its
folder anywhere and call its `Allrun`. Initial flow fields with the inflow
profile are set by `setExprFields` (`system/setExprFieldsDict`). The
self-contained cases (`buriedCylinder2D`, `rigidFloor`, `suspendedLoad/cases/*`,
`freeSurface/*`) run in place with their own `Allrun`.

## Bed level and area mesh

The height-field bed lives on a finite-area mesh of the floor patch:

* **Area mesh.** With `system/finite-area/faMeshDefinition` (as in the
  OpenFOAM finite-area tutorials) `makeFaMesh` writes it to
  `constant/finite-area`; in parallel, `decomposePar -no-finite-area` and
  `mpirun makeFaMesh -parallel` decompose it with the bed fields. Without the
  definition the solver builds it at run time (serial cases need nothing
  else). `reconstructPar` of OpenFOAM v2606 corrupts the heap on
  finite-area fields (also on the stock `liquidFilmFoam` tutorial), so the
  parallel cases reconstruct with `redistributePar -reconstruct -parallel`.
* **Initial bed** (`bed` dictionary): `initialElevation` is a constant or any
  `PatchFunction1` on the floor patch: an `expression` of the face centre,
  or `mappedFile` for scattered survey points
  (`constant/boundaryData/<floorPatch>/points` and `.../0/<fieldTable>`,
  2-D Delaunay interpolation with `mapMethod planar`); `initialSurface` takes
  an STL/OBJ bed surface from `constant/triSurface` instead. A field file
  `0/finite-area/zb` takes precedence. The non-erodible floor takes the same
  forms (`rigidElevation`, `rigidSurface`, `0/finite-area/zRigid`).

`bedInput` tests each option, serial and on 4 ranks.

## Requirements

* OpenFOAM v2606 with ibScourFoam built (`./Allwmake` at the repository root).
* Python >= 3.8 with numpy (`$IBS_PYTHON`, default `python3`). The comparison and plotting
  scripts of `inclinedBed`, `deformingBedHump` and `common/plot_bed_shear.py`, and
  `bodyFittedPier --bed-vtp`, also need matplotlib, vtk and scipy.

## Running

```bash
# Load OpenFOAM first (or set IBS_ENV to a script that does)
tests/run_tests.sh ibChannelFlat sandSlideRelax humpMigration
tests/run_tests.sh pipelineScour:2d_hydro,2d_scour

# Selected variants
tests/run_tests.sh ibChannelFlat:lam_f025,ras_tilt2

# Print the tables of existing results
tests/run_tests.sh --collect all

# A single case directly
tests/pipelineScour/Allrun 2d_scour
SPINUP=0.2 MORPH=0.3 tests/pipelineScour/Allrun 2d_scour   # shorter
```

Environment variables: `IBS_RUN_ROOT`, `IBS_ENV`, `IBS_PYTHON`,
`IBS_MAX_JOBS` (concurrent serial runs), `IBS_SOLVER` (default
`ibScourPimpleFoam`), `IBS_TIMEOUT` (solver wall-time limit), `IBS_MPI_BIND`
(default `--bind-to core`).

**Do not edit a case's `Allrun` while a job is executing it.** Bash reads
scripts incrementally, so editing a running script corrupts the job.

## Status semantics

`check.py` writes `metrics.json` with one entry per metric:

* `PASS` / `FAIL` – the metric has a tolerance (listed below).
* `REPORT` – informative only.
* overall `ERROR` – the run crashed or produced no output.

## Cases

| Case | Mechanism tested | Size | Runtime (one core) |
|---|---|---|---|
| **`pipelineScour` (primary)** | half-buried horizontal cylinder: IB bed meeting a body-fitted structure, hydro (fixed bed) and scour (mobile bed); regression and robustness case | 2D slab 24 k cells (serial); true 2D 12 k (hydro); 3D 0.19 M (8 ranks) | 2D: ≈ 15 min per variant |
| T1 `ibChannelFlat` | immersed plane bed: no-slip/wall law, mass leakage, τ_b accuracy and noise vs bed position in the cell | 22–36 k cells, serial, 10 IB + 10 body-fitted runs | 1–2 min per run |
| T2 `sandSlideRelax` | sand slide alone (bed-only run): relaxation to the angle of repose, volume conservation | 35–45 k cells, serial, 20 steps | seconds |
| T3 `humpMigration` | Exner + bed load + IB update on a migrating hump (quasi-2D) | 25 k cells, serial, 1800 steps, hump run + flat control | ≈ 5 min per run |
| `bedInput` | bed-level input: expression, scattered survey points (mappedFile), STL surface, non-erodible floor; area mesh from makeFaMesh; 4-rank decompose, run and reconstruct with the stock tools | 32 k cells, bed only, 2 steps per variant | ≈ 1 min |
| `buriedCylinder2D` | solid overlap of the height field: a half-buried immersed cylinder in a 2D bed-only channel (prescribed shear, no flow); the sand under the cylinder is released as bed change, not sediment; columns under a solid are inactive until exposed; deposition stops at the underside | 2.9 k cells, serial and 4 ranks, 400 steps | ≈ 5 s |

### `pipelineScour` (primary debugging case)

A half-buried pipeline: a horizontal cylinder of D = 0.05 m, axis across
the flow at the bed level, so half of it is below the immersed sand bed
(z = 0). Depth h = 4D = 0.2 m, U_b = 0.55 m/s (flat-bed θ ≈ 0.11 ≈ 2θ_c0,
d50 = 0.26 mm, live bed). The domain spans 10 D upstream to 20 D downstream
and the mesh reaches 1 D below the bed. The cylinder is body-fitted (wall
patch `cylinder`, wall functions); the bed is a height-field IB bed on the
floor patch `bottom` (Spalding wall law, Engelund–Fredsøe bed load,
Koch–Flokstra transport, sand slide, equilibrium-feed strip at the inlet).

* **2D mesh** (`2d_*`): blockMesh O-grid (9 blocks,
  128 cells around the cylinder, 3.1 mm cells across the bed band),
  in `system/blockMeshDict`, shifted up by Δz/4 so that no mesh face lies in
  the bed plane. `--span-cells N` gives a slab of N cells with `symmetry`
  sides (default 2, 24 k cells); `2d1_*` is the true 2D mesh (one cell,
  `empty` sides).
* **3D mesh** (`3d_*`): a finite cylinder (length 3 D = 0.15 m,
  capped) in a 0.3 m wide channel with `symmetry` sides. snappyHexMesh is
  used (base 10 mm, level 1 near the bed, level 2 at the cylinder; no
  layers): 0.19 M cells, run on 8 ranks. This variant adds 3-D effects
  (end vortices, lateral flow around the pipe ends).
* **Procedure** (`Allrun`): one run of `SPINUP` + `MORPH` seconds (1 + 2 s,
  Δt = 1 ms); the bed is fixed throughout (`*_hydro`) or mobile after the
  spin-up (`*_scour`, `morphStartTime`, morphological factor `MF` = 20).
* **Variants**: `2d_hydro`, `2d_scour`, `2d1_hydro`, `3d_hydro`, `3d_scour`.

| Metric | Tolerance |
|---|---|
| `q_imbalance` | ≤ 1e-3 |
| `leak_gross_rel` (Σ\|φ\| fluid–dead / Q_in, last write) | ≤ 1e-3 |
| `tau_cv_upstream` (IB τ_b on the undisturbed bed) | ≤ 0.05 |
| `log_volume_error_max` (Exner volume error of the bed report) | ≤ 1e-12 m³ |
| `spike_rel` (second-difference residual / D) | ≤ 0.02 |
| `span_spread_rel` (across-span bed range / D; REPORT in 3D) | ≤ 0.02 (2D) |
| `tau_inf`, `tau_amp_max`, `tau_amp_max_x_over_D`, `tau_lee_max` (τ_b / τ∞ near the pipe) | REPORT |
| `scour_depth_final`, `scour_over_D_final` (pipe flanks D/2 < \|x\| < 2D, span-mean), `volume_change_rel`, `slope_max_deg` | REPORT |
| `anomaly_order`, `s_per_step` | REPORT |

In 2D the pipe without a gap blocks the near-bed flow: τ_b falls to almost
zero at both toes and recovers to only 0.5 τ∞ in the lee by x ≈ 1.2 D, so
the bed next to the pipe stays near the threshold and no scour hole forms
in 40 s of morphological time (the far field is live). The 2D variants are
therefore a robustness and conservation test (both pass: discharge to
3e-9, no flow through the bed, volume error 5e-23 m³, spike 0.017 D). In
the 3D variant τ_b reaches 2.6 τ∞ around the pipe ends and the bed erodes
there (6 mm after 40 s of morphological time); `3d_scour` passes on 8 ranks
(volume error 4e-22 m³, spike 0.013 D).

Parallel runs: the cases have a `faMeshDefinition`, so `ibs_solve` uses
`decomposePar -no-finite-area`, `makeFaMesh -parallel` and
`redistributePar -reconstruct` (see "Bed level and area mesh"). For cases
without one it falls back to `common/decompose_zb.py` and
`common/reconstruct_zb.py`.

### T1 `ibChannelFlat`

A straight channel with a rigid slip lid, slip side walls and a flat bed.
**IB variant:** a box mesh from 4 cells below the bed to the lid; the bed is
an immersed plane solid (`constant/immersedBoundaryDict`) at `z_b = f·Δz`
(+ an optional cross-stream tilt). **Body-fitted (BF) reference:** the same
channel with the bed as a wall patch, solved with stock `pimpleFoam`. The
inflow profile (`constant/boundaryData/inlet`, zero below the bed) is also
used as the initial field.

* `lam_*`: laminar, h = 0.05 m, U_b = 0.01 m/s (Re_h = 500), Δz = h/20,
  analytical half-Poiseuille profile, τ_b/ρ = 3νU_b/h; laminar wall law.
* `ras_*`: kOmegaSST in both runs (Spalding wall law in the IB run, wall
  functions in the BF run), h = 0.1 m, U_b = 0.3 m/s, Δz = h/20 (first-cell
  y⁺ ≈ 35).
* Bed position: `f000`, `f025`, `f050`, `f075` (fraction of a cell
  height; at `f050` the cell centres lie exactly on the bed) and `tilt2`
  (2° cross-stream tilt around f = 0.5, so that the bed cuts the cells at
  every fraction across the span; tilting across the flow keeps the depth
  uniform along the channel).

| Metric | Tolerance | Meaning |
|---|---|---|
| `q_imbalance` | ≤ 1e-3 | \|Q_out/Q_in − 1\| (outlet flux through fluid + dead cells) |
| `Q_out_dead` | REPORT | outflow through below-bed (dead) outlet cells |
| `leak_gross_rel` | ≤ 1e-3 | Σ\|φ\| on fluid–dead faces / Q_in (flow through the bed) |
| `ib_ib_gross_rel` | REPORT | Σ\|φ\| on IB–IB faces / Q_in |
| `u_rms_err_analytic`, `u_rms_err_bf` | ≤ 0.05 | RMS (u − u_ref)/U_b over live cells, 0.5–0.9 L (lam) / 0.6–0.9 L (ras) |
| `u_ib_cells_rel_err_mean` | REPORT | mean (u_IB − u_BF)/U_b in IB cells |
| `tau_cv` | ≤ 0.05 | std/mean of τ_b over bed columns after removing a plane fit |
| `tau_ratio_err_analytic` / `_bf` | ≤ 0.10 (lam), 0.15 (ras) | \|τ̄_IB/τ̄_ref − 1\| |

τ_b (kinematic, m²/s²) for RAS is \|`ibWallShearStress`\| in IB cells. For
laminar flow it is ν ∂u/∂z at the bed from a quadratic through u(z_b) = 0
and the first two fluid cells with centres more than 0.1 Δz above the bed.

All ten variants pass: discharge to 1e-8, no flow through the bed, velocity
profiles within 0.6 % of U_b, τ_b within 2.3 % of the BF run (laminar: 2.1 %
of the analytical value) and uniform to 0.2 %.

### T2 `sandSlideRelax`

A conical pit (`pit45`, `pit60`) or mound (`mound45`) with rim radius
R = 0.05 m and side slope β > φ_r (μ_s = tan φ_r = 0.625) in a height-field
bed on a 5 mm grid (`initialElevation` expression). The flow is not solved
(`solveFlow false`) and the prescribed bed shear is zero, so only the sand
slide changes the bed. Analytical relaxed state (volume conserved,
surroundings flat): R′ = (R² d/μ_s)^{1/3}, d′ = μ_s R′. The slide limits the
slope across every bed edge, so the relaxed cone on the Cartesian bed is
slightly square in plan and deeper than d′.

| Metric | Tolerance |
|---|---|
| `slope_max_final` (steepest bed edge) | ≤ 1.02 μ_s |
| `volume_err_rel` = max\|V(t) − V(0)\|/V_cone | ≤ 1e-6 |
| `depth_err_rel` = \|d_final/d′ − 1\| | ≤ 0.10 |
| `log_volume_error_max` (Exner volume error of the bed report) | ≤ 1e-12 m³ |
| `n_slide_not_converged` | 0 |
| `steps_to_repose`, `slope_max_initial` | REPORT |

All three variants reach 1.01 μ_s in the first step with the volume
conserved to 2e-11; the depth is within 2.8 % of d′ (7.3 % for `pit60`).

### T3 `humpMigration`

A quasi-2D channel (4 cells across, slip sides), h = 0.1 m, U_b = 0.45 m/s,
flat-bed θ₀ ≈ 0.09 (≈ 2θ_c0; d50 = 0.26 mm), with a Gaussian hump
A = 0.01 m, σ = 0.05 m at x = 0.4 m in a height-field bed (5 mm faces,
`initialElevation` expression). One run: 3 s of flow on the fixed bed
(`morphStartTime`), then 15 s of mobile bed with morphological factor 20,
Engelund–Fredsøe bed load and Koch–Flokstra gravitational transport (which
damps the grid-scale bed modes); a fixed strip at the inlet feeds sediment at
equilibrium. A flat-bed control run (`*_flat`) is subtracted. Fields and bed
are written every 0.5 s.

**Analytical celerity** (quasi-steady rigid-lid flow, constant friction
factor, small amplitude): θ = θ₀(H/h)², so ∂θ/∂z_b = 2θ₀/H and

  c = MF/(1 − n) · dq_b/dθ · 2θ₀/H,

with dq_b/dθ = 18.74 √(R g d³) [√θ − 0.7√θ_c + (θ − θ_c)/(2√θ)] for the
Engelund–Fredsøe law. The comparison is REPORT only: the estimate leaves out
the slope correction of the threshold and the phase lag of the shear.

| Metric | Tolerance |
|---|---|
| `log_volume_error_max` (Exner volume error of the bed report) | ≤ 1e-12 m³ |
| `spike_rel` (max second-difference residual / A) | ≤ 0.05 |
| `span_spread_rel` (max across-span range / A on the hump: columns with z_b − z_control > 0.1 A) | ≤ 0.05 |
| `flat_noise_rel` (std of z_b / A on the flat bed upstream, linear trend removed) | ≤ 0.02 |
| `celerity`, `celerity_analytic`, `celerity_ratio`, `theta0`, `amplitude_ratio_final`, `volume_change_rel`, `flat_std_rel`, `span_spread_wake_rel` | REPORT |

The span spread is tested on the hump only. Behind the steep front the flow
separates, and on the narrow slip-sided span (4 cells) the wake varies across
the span, so the scour trough behind the front carries small 3-D bed waves.
Their spread (`span_spread_wake_rel`) is transient and depends on the
start-up: 0.041–0.052 A in three start-ups, against 0.015–0.022 A on the hump.
The check was restricted to the hump after the ready-made case peaked at
0.052 A in the trough (t_m = 200 s, x = 0.79 m).

`volume_change_rel` (hump volume over 0.15–1.05 m minus the control) is not
a conservation measure here: the hump lowers the shear back to the inlet, so
more of the equilibrium feed deposits upstream than in the control
(+35 % of the hump volume after 300 s of morphological time). Conservation
of the bed update is `log_volume_error_max` (2e-22 m³). `mf20` passes; the
celerity is 1.29 times the estimate, and the hump keeps 72 % of its height.
The hump does not translate as a block: the crest moves only from 0.41 to
0.50 m while the lee face steepens to the angle of repose and advances from
0.45 to 0.70 m, and the stoss side lengthens (the nonlinear bed load moves the
crest faster than the toe, the Exner analogue of wave steepening).

## Body-fitted and IB-vs-body-fitted cases (new library)

These cases use the shared sediment core (`bedSediment`) with stock
`pimpleFoam` (body-fitted, `deformingBed` function object) or
`ibScourPimpleFoam` (immersed boundary). They are generated by scripts and
run by hand (commands below) or have their own `Allrun`; they are not part
of `run_tests.sh`.

| Case | Purpose | Script(s) |
|------|---------|-----------|
| `bedRegionMesh` | the in-memory plan-view bed mesh gives the same faMesh as a flat floor patch (serial and parallel) | `Test-bedRegionMesh` |
| `suspendedLoad` (deposition) | still water over a flat immersed bed, uniform C0 = 0.01 in 0.1 m of fluid, ws = 0.01 m/s, no entrainment: the bed rises by C0 H/(1 − n) = 1.667 mm (1.6664 at 60 s), the sediment budget (bed + suspended) holds to 6e-8 and the Exner volume error to 1e-21 (`tests/suspendedLoad/cases/depositionFlat`); the suspended-load option follows libsuspendedSediment (`constant/sedimentProperties`, `ibBedExchange { mode }`) | `ibScourPimpleFoam` |
| `suspendedLoad` (Rouse profile) | still column Z in [0.05, 1] with the prescribed parabolic diffusivity kappa Z (1 − Z), Vs = 0.2 (Rouse number 0.488), Cheng (1984) Robin condition at the reference level: equilibrium C against the Rouse closed form to 5e-4 of C*, bed flux equilibrium 2e-6 (`tests/suspendedLoad/cases/rouseColumn`, 1 min) | `suspendedSedimentFoam` |
| `suspendedLoad` (settling column, body-fitted bed) | the settling column with the `deformingBed` function object (no flow, no bed load, no mesh motion) receiving D − E of the bed condition: bed rise C0 H/(1 − n) = 0.3333 mm to 1e-3 (`tests/suspendedLoad/cases/settlingColumnBed`) | `suspendedSedimentFoam` + `deformingBed` |
| `suspendedLoad` (settling column) | the settling column with `suspendedSedimentFoam` (now in `solvers/`): still water, H 0.2 m, C0 1e-3, Stokes settling, deposition through the `sedimentBedFlux` bed condition: Kynch deposit ws C0 t, final deposit C0 H, budget 7e-15 (`tests/suspendedLoad/cases/settlingColumn`, 10 s) | `suspendedSedimentFoam` |
| `suspendedLoad` (HRN column, immersed bed) | the turbulent flume flow (Fuhrman flume, k-epsilon, ks 9.9 mm) over a plane immersed bed on a cell face, Liu (2014) sediment 1, HRN linear with the first-cell diffusivity: equilibrium E = D, wall-cell C equal to the two-point closed form C_b* (1 + P2)(1 − Pe1) to 1e-3; against a body-fitted run u* 0.01929 vs 0.01934 m/s and C within 1.5 % over the depth (`tests/suspendedLoad/cases/hrnColumnIB`, serial, about 8 min) | `ibScourPimpleFoam` |
| `deformingBedHump` | body-fitted hump migration: coupling, celerity, adaptation length | `make_case.py`, `check.py`, `check_adaptation.py` |
| `bodyFittedPier` | body-fitted Roulund pier (O-grid): frozen bed from an IB `bed.vtp` or live bed (layered mesh motion, adaptation length) | `make_case.py`, `to_layered.py` |
| `inclinedBed` | 2D channel with a trapezoidal pit (down-ramp, bottom, up-ramp): bed shear stress of the IB vs a body-fitted reference on slopes | `make_case.py`, `compare.py`, `profiles.py` |

### `bodyFittedPier`

```bash
make_case.py setup --run-dir DIR --inflow REF --mode live \
    --dr 2.2e-3 --dz 1e-3 --n-r 30 --nz 45 --end-time 155 \
    --morph-start 5 --morph-factor 2 [--adaptation-length 2.6e-3]
blockMesh -case DIR
make_case.py init --run-dir DIR      # fields + motion zones
```

Live mode moves the mesh with `displacementLayeredMotion`: every vertical
column stretches linearly between its bed point and the lid. A Laplacian
displacement solve with cell-to-point interpolation smears the bed
displacement over neighbouring columns and collapsed the first cells on
deposition mounds. `to_layered.py CASE TIME` converts an existing
decomposed case (Laplacian, `pointDisplacementz`) to layered motion.
With layered motion and L = 10 d50 the coarse case (530 k cells, 48 ranks,
about 6.5 h) reproduces the measured upstream scour depth of Roulund et al.
(2005) within about 0.05 D up to t = 5 min.

### `inclinedBed`

```bash
make_case.py setup --variant ib|bf --slope-deg 30 --run-dir DIR \
    [--dx 2e-3 --dz 1e-3 --pit-depth 0.03]
blockMesh -case DIR
make_case.py init --run-dir DIR
decomposePar -case DIR                      # parallel runs
make_case.py decompose-zb --run-dir DIR     # IB only, see below
compare.py --ib IBRUN --bf BFRUN --out fig.png
profiles.py --ib IBRUN --bf BFRUN --time T --x 0.90 0.94
```

The bed is frozen (morphology never starts). Both variants write
`postProcessing/bed/<t>/bed.vtp`; `compare.py` tabulates the mean tau_x per
segment and the IB/BF ratio. Notes:

* `decomposePar` does not decompose `0/finite-area/zb` (the IB faMesh is
  built at run time), so a parallel IB run needs `decompose-zb`; frozen IB
  beds do not write `zb` at write times, so a restart needs it copied into
  the restart time.
* Findings (Oct 2026): the IB is accurate on the flat sections; on the
  up-ramp it gives 0.40 (15 deg) and 0.11 (30 deg) of the body-fitted
  shear stress. Near-wall IB cells on a slope (y_c -> 0) carry a huge
  omega that laminarises their neighbours; `wallModel { yTurbMin 0.5; }`
  removes this (0.49 and 0.23). The remaining deficit does not change with
  the sampling distance or with grid refinement (stair-stepped blocked
  faces).
* `compare_stair.py`: a body-fitted run on the IB staircase itself
  (`topoSet` fieldToCell Gamma > -0.5 on the IB mesh, `subsetMesh -patch
  bed`, wall law on all stair faces) reproduces the IB results (up-ramp
  0.38 of the smooth body-fitted shear, flat downstream identical to the
  IB): the deficit is caused by the stair-stepped blocked faces (form drag
  of the risers), not by the wall model.

## Results in ParaView

Every run directory keeps all of its time directories (nothing is pruned)
and contains an empty `<run>.foam` file. `common/make_view.py` builds
`<run>/view/` with links to `constant`, `system` and the time directories
that contain flow fields, plus `view.foam`; `ibs_view`
(common/runFunctions.sh) calls it at the end of every run. To rebuild views
by hand:

```bash
python3 tests/common/make_view.py <run> [<run> ...]
```

The height-field bed is written as `postProcessing/bed/<time>/bed.vtp`
(z_b, bed shear stress, bed-load flux, active and fixed flags; merged in
parallel) and as the area field `<time>/finite-area/zb` (read by the
checks). Parallel runs are reconstructed (`reconstructPar`, all times).
