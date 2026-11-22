# libsuspendedSediment: bed exchange of suspended sediment

OpenFOAM v2606, GPL-3.0-or-later. Used by `solvers/suspendedSedimentFoam`
(body-fitted bed) and by `libibScour` (`ibSuspendedLoad`, immersed bed);
cases in `tests/suspendedLoad/cases/`.

Build: `./Allwmake` (into `$FOAM_USER_LIBBIN/libsuspendedSediment.so`).

## Contents

| Item | Type | File |
|---|---|---|
| `sedimentBedFlux` | `mixed` patch condition for C | `derivedFvPatchFields/sedimentBedFlux/` |
| `sedimentBedExchange` | function object: E − D per bed face, budget | `functionObjects/sedimentBedExchange/` |
| `settlingVelocityModel` | run-time selectable: `constant`, `stokes`, `vanRijn1984`, `soulsby` | `sedimentModels/settlingVelocityModel.*` |
| `referenceConcentrationModel` | run-time selectable: `zysermanFredsoe` (default), `vanRijn1984`, `constant` | `sedimentModels/referenceConcentrationModel.*` |

The Zyserman-Fredsøe constants reproduce the C_b* of Liu (2014, Table 3).

## sedimentBedFlux

Net upward sediment flux through the bed = entrainment − deposition:

    Gamma_f snGrad(C) − w_n C_f = E − beta C_f

with `Gamma_f` the bed-face value of the diffusivity `DC` of the C equation,
`w_n = ws·n` the settling velocity through the patch and `C_f` the patch
value. The left side is exactly what the finite-volume C equation carries
through the face (diffusion in, drift flux `phi + ws·Sf` out), so the
discrete bed flux is `E − beta C_f`. As an OpenFOAM `mixed` condition
with `K = beta − w_n ≥ 0`:

    valueFraction = K/(K + Gamma_f deltaCoeffs),  refValue = 0,  refGradient = E/Gamma_f

which is the Robin family of Cheng (1984): his Eq. 3,
`Q_0 = −[k ∂c/∂z + w c]_(z=0) = −beta (c − c_*) = E − D` with `E = beta c_*`
and `D = beta c_0` (J. Geophys. Res. 89(C5), p. 8209), is the left side
above with the outward (downward) bed normal. For `beta = w_n` (deposition
= settling flux) `valueFraction = 0` and the condition is the fixed
gradient of Liu (2014, 2018).

| `mode` | E | beta | intended patch location |
|---|---|---|---|
| `referenceLevel` | `w_n C_b*` | `w_n` | mesh bottom at y = delta (Liu 2014 Eq. 10); at the wall it is the "naive" condition |
| `HRN` | `w_n C_b* (1 + P2)` (`form linear`, Eq. 18) or `exp(P2)` (`exponential`, Eq. 17), `P2 = w_n delta/eps_s` | `w_n` | wall, wall functions (Liu 2014) |
| `LRN` | `w_n C_b* exp(∫_0^delta w_n/DC dy)`; `integral computed` (model nut along the cell column, default) or `table` (I1(ks+, delta+) of Liu 2018 Eq. 10, `I1Model vanDriestCC` (default), `SST` or `file`) | `w_n` | wall, y+ < 1 (Liu 2018) |
| `Robin` | `beta C_b*`, `beta` given | `beta` | reference level (GITT benchmarks) |

`diffusivity` (all modes) sets the bed-face eps_s (= Gamma_f):
- `firstCell` (default for HRN) uses the
  wall-cell value, the constant eps_s of Liu (2014); it reproduces Liu's
  reference-height equivalence to 0.01 % (Z = 0.06) and 3 % (Z = 0.6);
- `patch` (default for the other modes) keeps the patch value of DC
  (wall-function nut_w for HRN, nu/ScC for LRN); with HRN it gives negative
  C in the wall cell when `w_n y1/Gamma_f > 1` (Liu 2014 sediment 2);
- **HRN with `firstCell` stops with a FatalError when the first-cell
  Péclet number `sigma_c w_n/(kappa u*)` is ≥ 1 on any bed face** (reduced
  over the processors; faces with u* = 0 are skipped). With wall functions
  Gamma_1 = kappa u* y1/sigma_c, so `w_n y1/Gamma_1` equals this Rouse
  number for every y1, and for values ≥ 1 the two-point wall flux drives C
  negative (van Rijn 1981 case, Z' = 1.19). The message names
  `diffusivity logLaw` as the remedy. The other options and modes
  have no such stop;
- `logLaw` is **an extension, not part of Liu (2014)**:
  `Gamma_f = w_n y1/(1 − (delta/y1)^Z')`, `Z' = sigma_c w_n/(kappa u*)`,
  which makes the two-point wall flux exact for the equilibrium profile of a
  log layer.

The solver applies the value (`faceDiffusivity()`) before assembling the C
equation. Findings of the verification runs:
- HRN is **not mesh independent inside the guideline 0.2 k_s < y1 < 0.1 H**
  of Liu (2014) when Z is of order 0.5 or more: C(y) scales with y1^Z
  (spread over y1 = 2-6 mm: 7.9 % for Z = 0.06, 77 % for Z = 0.6 with
  `firstCell`; 1.4 % and 19 % with `logLaw`; turbulent flume of Fuhrman et al.);
- with `firstCell` the first-cell Peclet number `w_n y1/Gamma_1` equals the
  Rouse number `Z' = sigma_c w_s/(kappa u*)` in a log layer for every y1;
  for `Z' > 1` the two-point wall flux drives C negative (van Rijn 1981 case,
  Z' = 1.2); `logLaw` has no such limit;
- the LRN condition needs a near-wall mesh with `w_s Delta/nu` of about 0.1
  or less near the wall (the equilibrium profile decays over nu/w_s; low-Reynolds-number flume:
  first cell 0.12, grading 1.03, C(delta)/C_b* = 0.98).

**LRN tables.** `integral table` uses `I1(ks+, delta+) = ∫_0^delta+ dy+/(1 +
nut+)`, `ks+ = ks u*/nu` (entry `ks`), from a built-in table
(`derivedFvPatchFields/sedimentBedFlux/I1Tables.H`) or a file:
- `I1Model vanDriestCC` (**default**): van
  Driest mixing length with the Cebeci & Chang (1978) roughness shift
  (Eqs. 8–9, valid for 4.535 < ks+ < 2000); for ks+ ≤ 4.535 it is the
  smooth van Driest table exactly;
- `I1Model SST` (option): k-omega SST columns of OpenFOAM, smooth
  (omegaWallFunction) and rough (omega_w = u*² S_R/nu). **Caveats:** its
  smooth row misses the smooth van Driest table by 29 %: the low-Re SST near-wall state with omegaWallFunction
  is case dependent (nu_t+(15) = 0.78 in the table column, 2.5 in
  flume columns with the same wall-unit grading); the Wilcox
  (2006) S_R constants of the rough rows were **taken from the SU2 source
  and are not yet checked against Wilcox (2006)**;
- `I1Model file`: sub-dictionary `I1Table` (an interpolation2DTable;
  chosen automatically when `I1Table` is given).

The tables are **a recomputation by the authors**, not the
rough-wall table of Liu et al. (2013, IAHR), which is not available: a
stated limitation. Liu (2018) writes the integral with nu + nut and no
sigma_c; the `computed` mode uses DC = nu/ScC + nut/sigma_c.

u* is `sqrt(nuEff_f |snGrad(U)_t|)` on the patch (wall shear of the flow
solution). It sets the HRN P2, the logLaw Z', the LRN wall units and the
firstCell check. C_b* and the reference level delta come from
`referenceConcentrationModel` in `constant/sedimentProperties`, with the
shear velocity of its Shields number (next section).

## Reference concentration: which shear velocity

Entry `shieldsNumber` in `constant/sedimentProperties`:

| Model | `skinFriction` (**default**) | `total` |
|---|---|---|
| `zysermanFredsoe` | u*' of **Zyserman & Fredsøe (1994) Eqs. 9–10**: V/u*' = 6 + 2.5 ln(D'/k_N), u*' = (g D' I)^½, k_N = 2.5 d50 (entry `kN`), as in their procedure (step 1, p. 1038); with g D I = u*² this is V = u*' [6 + 2.5 ln(D u*'²/(u*² k_N))], solved per face | the wall-shear u* (Liu 2014, 2018) |
| `vanRijn1984` | grain shear velocity of **van Rijn (1984a,b)**: u*' = g^½ V/C', C' = 18 log10(12 D/k_N), k_N = 3 d90 (entry `d90`, or `kN`), R_b = D | the wall-shear u* |
| `constant` | — (no shear velocity) | — |

V is the depth-mean tangential velocity and D the depth of the cell column
above each bed face (bed-normal walk through hexahedral or prismatic cells
up to the top boundary of the mesh, built once; the column must lie on
one processor, otherwise the run stops with a message). In both models
u*' ≤ u* (skin friction is part of the total; a bed no rougher than k_N
gives u*' = u*). Use `total` to reproduce Liu (2014, 2018) or on a smooth
flat bed where the wall shear is the skin friction. The library value
agrees with an independent Python implementation to 10 digits (turbulent
flume, sediment 2 of Liu 2014: u*' = 0.012080 m/s against u* =
0.019344 m/s, C_b* 0.00306 against 0.0271).

Debug output (u*, u*', C_b* per update): `DebugSwitches { sedimentBedFlux 1; }`
in `system/controlDict`.

Example (HRN):

    bottom
    {
        type        sedimentBedFlux;
        mode        HRN;
        diffusivity firstCell;      // default for HRN
        form        linear;
        value       uniform 0;
    }

Needs the registered fields `DC` (diffusivity) and `ws` (uniform settling
velocity) and the dictionary `sedimentProperties` of suspendedSedimentFoam.
Accessors `E()`, `D()`, `Cstar()`, `uStar()` serve the function object and an
Exner solver (ibScourFoam interface).

## sedimentBedExchange

    bedExchange
    {
        type        sedimentBedExchange;
        libs        (suspendedSediment);
        writeFields yes;      // bedExchange (E − D) and bedDeposit (∫(D − E) dt) on the bed faces
    }

Writes per time step `postProcessing/<name>/<t>/sedimentBedExchange.dat`:
∫E dA, ∫D dA, ∫(E − D) dA, the same from the face fluxes (check), outflow
through the other patches, V_C = ∫C dV, the cumulative deposit and the
budget error Σ (dV_C/dt)ⁿ Δtⁿ + deposit + other. The storage term uses the
time scheme of the C equation (`ddtSchemes` entry for `ddt(C)`): Euler (and
steadyState) `(Vⁿ − Vⁿ⁻¹)/Δt`; backward with the variable-step
coefficients of OpenFOAM's `backwardDdtScheme`, Euler on the first step of
a run (no old-old field yet, as in the scheme); other schemes fall back to
Euler with a warning. The budget therefore closes to solver tolerance with
both schemes (Rouse column with backward: ≤ 3e-13 of the suspended volume).
