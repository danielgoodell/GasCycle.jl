# NPSS Cross-Validation Plan: BRU3.mdl ↔ GasCycle

Written 2026-06-10, reconstructing the validation campaign Daniel intended
before the previous session ended (it was never documented; several
prerequisites were quietly built: the `η_type=:isentropic` option matching
NPSS `effDes`, `set_boundary!` for the HeatSinkHx coolant side, `Oil.fpt` /
`H2O.fpt` readers, and the item-6 cold-end elements).

## Goal

Run GasCycle against the **identical** NPSS model
(`reference/BRU3.mdl`) and either:

1. match every station quantity to within solver precision (both solvers at
   ~1e-6 relative tolerance), or
2. attribute every digit of remaining difference to a specific, named cause.

The current state (README, `examples/bru_10kw.jl`) is ~1.5 % on temperatures
and ~7 % on heater input with *hypothesized* causes. That is not validation;
it is consistency. This campaign replaces hypotheses with measurements.

## Why the current comparison cannot be tight: known confounders

Each of these alone is larger than solver precision. They must be removed or
measured one at a time.

1. **Different fluid model.** BRU3.mdl runs `fs.comp = "CEAT"` (a
   CEA-generated table we do not have); GasCycle runs `HeXe84.fpt`. Estimated
   ~14 K at comp outlet, ~6 % net power. **Removal:** rerun NPSS with the
   commented-out `fs.comp = "HeXe84"` line (one-line change — the model was
   originally built for it), or obtain the CEAT table file.
2. **Different efficiency semantics.** NPSS `effDes` is adiabatic
   (isentropic); `bru_10kw.jl` uses polytropic. **Removal:** the
   `η_type=:isentropic` option now exists — the replica must use it.
3. **Different solver constraint sets.** This one was previously
   unrecognized. BRU3.mdl is *not* a self-consistent closed loop with fixed
   TIT and pressure closure (what `bru_10kw.jl` solves). Reading its
   independents/dependents (`autoSetup` flags + `autoSolverSetup()`):
   - `HotStart` tears the loop at station 6/7 with **fixed Tt = 1701 °R and
     W = 1.32 lb/s**, only its Pt floating (`ind_Pt.autoSetup = TRUE`).
   - `HeatSourceHX.Tout` is an **independent** (`ind_HeatSource.autoSetup =
     TRUE`): TIT floats so the turbine outlet matches the 1701 °R tear.
     2060 °R is only its initial value.
   - `Turb.PRdes` is an **independent** (`ind_PRdes.autoSetup = TRUE`):
     turbine PR floats to balance the shaft against the `HPX` load at fixed
     36 000 rpm. 1.75 is only its initial value.
   - `Start` pins the comp inlet at 540 °R / 23.7 psia, W fixed at 1.32
     (its `ind_StartW` exists but is not auto-enabled).
   - The oil-loop independents (`Dow200Start.Wflow`, `HeatSinkHx.effect`,
     radiator area) exist but are **not** auto-enabled; `Dow200Start.dep_Tt`
     is disabled — the oil enters the HeatSinkHx at fixed 527 °R.
   So NPSS solves: unknowns {TIT, PR_t, Pt₆} s.t. {Tt_turb_out = 1701 °R,
   shaft power balance with HPX extraction, Pt tear closure}. The replica
   must mirror this exact constraint set (or both models must be reposed
   into the same physically-clean closed loop — preferable long-term, since
   the .mdl header itself flags the 1701 °R / 1.32-vs-1.28 lb/s
   inconsistencies).
4. **Interpolation scheme.** `FPTFluid` uses bilinear `Gridded(Linear())`
   (the README previously claimed bicubic — fixed). NPSS's FPT lookup scheme
   must be confirmed (likely linear on the same NEO grid, possibly with
   different unit constants). Same table + same scheme ⇒ this term should be
   ~0; measure it pointwise first (Rung 0 below).
   **RESOLVED on the GasCycle side (2026-06-10, see RESULTS.md):** linear-in-P
   interpolation of s ∝ −R·ln P was distorting ∂s/∂P by 14 % mid-cell and
   caused +14 °R at the compressor outlet — the bulk of the historical
   "fluid model" gap. Fixed by pressure-detrended interpolation
   (`s_interp=:log_pressure`, default); `s_interp=:linear` retained for
   NPSS-compat comparison runs.
5. **Bleed path details.** 2 % comp-exit bleed (`fracBldP = 1`), forced to
   559 °R in the bearing housing, reinjected at the turbine **inlet**
   (`Bld.Pfract = 0`). GasCycle replica: Splitter → fixed-TtExit HeatSource →
   Mixer (enthalpy-balance, Pt = min). NPSS's bleed-port mixing rule (does it
   charge a mixing pressure loss? mass-average Pt?) must be checked against
   our Mixer once outputs are in hand.
6. **Bleed vs tear inconsistency in BRU3.mdl itself.** The 1701 °R HotStart
   tear value equals the *no-bleed* isolated turbine output (GasCycle
   reproduces it to 0.2 °R), but the model also injects the 2 % bleed
   (cooled to 559 °R) at the turbine inlet, which lowers the self-consistent
   turbine outlet to ~1673 °R. The .mdl header flags exactly this ("Current
   issue seem the 1701 R, look at Interstage Bleed Port"). The NPSS run
   listing will show which constraint wins in NPSS; comparing no-bleed runs
   of both models sidesteps it entirely.
7. **HPX units ambiguity.** `HPX = 10.9/.92 + 1.572` with the comment
   "horsepower" but a value that reads like kW (10.9 ≈ alternator kW /
   0.92 η_alt + parasitic kW). If HPX is in hp the shaft extracts 10.0 kW;
   if kW, 13.4 kW. `bru_10kw.jl` currently assumes 13.4 kW. The NPSS output
   listing settles this immediately (NPSS shaft power is natively hp; the
   .mdl's own `pwrkW = 0.7457*pwr` suggests hp).

## The precision ladder

Compare bottom-up; do not debug a rung until the one below is at solver/
table precision. Each rung has a harness in `validation/`.

- **Rung 0 — property lookups.** Same FPT file on both sides
  (`HeXe84.fpt`, `Oil.fpt`): h, s, cp, γ, ρ at the exact station (Tt, Pt)
  points, GasCycle `FPTFluid` vs NPSS `FlowStation` prints. Confirms parsing,
  unit constants (Rankine/BTU!), reference states, and interpolation scheme.
  Target: exact to output-print precision, or a named interpolation
  difference.
- **Rung 1 — single elements in isolation.** BRU3.mdl already contains the
  exact harness (its diagnostic section, lines ~288–311): turbine at
  (2060 °R, 43.2 psia, W = 1.31, PR = 43.2/24.69, effDes = 0.87); recup at
  (737 °R/45.03 + 1701 °R/24.68); HeatSinkHx at (786 °R gas + 527 °R oil).
  `validation/bru3_isolation.jl` mirrors it. Compare outlet Tt/Pt digit by
  digit. Isolates element physics + property model with no solver coupling.
- **Rung 2 — full model, identical constraint set.** Replica solving exactly
  the NPSS unknowns/residuals from confounder 3. Station-by-station diff
  table (Tt, Pt, W, ht at stations 0–9), component powers, Q's.
- **Rung 3 — factor attribution.** For whatever residual difference remains
  (and for the historical record): toggle one factor at a time —
  η semantics, bleed on/off, Mixer Pt rule, cold-end modeled vs pinned,
  fluid backend — and tabulate the per-station effect of each. This is the
  "understand EXACTLY where the difference comes from" deliverable even if
  rung-2 convergence to solver precision is blocked.

## NPSS-side artifacts needed (blocked on NPSS access / collaborator)

1. Full design-point output listing of BRU3.mdl as-is (`ncpView` station
   table: Tt, Pt, W, ht, s at stations 0–13; comp/turb `pwr`; HX Q's;
   converged independents TIT/PR_t/Pt₆; solver tolerance used).
   **Still the top request** — it pins the last ~1 °R / 2 % residuals.
2. The four `cout` lines from the isolation diagnostic section.
3. ~~A rerun with `fs.comp = "HeXe84"` or the CEAT FPT file.~~ **RESOLVED
   2026-06-10:** CEAT.fpt received (`reference/CEAT.fpt`) — it is a live
   CEA passthrough, not a table; for He-Xe it is exactly ideal monatomic
   with M = 83.328 (w_He = 0.0181), replicated to machine precision by
   `IdealGasFluid(M_molar = 83.328)`. See RESULTS.md.
4. Confirmation of NPSS's FPT interpolation scheme (or just rely on Rung 0
   pointwise comparison). Less urgent: BRU3.mdl never interpolates a table
   (CEAT computes live), so this only matters for future HeXe84.fpt runs.
5. The FPT-generation Python (already wanted for roadmap item 4).
6. ~~`solver.bad` diagnostic file.~~ **RESOLVED 2026-06-10:**
   `reference/solver.bad` received — setup dump only (no converged values),
   but it exposed the active constraint set and proved the as-dumped 5×5
   system is structurally unsatisfiable (no turbine-PR independent), i.e.
   it documents a failed mid-debug configuration. See RESULTS.md.
7. (New) If the collaborator can rerun: a converged listing **with the
   bleed actually flowing** and HPX units confirmed — the rung-2 replica
   shows the historical comments correspond to a no-bleed, HPX-as-kW run.

## What is already runnable today (no NPSS artifacts)

- Rung 1 GasCycle side: `validation/bru3_isolation.jl` (NPSS columns marked
  pending).
- The η-semantics factor of Rung 3: isentropic vs polytropic comp/turb
  outlet on HeXe84.fpt at the BRU state points.
- The fluid-backend factor: HeXe84.fpt vs ideal-gas M = 83.8 at the same
  points.
- Known anchors from the .mdl comments for coarse checks: station 1 =
  737 °R / 45.03 psia, station 6 = 1701 °R / 24.69 psia, station 8 =
  786 °R, heater Q ≈ 33.1 kW.

## Success criteria

- Rung 0: agreement to print precision (or a named, quantified
  interpolation difference).
- Rung 1: outlet temperatures within 0.1 °R with the same fluid table;
  every larger deviation traced to a named cause.
- Rung 2: all station Tt/Pt within combined solver tolerances with the same
  fluid table; net shaft power within 0.1 %.
- Rung 3: a table in `validation/RESULTS.md` attributing 100 % of the
  historical 1.5 %/7 % gaps to measured factors.
