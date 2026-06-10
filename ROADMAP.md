# GasCycle.jl Roadmap

Status assessment as of 2026-06-09. Supersedes the "Suggested Next Steps" in
REVIEW.md (several of those items have since been addressed: NonlinearSolve
migration, polytropic component models, partial turbine-PR display fix).

## Where the project stands

Design-point analysis is solid and validated: thermo layer (ideal-gas + FPT
backends), full element library, flow network with branch support,
NonlinearSolve-based solver, ForwardDiff gradients verified against finite
differences, BRU 10 kW validation case. Full test suite passes.

Off-design is scaffolded but **not yet functional** — that is the main gap.

## Ordered next steps

### 1. Fix the off-design map residual formulation  ✅ DONE (87239c0)
The High-severity finding from REVIEW.md, plus a related degeneracy:

- `Compressor`/`Turbine` set `Wc_map = Wc_act` in `compute!`, then the
  residual computes `Wc_map - Wc_act` — identically zero.
- Worse: `indep_vars` for Compressor and Turbine both return `N_shaft`,
  while `Shaft` also returns `N` — three copies of the same unknown,
  masked by the two zero residuals.

Fix (NPSS-style formulation):
- The map flow coordinate `Wc_map` becomes the solver independent variable
  for each off-design turbomachine (seeded from actual corrected flow on
  the first pass).
- `compute!` queries the map at `(Nc, Wc_map)` — not at the actual flow —
  and applies the resulting PR and η.
- Residual: `(Wc_map - Wc_act) / Wc_map = 0` (flow continuity onto the map).
- Shaft keeps `N` as its independent with the power-balance residual.
- Add an end-to-end off-design test: design-point reproduction (off-design
  solve at design boundary conditions must recover the design N, PR, η)
  plus a perturbed-TIT case.

### 2. Two quick wins  ✅ DONE (a1ede04)
- `Pkg.test()` works: removed `Pkg.activate` from test files; Test and
  ForwardDiff declared via `[extras]`/`[targets]`.
- `pressure_ratio(::Turbine)`/`(::Compressor)` accessors added and used by
  all four examples (fixes stale PR in simple_brayton and bru_tit_sweep).

### 3. Off-design TIT sweep validation (plan Phase 7)  ✅ DONE (a661a48)
- Back-edge (Tt, Pt) states folded into the off-design Newton vector
  (replacing stale fixed-point seeding); unknowns normalized (raw scaling
  gave cond(J) ~1e6); TrustRegion for robustness across map-cell kinks.
- `Shaft` gained `P_load` (generator extraction) for the off-design power
  balance.
- Closed-loop BRU-like design-point reproduction + constant-speed TIT sweep
  tests; `examples/bru_tit_sweep_offdesign.jl` sweeps the FPT-fluid BRU
  cycle 2060→1236 °R at 36 krpm (21/21 converged, design power reproduced,
  self-sustain threshold near 1330 °R).

### 3b. NPSS cross-validation campaign (BRU3.mdl)  ✅ DONE (2026-06-10)
Goal achieved same day: with `reference/HeXe.out` (NPSS 3.3 on HeXe84.fpt),
GasCycle matches every station/power within NPSS print precision, TIT
included (each code solving its own shaft balance; Δ = 0.025 °R). All
historical gaps attributed: GasCycle s-interpolation artifact (fixed;
NPSS confirmed to share it — `s_interp=:linear` reproduces it), isentropic
effDes, CEAT = ideal monatomic M = 83.328, NPSS bleed = workless bypass
(dhb/dh = 0), HPX kW-tagged (18 hp internal), oil cp = Oil.fpt's 0.8.
Ledger: `validation/RESULTS.md`. Open: sink-HX effect definition
(behaves as ε = 0.9405 vs 0.946 setting; gas loop unaffected).
Validate GasCycle against the identical NPSS model to solver precision, or
attribute every difference exactly. Full plan in `validation/PLAN.md` —
precision ladder (property lookups → element isolation → full loop with
NPSS's actual constraint set → factor attribution), the six known
confounders (fluid table CEAT vs HeXe84, effDes isentropic semantics,
BRU3.mdl's tear-based solver setup with floating TIT and turbine PR,
interpolation, bleed mixing rules, HPX units), and the NPSS-side artifacts
needed from the collaborator. Rung-1 GasCycle harness:
`validation/bru3_isolation.jl`.

### 4. Direct noble-gas property module  (waiting on collaborator's Python)
Replace the FPT-table dependency with a native property backend. ON HOLD
(2026-06-09) until Daniel obtains the collaborator's FPT-generation Python —
plan is to three-way compare code vs. papers vs. FPT output and resolve any
discrepancies before implementing.

Design note: implement as a general `NobleGasMixture(gas1, gas2, x1)`, not
hard-coded He-Xe. The El-Genk paper's correlations cover all 5 noble gases
(He, Ne, Ar, Kr, Xe) and their 10 binary mixtures with per-gas constants in
its Table 1 — so Ar and Kr (the cheap stand-in test gases for He-Xe systems,
including He-Ar and He-Kr mixtures) come for free from the same
implementation.

Method papers are in `reference/` (added 2026-06-09) and are sufficient to
implement directly; the collaborator's Python adds cross-checking:

- **Johnson, NASA/CR-2006-214394** — dilute transport (μ, k): Hirschfelder
  first-order binary Chapman-Enskog with LJ parameters (He σ=2.576 Å,
  ε/k=10.22 K; Xe σ=4.055 Å, ε/k=229 K), Ω(2,2)*/A*/B* tables, Singh
  third-order conductivity correction. Tables 4–6 are an exact test oracle
  (μ and k at M = 20.183/39.94/83.8, 400–1200 K); Table 7 has Prandtl data.
- **Tournier/El-Genk/Gallo, AIAA 2006-4154** — real-gas thermo + dense
  transport: virial EOS Z = 1 + Bρ̂ + Cρ̂² with corresponding-states
  correlations for B(θ), C(θ) and He-specific B (Eqs. 8–11), mixture
  combining rules (Eqs. 16–22), enthalpy/Cp/Cv from B(T), C(T) derivatives
  (Eqs. 12–15), dilute-mixture μ/k via Sutherland-Wassiljewa with
  tabulated coefficients (Eqs. 30–32, 36, Tables 1–3), dense-gas excess
  corrections Ψμ(ρr), Ψλ(ρr) (Eqs. 4, 23–28, 33). This matches the
  real-gas behavior observed in HeXe84.fpt.

Validation oracles: Johnson Tables 4–7, El-Genk spot values (e.g. Xe at
2 MPa/400 K: Cp +10.7%, γ=1.786), HeXe84.fpt itself, and Taylor-1988
experimental Prandtl numbers.

Decided 2026-06-09, superseding the earlier "FPT AD" item — rationale:

- Checked HeXe84.fpt against ideal monatomic gas: within ~0.5% at cycle
  conditions (300–1150 K, ~100–300 kPa), but real-gas corrections are
  genuinely present near the xenon-critical corner (at 260 K / 1.5 MPa:
  Cp −9%, density +4.8% vs ideal). The generator includes a virial-style
  EOS, not just kinetic theory.
- A direct module gives any mixture ratio without regenerating files, no
  interpolation error, smooth functions, and free ForwardDiff support if
  written generically — which makes AD-through-tables mostly moot.
- The transport-property part (μ, k, Pr) is the piece that can't be
  trivially rederived, and will be needed for real heat-exchanger sizing
  (NTU from UA) and loss models. The `FluidProperties` interface doesn't
  expose transport yet — add that alongside.

Keep the FPT reader regardless: reading the literal file the collaborator
feeds NPSS is the cleanest apples-to-apples cross-validation available.

### 5. Additional working fluids: air, N₂, CO₂
Real systems are commonly first tested with cheap fluids: air (low-temp
checkout, motoring), N₂, or CO₂, before committing the expensive xenon
inventory. Ar/Kr are covered by item 4; air/N₂/CO₂ are di-/triatomic with
temperature-dependent Cp, so monatomic kinetic theory does not apply.
Implement a `ThermallyPerfectGas` backend:
- Cp(T) from NASA 7- or 9-coefficient polynomials (h, s by closed-form
  integration — keeps AD support exact like IdealGasFluid).
- Ideal-gas EOS is adequate for checkout/motoring conditions; revisit only
  if someone needs sCO₂-style operation near the CO₂ critical point.
- Transport via Sutherland law or polynomial fits per gas.
- Air as a fixed-composition pseudo-species (standard dry air), N₂ and CO₂
  as pure species.

## Capability expansion (added 2026-06-09)

Priority driver: mission needs (SR-1 Freedom) put **off-design fidelity,
inventory control, and slow transients** at the top of the physics track.
The convenience batch is cheap and can be interleaved anytime.

### 6. Cold-end physics: Cooler / Radiator element  ✅ DONE (2026-06-10)
There is currently no heat-rejection element — loop closure leans on
`set_state!` pinning the compressor inlet. Add a `Cooler` (ground test:
water-cooled HX with ε or UA) and `Radiator` (space: T⁴ rejection with area
and sink temperature) so compressor inlet temperature responds to the
operating point. Prerequisite for meaningful off-design and transients.

Implemented as four pieces (see `test/test_cooler.jl` for the end-to-end
closed-loop demos):
- `ConstantPropertyLiquid` coolant backend (constant cp/ρ, closed-form
  inversions, full AD); reads the collaborator's function-style FPT files
  (`H2O.fpt`, `Oil.fpt` — constant-return Cp/rho in BTU/lbm·R and lbm/ft³,
  including Oil.fpt's one-level `Cp → Cpt` delegation).
- Cooler = `HeatExchanger` with the new `UA` keyword (`mode=:UA`):
  counter-flow ε-NTU each pass, so effectiveness responds to off-design
  flows/properties.  Coolant side enters via `set_boundary!`.  This is also
  the machinery item 8 needs (only the UA(ṁ) scaling law remains).
- `Radiator` element: segmented σεA(T⁴−Tsink⁴) integration (Heun marching);
  `mode=:fixed_TtExit` sizes A at design, `mode=:fixed_area` is the
  off-design response mode.  Two-way flux (warms toward sink if entering
  below it).
- Loop closure through the cold end: `connect_port!(...; back_edge=true)`
  marks cooler/radiator-outlet → compressor-inlet as a back-edge, so the
  compressor inlet (Tt, Pt) joins the Newton unknowns and `set_state!` only
  supplies W, fluid, and the initial guess.  Works in design mode when the
  turbine uses `:pressure_closure` and in map-based off-design (where the
  maps give pressure closure a free variable); a design loop with every PR
  fixed is structurally degenerate and still needs the pinned inlet.

### 7. Inventory control and per-component volume bookkeeping
Real CBC power control is charge-pressure (inventory) control. Plan:
- Make loop inventory (total gas mass/moles) a model quantity: assign a
  volume to every component; where components are currently joined
  directly, connecting ducts may need to be added as the volume holders.
- Steady-state inventory mode: loop mass W becomes a solver independent
  with a pressure-closure residual (replaces the fixed-seed-Pt assumption;
  resolves the dangling-pressure caveat from the off-design TIT sweep).
- Mass accounting in components with large T/P gradients is the hard part:
  m = ∫ρ dV over a strong axial temperature gradient (recuperator!) is not
  well-approximated by inlet/outlet averages. Use segmented volumes or an
  analytic mean density per component (e.g., ideal-gas with linear-T duct
  has closed-form m = PV/(R·T_lm) using log-mean temperature); validate
  segment-count convergence.
- Sets up transients: volume + mass state per component is exactly the
  capacitance structure the transient extension needs.

### 8. Off-design recuperator effectiveness (ε-NTU from scaled UA)
ε is currently a fixed parameter; physically it varies as flows and
properties deviate from design. Compute ε from NTU with UA scaled from the
design point (UA ∝ ṁ^0.8-type scaling or Colburn-based). Fold in the
enthalpy-based energy balance from the backlog. Highest realism-per-effort
in the physics track since recuperator performance dominates cycle η.

### 9. Transients — shaft dynamics first, then thermal
The original mission of this architecture (reactor coupling). Stage it:
- Shaft speed dynamics: I·ω·dω/dt = unbalanced shaft power. The power-
  balance machinery already computes the imbalance; startup, load-step,
  and loss-of-load events become solvable with just shaft inertia.
- Slow thermal transients: thermal capacitance (metal + gas) per component
  using the volumes from item 7, wrapped with DifferentialEquations.jl.
  Reactor side couples here (point kinetics or supplied Q(t)).

### 10. Convenience batch  ✅ DONE (ef8c91d…354fa66, 2026-06-10)
- `summary(sol)`: NPSS-style station table in physical flow order
  (back-edges walked so the recuperator hot side reads in cycle position,
  bleed branches after) + component table + cycle totals.  `stations(sol)`
  returns the ordered station data.
- NPSS `.map` reader: `read_npss_map` (generic NEO Table parser) +
  `to_performance_map` (Rline speed lines → rectangular (Nc, Wc) map,
  exact along each tabulated line).  Validated against a synthetic
  format-faithful fixture only — **cross-check against a real NPSS map
  when one is available**; maps operated near choke (vertical speed-line
  segments) eventually want an Rline-native solver coordinate instead of
  the rectangular resample.
- Plot recipes (RecipesBase, no Plots dep): `tsdiagram(sol)`,
  `mapplot(map_or_turbomachine)` with operating point.
- Exported unit helpers (`R_to_K`, `psia_to_Pa`, `lbps_to_kgps`, BTU,
  density, rpm, hp, all with inverses) — examples no longer redefine them.
- CI workflow (GitHub Actions): `Pkg.test()` on latest stable; pinned
  Manifest dropped in CI so runs resolve fresh against [compat].

### 11. Optimization showcase
The never-built verification item 6 from the original plan: maximize cycle
η over (PR, ε, mixture MW) with ForwardDiff + Optim.jl. Cheap once the
ideal-gas path is the backend, and it is the demonstration of why this
tool beats NPSS (exact gradients).

### 12. OTAC-style meanline turbomachinery
Geometry-based component models instead of scaled generic maps, à la
OTAC (Jones, NASA Glenn). Scope notes:
- Space Brayton machines are radial: use NASA Glenn radial correlations
  (Galvas for centrifugal compressors; Glassman/Wasserbauer/Rohlik for
  radial-inflow turbines) rather than general axial meanline.
- Pragmatic first step: a meanline *design* tool that generates a map
  from geometry once, consumed by the existing map machinery — most of
  the value without velocity triangles inside every Newton iteration.
- Full OTAC-equivalent (meanline inside the solve) and AD-through-geometry
  sensitivities come after.

### Backlog (lower priority)
- **FPT AD via implicit-function rule** (deferred with trigger): dT = dh/cp
  on the bisection inversions is ~half a day, but gradients through
  `Gridded(Linear())` tables are piecewise-constant — useful AD also needs
  cubic interpolation. Do this only if a table-only fluid enters the model
  (e.g., Oil.fpt / H2O.fpt for a heat-rejection coolant loop, which have no
  closed form to escape to).
- Recuperator uses `cp·ΔT` instead of enthalpy differences (fine for
  monatomic He-Xe, approximate in general).
- Network validation: error if `one_pass!` leaves elements unprocessed;
  convergence checks watch only outlet temperatures (not Pt, W, power).
- `Splitter`/`Mixer` branch networks: add end-to-end test for the BRU bleed
  configuration.
