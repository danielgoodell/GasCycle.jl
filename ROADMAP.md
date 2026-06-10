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

### 4. FPT AD via implicit-function rule  ← NEXT
`T_from_h` / `T_from_s` in FPTFluid use bisection that drops Dual derivatives,
so AD currently only fully works with the ideal-gas backend. Add a custom
ForwardDiff rule (implicit function theorem: dT = dh / cp etc.) so
gradient-based optimization works with real He-Xe table data — a headline
goal over NPSS.

### Backlog (lower priority)
- Recuperator uses `cp·ΔT` instead of enthalpy differences (fine for
  monatomic He-Xe, approximate in general).
- Network validation: error if `one_pass!` leaves elements unprocessed;
  convergence checks watch only outlet temperatures (not Pt, W, power).
- `Splitter`/`Mixer` branch networks: add end-to-end test for the BRU bleed
  configuration.
