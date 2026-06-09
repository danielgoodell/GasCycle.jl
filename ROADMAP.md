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

### 1. Fix the off-design map residual formulation  ← IN PROGRESS
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

### 2. Two quick wins
- **`Pkg.test()` fails**: `test/runtests.jl` does `using Pkg; Pkg.activate(...)`,
  which breaks the isolated test env. Remove it and declare test deps via
  `[extras]`/`[targets]` in Project.toml.
- **Stale turbine PR display**: `examples/bru_tit_sweep.jl:79` still prints
  `turb.PR`, which is left unchanged in `:pressure_closure` mode. Add a
  `pressure_ratio(::Turbine)` accessor and use it in both examples
  (recuperated_brayton.jl already has a local workaround).

### 3. Off-design TIT sweep validation (plan Phase 7)
Once item 1 is in, build the turbine-inlet-temperature sweep (100% → 60% of
design) with scaled maps on the BRU cycle: confirm shaft power balance holds
and the maps track. This is the test that proves off-design capability is real.

Known issue to watch: in off-design mode the solver handles closed-loop
back-edges by legacy fixed-point seeding inside `one_pass!` (stale state
between Newton iterations). For closed loops this may need the back-edge
states folded into the outer Newton vector.

### 4. FPT AD via implicit-function rule
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
