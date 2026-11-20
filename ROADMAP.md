# GasCycle.jl Roadmap

Status as of 2026-11-19. Completed work (off-design solver, NPSS
cross-validation campaign, NobleGasMixture backend, cold-end elements,
off-design HX effectiveness from flow-scaled UA, PolynomialLiquid coolants,
NPSS-syntax R-line map reader with digitized BRU argon maps, gradient/Hessian
optimization example) is documented in the README, `validation/RESULTS.md`,
and git history — this file tracks only what's still open.

## Where the project stands

Design-point and map-based off-design analysis of the gas loop are complete
and validated against NPSS 3.3 to print precision. The loop closes physically
through a cooler or radiator, with gas and liquid property backends that are
all AD-exact. `examples/optimize_cycle.jl` demonstrates exact gradients and
Hessians through the full solve.

What is missing is everything *around* the thermodynamic cycle: the shaft
ends in a scalar `P_load` (no electrical system, so no DC output figure), and
losses are scattered across elements with no system-level accounting. User
targets and free variables (item 1) are in place, so outputs such as a DC
power figure become solvable requirements as soon as they exist. Real
space power work is specified at the DC bus ("deliver N kWe"), sized by
specific mass, and controlled by inventory and parasitic load.

Priority driver: mission needs (SR-1 Freedom) — DC power delivery, inventory
control, and slow transients.

## Design principles for new work

These hold for every item below; they are what make the tool better than NPSS
rather than a copy of it.

- **AD through everything.** New elements are parametric in `T<:Real`;
  fields the solver writes are typed `Real` so they accept an outer
  gradient's Duals (see CLAUDE.md).
- **Every element reports its losses** through one interface (item 2), so a
  system loss ledger and energy-closure check come for free.
- **Multi-domain ports.** Fluid ports exist; electrical and heat ports are
  added alongside them (items 3, 4), sharing the same network graph and
  Newton vector.
- **Composable assemblies.** Higher-level components (a turbo-alternator-
  compressor, a power conversion unit) are built from elements, not as
  monoliths (item 5).
- **Discretized components solve their own interiors.** A 1D/2D heat
  exchanger or a heat-pipe radiator iterates its internal nodes inside its
  own compute step (marching or a local Newton, with implicit-function-theorem
  derivatives) and exposes only boundary states to the cycle solve. The
  global Newton vector stays at a few dozen unknowns, so the small built-in
  dense solvers (`src/solver/NonlinearSolvers.jl`) remain the right tool.
  Tightly coupled cases that genuinely need one large sparse system get a
  sparse or NonlinearSolve.jl backend as a package extension behind the
  `NonlinearBackend` interface — never a hard dependency.
- **Lightweight core dependencies.** Core deps are ForwardDiff,
  Interpolations, RecipesBase, and stdlibs. Heavy stacks (SciML,
  DifferentialEquations, Plots, FMI) enter only as package extensions.

---

## Phase A — Steady-state system modeling

### 1. User-defined targets and free variables (solver) — DONE 2026-11-19
`add_free!` (element/shaft field, seed-state field, or get/set closures, with
`ref` scaling and box bounds) and `add_target!` (any `net -> Real`) join the
element unknowns and back-edges in one Newton vector, in design and
off-design. Off-design solves became differentiable at the same time
(solver-written fields typed `Real`; Jacobians on primal values with
derivative parts converged after the primal). Tested against finite
differences, including second derivatives and d/d(target value).
Follow-ups:
- A free variable on a concretely typed design field (e.g.
  `HeatSource{Float64}.TtExit`) cannot hold an outer gradient's Duals; the
  user must construct that element in the gradient's number type (clear
  error message). If this proves annoying, type more parameter fields `Real`
  rather than retyping networks. Cost reference: making the solver-written
  fields `Real` plus type-generic residual assembly slowed the benchmark
  solves by 8–36% (off-design 16.5 → 22.4 µs).
- The coupled solve uses finite-difference Jacobians (the unknowns live in
  typed fields, so ForwardDiff can't seed them). Fine at current sizes; an
  exact alternative is a ForwardDiff Jacobian over a retyped copy of the
  network, if conditioning ever demands it.
- Report free variables/targets in `summary(sol)` (currently only
  `solve!(net; verbose=true)` prints them).

### 2. Loss and energy accounting
One interface, `losses(el) -> Vector{Loss}`, where each entry has a name, a
category (`:pressure`, `:heat_leak`, `:mechanical`, `:electrical`,
`:heat_rejection`, `:aero` …), a power-equivalent value [W], and where the
heat goes (ambient, a fluid stream, a radiator).
- System ledger: `loss_ledger(net)` tables every loss; `energy_closure(net)`
  checks Q_in = P_out + Σ losses + Q_rejected to solver tolerance (becomes a
  standard test on every example cycle).
- Second-law view: entropy generation / exergy destruction per element
  (T₀·S_gen, with a user-chosen dead state — radiator sink temperature is the
  natural T₀ for space). Cheap because `entropy` already exists everywhere.
- Waterfall output in `Summary.jl` and a plot recipe: reactor heat → shaft
  power → alternator → cable → PCU → DC bus, with each loss called out.
- Build this early: every element added after it reports losses from day one.

### 3. Electrical domain: alternator → cable → power conversion → DC bus
The shaft currently ends in a scalar `P_load`. Replace that with an
electrical chain whose output is the DC power the mission actually specifies.
- **Electrical port** carrying power, voltage, and current (DC) or power,
  RMS voltage, frequency, and power factor (AC). Start with power-flow-level
  modeling (no waveforms); NPSS 3.3's electric ports are the same idea.
- **Alternator** (on a `Shaft`, replacing `P_load`): electrical frequency
  from poles and speed; efficiency from a loss model (copper I²R, core loss
  ∝ f^α, windage, bearing) or an efficiency map vs. speed and load fraction.
  Its losses are *heat*, delivered through a heat port (item 4) to whatever
  cools it: the working fluid (gas-cooled alternator in the cycle, e.g. a
  compressor-discharge bleed), the main coolant loop, or a secondary loop.
- **Cable** (the electrical "duct"): resistance from length, gauge, material,
  and conductor temperature (R(T) feedback); I²R loss and voltage drop; AC
  or DC.
- **Power conversion unit** (rectifier / DC-DC / inverter): AC→DC with a loss
  model P_loss = a + b·P + c·P² (standby + switching + conduction) or an
  efficiency curve; output voltage regulation; its losses also go to a heat
  port (electronics cold plate on a coolant loop).
- **DC bus and loads**: user load, plus a **parasitic load radiator (PLR)**,
  the standard space-Brayton speed control — the PCU dumps surplus power
  into the PLR to hold shaft speed. Bus power is the canonical target for
  item 1.
- **Shaft mechanical losses** (bearing, windage, seal drag) as functions of
  speed and gas density — significant at the kW scale.
- Validate against BRU/SNAP-era alternator and PCU data (NASA reports
  already in `reference/` and the BRU literature).

### 4. Heat ports and thermal coupling between subsystems
Losses in the alternator, PCU, bearings, and cables end up as heat in a
fluid stream, and that coupling is what closes the system-level energy
balance.
- **Heat port**: a heat flow [W] with a source temperature, connecting a
  loss-producing element to a fluid stream element (generalized
  `HeatSource` with Q supplied by another element rather than fixed).
- Heat-flow links can create new algebraic loops (alternator losses depend
  on coolant temperature through R(T), coolant temperature depends on losses)
  — they become back-edges / Newton unknowns like the existing recuperator
  loop.
- **Pump element** for liquid coolant loops: head-flow curve or fixed
  efficiency, with pump electrical power drawn from the DC bus as a
  parasitic. Net power = bus power − parasitics.
- **Valves / bypass / pressure-loss elements** (Cv or K-factor): needed for
  coolant loop balancing and for turbine-bypass power control (item 6).

### 5. Assemblies (high-level components)
NPSS-style hierarchical assemblies: a named group of elements with exposed
ports, defaults, and its own targets.
- `Assembly` type that the network flattens at build time (no runtime cost,
  same Newton vector).
- Shipped assemblies as building blocks and examples:
  - **TAC** (turbine–alternator–compressor on one shaft, with bearing and
    windage losses and alternator cooling path),
  - **PCU** (TAC + recuperator + gas cooler, i.e. the BRU power conversion unit),
  - **Heat rejection subsystem** (pumped coolant loop + radiator),
  - **Power management and distribution** (cable + PCU + bus + PLR).
- Assemblies make multi-unit systems (e.g. 2 PCUs sharing a reactor and
  radiator) a one-liner, and are the natural unit for mass estimation (item 7).

### 6. Inventory control, volumes, and control elements
Real CBC power control is charge-pressure (inventory) control, plus PLR and
turbine bypass.
- Make loop inventory (total gas mass/moles) a model quantity: assign a
  volume to every component; where components are currently joined
  directly, connecting ducts may need to be added as the volume holders.
- Steady-state inventory mode: loop mass becomes a free variable with a
  pressure-closure residual (replaces the fixed-seed-Pt assumption;
  resolves the dangling-pressure caveat from the off-design TIT sweep).
  Uses the item 1 machinery rather than a special case.
- Mass accounting in components with large T/P gradients is the hard part:
  m = ∫ρ dV over a strong axial temperature gradient (recuperator!) is not
  well-approximated by inlet/outlet averages. Use segmented volumes or an
  analytic mean density per component (e.g., ideal-gas with linear-T duct
  has closed-form m = PV/(R·T_lm) using log-mean temperature); validate
  segment-count convergence.
- **Schedules and limiters**: setpoint schedules (e.g. charge pressure vs.
  demanded power) and min/max limiters on free variables (speed, TIT,
  surge margin) with active-limit switching. Limiter switching is a
  non-smooth branch — implement as smooth min/max or as an outer active-set
  loop, never as a branch on Dual values.
- Sets up transients: volume + mass per component is exactly the
  capacitance structure item 9 needs.

### 7. System mass and sizing metrics
For space power the real objective is usually specific mass (kg/kWe), not
efficiency.
- Per-element mass models: radiator (area × areal density, from the sink
  temperature and rejected heat), heat exchangers (from UA and effectiveness),
  turbomachinery and alternator (scaling laws from power and speed), cable
  (length × gauge), PCU (W/kg), coolant inventory.
- `system_mass(net)`, `specific_mass(net)`, radiator area, and
  footprint/volume as targetable and optimizable quantities.
- Extend `examples/optimize_cycle.jl` to the mission trade:
  minimize kg/kWe over (PR, recuperator ε, mixture MW, TIT, radiator
  temperature) subject to a DC power target (item 1). This also closes out
  the original "optimize over PR, ε, mixture MW" goal.

### 8. Multi-point design
Size the cycle at one point while satisfying constraints at others
(beginning/end of mission, hot/cold sink, part power, a degraded case).
- One shared design vector, several operating-point networks solved
  together; design-point sizing plus off-design targets in one Newton or
  optimizer problem (pyCycle's multi-point design).
- AD through both design and off-design makes the coupled Jacobian exact —
  a capability NPSS approximates with nested finite-difference loops.

### 9. Solver robustness and case management
- Continuation for sweeps: warm-start each point from the last converged
  state, with step halving on failure; arc-length continuation through turning
  points (e.g. near surge).
- Diagnostics on non-convergence: report which residuals are worst and which
  free variables hit bounds, in engineering units.
- Network validation: error if `one_pass!` leaves elements unprocessed;
  convergence checks currently watch only outlet temperatures (not Pt, W,
  power).
- Case runner: run a table of cases (CSV in), collect named outputs
  (CSV / DataFrame out). NPSS 3.3 added native CSV I/O for the same reason.

---

## Phase B — Dynamics and controls

### 10. Transients — shaft and electrical first, then thermal
The original mission of this architecture (reactor coupling). Stage it:
- Shaft speed dynamics: I·ω·dω/dt = turbine power − compressor power −
  alternator load − mechanical losses. With the electrical chain (item 3)
  and PLR, startup, load-step, and loss-of-load events become solvable
  with just shaft inertia.
- Slow thermal transients: thermal capacitance (metal + gas) per component
  using the volumes from item 6, wrapped with DifferentialEquations.jl
  (as a package extension, so steady-state users don't load SciML).
  Reactor side couples here (point kinetics or supplied Q(t)).
- Mixed time scales: fast electrical/shaft dynamics alongside slow thermal
  dynamics. Use stiff implicit integrators on the whole DAE first; consider
  nested/multirate stepping (as in NPSS 3.3's nested transient solver) only
  if profiling demands it.
- Controller elements (PI on PLR for speed, inventory valve for power) using
  the schedules and limiters from item 6.

### 11. Linearization for control design
- `linearize(net; states, inputs, outputs)` → state-space (A, B, C, D) at
  any converged operating point, from exact AD Jacobians (NPSS uses
  finite-difference perturbation).
- Steady-state version (output/input gain matrices) is available as soon as
  item 1 lands; full dynamic A/B/C/D after item 10.
- Export to ControlSystems.jl; a linear-model-vs-nonlinear-transient step
  response comparison as the validation test.

### 12. Co-simulation interfaces (lower priority)
- FMU export/import (FMI.jl) so GasCycle models can couple to Simulink or
  Modelica controller models — NPSS 3.3 added FMI for exactly this.
- External function hooks (Julia `ccall` covers C/Fortran natively), e.g. a
  reactor kinetics code supplying Q(t).

---

## Phase C — Component fidelity

### 13. Geometry-based heat exchangers
- UA and pressure drop from geometry plus correlations (Kays & London
  surfaces; plate-fin and tube-in-shell for the recuperator and gas cooler),
  replacing specified ε/ΔP.
- Metal mass from geometry → feeds item 7 (mass) and item 10 (thermal
  capacitance).
- Radiator: fin efficiency, view factor, segmented radiator with axial
  coolant temperature drop.

### 14. OTAC-style meanline turbomachinery
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
- Turbomachinery mass scaling from geometry feeds item 7.

### 15. Additional working fluids: air, N₂, CO₂
Real systems are commonly first tested with cheap fluids: air (low-temp
checkout, motoring), N₂, or CO₂, before committing the expensive xenon
inventory. Ar/Kr (and He-Ar, He-Kr) are already covered by
`NobleGasMixture`; air/N₂/CO₂ are di-/triatomic with temperature-dependent
Cp, so monatomic kinetic theory does not apply. Implement a
`ThermallyPerfectGas` backend:
- Cp(T) from NASA 7- or 9-coefficient polynomials (h, s by closed-form
  integration — keeps AD support exact like IdealGasFluid).
- Ideal-gas EOS is adequate for checkout/motoring conditions; revisit only
  if someone needs sCO₂-style operation near the CO₂ critical point.
- Transport via Sutherland law or polynomial fits per gas.
- Air as a fixed-composition pseudo-species (standard dry air), N₂ and CO₂
  as pure species.
- Prerequisite: fix the recuperator's `cp·ΔT` (see backlog) — it is only
  exact for constant-cp gases.

---

## Suggested sequencing

1 (targets) → 2 (loss accounting) → 3 + 4 (electrical chain with heat
coupling) → 5 (assemblies) → 6 (inventory/controls) → 7 (mass) → 8, 9 →
Phase B. Items 1 and 2 are small, cross-cutting, and change how every later
item is written, so they go first. Item 9's network-validation and
continuation pieces can be picked up opportunistically at any point.

## Backlog (lower priority)

- **Native fluid-definition format + conversion utilities** (deferred with
  trigger, 2026-06-12): one declarative TOML per fluid (kind, SI
  coefficients, validity range, machine-readable provenance/fit residuals)
  loaded by a single `load_fluid(path)`; a generic sample-and-fit importer
  (evaluate any oracle — CoolProp.jl, REFPROP dump, FPT table, collaborator
  script — on a grid, fit to the native representation, stamp residuals);
  FPT demoted to a pure export target regenerated for NPSS
  (`data/generate_liquids.jl` becomes the export half). CoolProp stays a
  generation-time oracle only, never a runtime dependency, so AD-exactness
  is preserved. Build when any trigger hits: coolant library grows past
  ~5 fluids, a third consumer tool appears, or a fluid needs re-fitting and
  the anchors aren't machine-readable.
- **NPSS map reader cross-check**: `read_npss_map` now follows NPSS R-line
  map syntax and reads the digitized BRU argon maps (NASA TM X-2129), but
  those files were generated here, not by NPSS — cross-check against a map
  file NPSS itself has consumed when one is available. Maps operated near
  choke (vertical speed-line segments) eventually want an Rline-native
  solver coordinate instead of the rectangular resample. Compressor map
  efficiency is a first-pass digitization (see file header).
- **Sink-HX effect definition** (open item from the NPSS campaign): NPSS's
  HeatSinkHX behaves as ε = 0.9405 vs its 0.946 setting; gas loop
  unaffected. See `validation/RESULTS.md`.
- **FPT AD via implicit-function rule** (deferred with trigger): dT = dh/cp
  on the bisection inversions is ~half a day, but gradients through
  `Gridded(Linear())` tables are piecewise-constant — useful AD also needs
  cubic interpolation. Do this only if a table-only fluid enters the model
  (e.g., Oil.fpt / H2O.fpt for a heat-rejection coolant loop, which have no
  closed form to escape to).
- Recuperator uses `cp·ΔT` instead of enthalpy differences (fine for
  monatomic He-Xe, approximate in general; must be fixed before item 15).
- `Splitter`/`Mixer` branch networks: add end-to-end test for the BRU bleed
  configuration (also the natural flow path for a gas-cooled alternator,
  item 3).
