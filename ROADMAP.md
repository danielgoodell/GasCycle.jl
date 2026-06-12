# GasCycle.jl Roadmap

Status as of 2026-06-11. Completed work (off-design solver, NPSS
cross-validation campaign, NobleGasMixture backend, cold-end elements,
off-design HX effectiveness from flow-scaled UA, output/plotting/CI batch)
is documented in the README, `validation/RESULTS.md`, and git history —
this file tracks only what's still open.

## Where the project stands

Design-point and map-based off-design analysis are complete and validated
against NPSS 3.3 to print precision. The analytic `NobleGasMixture` backend
(thermo + transport) replaces FPT tables for He-Xe work. The loop closes
physically through a cooler or radiator.

Priority driver: mission needs (SR-1 Freedom) put **inventory control and
slow transients** at the top of the physics track.

## Ordered next steps

### 1. Inventory control and per-component volume bookkeeping
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

### 2. Transients — shaft dynamics first, then thermal
The original mission of this architecture (reactor coupling). Stage it:
- Shaft speed dynamics: I·ω·dω/dt = unbalanced shaft power. The power-
  balance machinery already computes the imbalance; startup, load-step,
  and loss-of-load events become solvable with just shaft inertia.
- Slow thermal transients: thermal capacitance (metal + gas) per component
  using the volumes from item 1, wrapped with DifferentialEquations.jl.
  Reactor side couples here (point kinetics or supplied Q(t)).

### 3. Additional working fluids: air, N₂, CO₂
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

### 4. Optimization showcase
Maximize cycle η over (PR, ε, mixture MW) with ForwardDiff + Optim.jl.
Cheap now that the property backends are AD-complete, and it is the
demonstration of why this tool beats NPSS (exact gradients).

### 5. OTAC-style meanline turbomachinery
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

## Backlog (lower priority)

- **NPSS map reader cross-check**: `read_npss_map` was validated against a
  synthetic format-faithful fixture only — cross-check against a real NPSS
  map when one is available. Maps operated near choke (vertical speed-line
  segments) eventually want an Rline-native solver coordinate instead of
  the rectangular resample.
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
  monatomic He-Xe, approximate in general).
- Network validation: error if `one_pass!` leaves elements unprocessed;
  convergence checks watch only outlet temperatures (not Pt, W, power).
- `Splitter`/`Mixer` branch networks: add end-to-end test for the BRU bleed
  configuration.
