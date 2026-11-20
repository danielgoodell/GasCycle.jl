# GasCycle.jl — API Reference

Everything exported by `using GasCycle`, grouped by layer. All quantities are
SI (K, Pa, kg/s, W, J/kg) unless noted. See [GUIDE.md](GUIDE.md) for the
workflow these pieces fit into.

## Fluid property backends

All backends implement the same `FluidProperties` interface; elements never
care which one they're given.

| Constructor | Backend | Notes |
|---|---|---|
| `HeXe(M_molar)` | analytic noble-gas mixture | He-Xe by molecular weight in kg/kmol, e.g. `HeXe(83.8)` |
| `NobleGasFluid(gas)` | analytic, pure gas | `gas` ∈ `HELIUM`, `NEON`, `ARGON`, `KRYPTON`, `XENON` |
| `NobleGasMixture(gas1, gas2, x1; name)` | analytic, any binary | `x1` = mole fraction of `gas1` |
| `FPTFluid(path; bounds=:error, s_interp=:log_pressure)` | NPSS FPT table file | `s_interp=:linear` reproduces NPSS's linear-in-P entropy interpolation (for exact NPSS replication only) |
| `IdealGasFluid(; M_molar)` | constant-cp monatomic ideal gas | exact closed-form inversions |
| `HeXeIdealGas(x_He)` | same, He-Xe by He mole fraction | |
| `PolynomialLiquid(; cp, rho, mu=[], k=[], name, T_min, T_max)` | coolant with polynomial T-dependence | coefficients ascending in T [K], scalars = constant; `PolynomialLiquid(path)` reads `data/Water.fpt`, `data/DC200.fpt`, `data/WaterEG50.fpt` |
| `ConstantPropertyLiquid(; cp, rho, name="liquid")` | constant-property coolant | also `ConstantPropertyLiquid(path)` for the legacy NPSS-replication files (`data/H2O.fpt`, `data/Oil.fpt`) |

The analytic `NobleGasMixture` family (virial EOS + corresponding-states
transport, Tournier/El-Genk AIAA 2006-4154) is the default choice for He-Xe:
real-gas thermo, transport properties, full AD including d/d(mixture ratio),
no table file. `FPTFluid` is kept for apples-to-apples NPSS cross-validation.

### Property functions

Forward, on any backend (`fl`, `T` [K], `P` [Pa]):

| Function | Returns |
|---|---|
| `GasCycle.cp(fl, T, P)` | specific heat [J/(kg·K)] — **not exported**, import explicitly |
| `enthalpy(fl, T, P)` | h [J/kg], referenced to 0 K |
| `entropy(fl, T, P)` | s [J/(kg·K)] |
| `density(fl, T, P)` | ρ [kg/m³] |
| `gamma(fl, T, P)` | cp/cv [-] |
| `viscosity(fl, T, P)` | μ [Pa·s] — `NobleGasMixture`; `PolynomialLiquid` when the file has μ data |
| `conductivity(fl, T, P)` | k [W/(m·K)] — `NobleGasMixture`; `PolynomialLiquid` when the file has k data |
| `prandtl(fl, T, P)` | Pr [-] — same availability as μ and k |

Inversions (Newton/closed-form per backend; pass `T_guess` when you have one —
element code always does):

```julia
T_from_h(fl, h, P; T_guess=500.0)
T_from_s(fl, s, P; T_guess=500.0)
h_from_s(fl, s, P)        # isentropic: h at entropy s and pressure P
```

`enthalpy(state)` / `entropy(state)` etc. also accept a `FluidState` directly.

### States and ports

```julia
s = FluidState(Pt, Tt, W, fluid)     # total pressure [Pa], total temp [K], flow [kg/s]
update(s; Tt=..., Pt=..., W=...)     # copy with fields replaced
p = Port(s); p[]                     # Port wraps a state; p[] dereferences
```

## Elements

Constructors with all keywords and defaults. Every element is mutable — you
may set fields (`heat.TtExit = ...`, `rad.mode = :fixed_area`) between solves.

### Compressor
```julia
Compressor(name; PR=2.0, η_poly=0.87, η_type=:polytropic,
           map=nothing, reynolds=nothing, mode=:design)
```
- `mode=:design` — `PR` fixed by the user.
- `mode=:off_design` — PR and η from `map` evaluated forward at `(Nc, Rline)`;
  the R-line coordinate `Rline` is the solver unknown, closed by a
  flow-continuity residual (map corrected flow = actual). `map` is a
  `TurbomachineMap` (e.g. `CompressorMap` or `FunctionMap`).
- `reynolds` — optional `ReynoldsModel` feeding the map's Reynolds tables
  (default `nothing` ⇒ no correction).
- `η_type` — `:polytropic` (small-stage) or `:isentropic` (NPSS `effDes`).

### Turbine
```julia
Turbine(name; PR=2.0, η_poly=0.90, η_type=:polytropic, map=nothing,
        reynolds=nothing, mode=:design, P_exit=101325.0)
```
- `mode=:design` — `PR` fixed.
- `mode=:pressure_closure` — PR computed so the exhaust hits `P_exit`
  (standard for closed-cycle design solves).
- `mode=:off_design` — PR (the map's native input) is the solver unknown;
  the map yields flow `Wp` and η at `(Np, PR)`, closed by a flow-continuity
  residual. Loop pressure closure is the network back-edge. Seed `PR` at the
  design expansion ratio. `reynolds` as Compressor.

### Duct
```julia
Duct(name; dPqP=0.02)        # isenthalpic, Pt_out = Pt_in·(1−dPqP)
```

### HeatSource
```julia
HeatSource(name; Q=0.0, TtExit=1200.0, dPqP=0.02, mode=:fixed_TtExit)
```
- `:fixed_TtExit` — exit temperature imposed, `Q` computed (reactor at set TIT).
- `:fixed_Q` — heat input imposed, exit temperature computed.

### HeatExchanger (recuperator / cooler)
```julia
HeatExchanger(name; ε=0.90, UA=nothing, UA_exp=0.8,
              W_hot_des=nothing, W_cold_des=nothing,
              dPqP_hot=0.01, dPqP_cold=0.01)
```
Counter-flow, ε-NTU, `Q = ε·Q_max` with Q_max from endpoint enthalpies. Mode
is inferred from the keywords:
- `:effectiveness` (default) — fixed ε.
- `:UA` (give `UA` [W/K]) — ε recomputed each pass from NTU = UA/C_min.
- `:scaled_UA` (give `UA` + `W_hot_des` + `W_cold_des`, or call `size_UA!`) —
  additionally scales UA with stream flows, Colburn-style per side:
  `UA = 2·UA_des / [(W_hot_des/W_hot)^n + (W_cold_des/W_cold)^n]`, `n = UA_exp`.

```julia
size_UA!(hx; UA_exp=hx.UA_exp) -> UA
```
Call after a converged design solve: inverts ε-NTU at the design state
(from `:effectiveness`) or keeps the given UA (from `:UA`), records design
flows, switches to `:scaled_UA`. Off-design at design conditions then
reproduces the design ε exactly.

Ports: `:hot_inlet`, `:cold_inlet` (alias `:inlet`), `:hot_outlet`,
`:cold_outlet` (alias `:outlet`). The serial `connect!` chain runs through the
cold side; the hot side is wired by `add_hx_pair!` (recuperator) or
`connect_port!` (cooler), and a coolant from outside the loop by `set_boundary!`.

### Radiator
```julia
Radiator(name; A=0.0, emissivity=0.85, T_sink=200.0, TtExit=400.0,
         dPqP=0.01, mode=:fixed_area, N_seg=50)
```
Segmented `σ·ε_r·A·(T⁴−T_sink⁴)` rejection (Heun marching, `N_seg` segments).
- `:fixed_TtExit` — design sizing: required `A` computed and stored.
- `:fixed_area` — off-design: exit temperature responds to the operating point
  (requires `A > 0`). Typical flow: size with `:fixed_TtExit`, then flip
  `rad.mode = :fixed_area`.

### Splitter / Mixer (bleed branches)
```julia
Splitter(name; fracs=[0.98, 0.02])   # fracs must sum to 1; W split, Tt/Pt unchanged
Mixer(name; n_inlets=2)              # mass + enthalpy balance; Pt_out = min(Pt_i)
```
Splitter outlets: `:outlet` (fracs[1]), `:bleed_outlet`, `:outlet_3`, …
Mixer inlets: `:inlet`, `:bleed_inlet`, `:inlet_3`, … All mixer inlets must
share the same fluid backend.

### Shaft
```julia
Shaft(name; N=10000.0, mode=:design, P_load=0.0)
```
- `:design` — no residual; `power_balance(shaft)` reports net power.
- `:off_design` — `N` [rpm] is a solver unknown; residual is the power balance
  `ΣW_turb − ΣW_comp − P_load = 0` (`P_load` = generator extraction [W]).
Added to the network with `add_shaft!`, which broadcasts `N` to the linked
turbomachines for their corrected-speed lookups.

### Element queries

| Function | Meaning |
|---|---|
| `specific_work(comp_or_turb)` | Δh across the machine [J/kg] |
| `pressure_ratio(comp_or_turb)` | actual Pt ratio from solved states (correct in `:pressure_closure`) |
| `Q_transferred(hx)` | hot→cold heat [W] |
| `Q_rejected(rad)` | radiated heat [W] |
| `power_balance(shaft)` | turbine − compressor power [W] |

## FlowNetwork

```julia
net = FlowNetwork()
add!(net, el1, el2, ...)                       # register elements (Shafts auto-sorted)
connect!(net, a => b => c => a)                # serial chain via :outlet → :inlet;
                                               # repeating the first element closes the loop
connect_port!(net, src, :port, dst, :port; back_edge=false)   # explicit edge (branches)
add_shaft!(net, shaft; drives=comp, driven_by=turb)           # accept elements or vectors
add_hx_pair!(net, hx; hot=turb)                # back-edge: hot.outlet → hx.hot_inlet
set_state!(net, first_el; Pt, Tt, W, fluid)    # entry state (or initial guess + W if loop closed)
set_boundary!(net, el, :port; Pt, Tt, W, fluid)  # fixed external stream (e.g. coolant)
one_pass!(net)                                 # single propagation (rarely needed directly)
```

`back_edge=true` on `connect_port!` makes the destination state a solver
unknown seeded each iteration — used to close the loop through a cooler or
radiator so the compressor inlet floats.

## Solver

```julia
sol = solve!(net; tol=1e-6, maxiter=100, verbose=false, alg=nothing)   # -> SolveResult
```

Mode is automatic: with no element unknowns and no free variables, a back-edge
Newton (`NewtonBacktracking`, ForwardDiff Jacobian) on the back-edge
`[Tt, Pt]` states; otherwise (off-design map line coordinates — compressor
`Rline`, turbine `PR` — shaft `N`, and user free variables) everything joins one
normalized vector solved by `DoglegTrustRegion` (finite-difference Jacobian).
`tol` is on the 2-norm of the normalized residuals. Non-convergence warns and
returns `status = :failed`. The system must be square (element unknowns + free
variables = element residuals + targets), otherwise `solve!` throws.

Both solves are differentiable: an outer `ForwardDiff` gradient/Hessian flows
through the converged solution (implicit function theorem), in design and
off-design, with or without targets.

Solver backends (subtypes of `NonlinearBackend`; pass one as `alg` to override
the mode default):

| Backend | Options (defaults) |
|---|---|
| `NewtonBacktracking` | `jacobian=:forwarddiff` (or `:finitediff`), `max_backtracks=10` |
| `DoglegTrustRegion` | `jacobian=:finitediff` (or `:forwarddiff`), `initial_radius=1.0`, `max_radius=1e3` (both relative to `max(‖u0‖,1)`) |

`jacobian=:forwarddiff` only works where the unknowns are not stored in
concretely typed element fields, i.e. back-edge-only design solves (the
default there).

`SolveResult` fields: `status` (`:success`/`:failed`), `iterations`,
`residual_norm`, `net`. `sol["Comp"]` looks an element up by name.

```julia
net_power(sol)          # Σ turbine − Σ compressor power [W]   (also net_power(net))
cycle_efficiency(sol)   # net power / Σ heat input             (also cycle_efficiency(net))
```

### Targets and free variables

NPSS-style solver pairs: free variables are model inputs the solver varies;
targets are model outputs it drives to a value. They join the same Newton
vector as the element unknowns and back-edges, in design or off-design mode.

```julia
add_free!(net, el, :field; name, ref, lower, upper)   # element (or Shaft) field
add_free!(net, :W; ...)                               # seed-state field: :W, :Tt, :Pt
add_free!(net, "name"; get=() -> x, set! = v -> ..., ...)   # anything, via closures
add_target!(net, f, value; name, scale)               # drive f(net) → value
clear_targets!(net)                                   # remove all pairs
```

| Keyword | Meaning (default) |
|---|---|
| `ref` | normalization of the free variable (`\|initial value\|`, or 1 if zero) |
| `lower`, `upper` | box bounds; trial values are clamped (`-Inf`, `Inf`) |
| `scale` | target residual is `(f(net) − value)/scale` (`\|value\|`, or 1 if zero) |

The current field value is the initial guess. `f` receives the network after
each pass, so `net_power`, `cycle_efficiency`, or any closure such as
`n -> comp.outlet[].Tt` works. `value` may carry Duals (differentiate w.r.t. a
requirement). Under an outer gradient, a free variable on a concretely typed
field (e.g. `HeatSource{Float64}.TtExit`) must be constructed in the gradient's
number type — `solve!` explains how if it is not. Solver-owned fields
(compressor `PR`/`η_poly`/`Rline`, turbine `PR`/`η_poly`, shaft `N`/`P_load`)
are `Real`-typed and accept anything.

## Performance maps

Maps are evaluated **forward in their native coordinates** (no inversion) through
a single interface, so the element never touches table data — a map can be a
tabulation from any file format or a plain function ("map as a script").

```julia
eval_map(m, Nc, Rline[, rc])    # CompressorMap -> (; Wc, PR, eff)
eval_map(m, Np, PR[, rc])       # TurbineMap    -> (; Wp, eff)
corrected_speed(N, Tt)          # N / √(Tt/288.15)
corrected_flow(W, Tt, Pt)       # W·√(Tt/288.15) / (Pt/101325)
scale_map(base::CompressorMap; Nc_des, Wc_des, PR_des, eta_des)
scale_map(base::TurbineMap;    Np_des, Wp_des, PR_des, eta_des)
```

`scale_map` solves the four scale factors so the (unscaled) base map passes
exactly through the design point at the map's own design anchors. `rc` is the
value looked up in the map's Reynolds-correction tables (`rc === nothing` ⇒ no
correction) — supplied on the element by a `ReynoldsModel`
(`ReDesIndex(Re_des)`, `RawRe()`, or `FunctionReynolds(f)`).

`FunctionMap(f; line_des=2.0)` wraps a closure `f(speed, line, rc) -> NamedTuple`.

NPSS `.map` (NEO `Subelement`/`Table`) files — the loader is fully positional
(NPSS cares about hierarchy, not names):

```julia
read_npss_map("file.map")                       # -> ParsedSubelement tree
cm = compressor_map("compressor.map";           # -> unscaled CompressorMap
                    flow="TB_Wc", pr="TB_PR", eff="TB_eff")
tm = turbine_map("turbine.map"; flow="TB_Wp", eff="TB_eff")   # -> TurbineMap
```

The builders read the file's design anchors and any nested Reynolds (`S_Re`)
subelement; run the result through `scale_map` at your design point.

## Output

```julia
summary(sol)                    # print NPSS-style station/component/cycle tables
sprint(summary, sol)            # same, captured as a String
stations(sol; branches=true)    # Vector of "Element.port" => FluidState, flow order
tsdiagram(sol)                  # Plots.jl recipe: T-s diagram with stations
mapplot(m); mapplot(comp)       # map speed lines; with operating point if element
losswaterfall(L)                # Plots.jl recipe: exergy waterfall of a LossLedger (or sol)
```

`stations` walks the main loop in physical order (through HX hot sides via
back-edges), then bleed branches.

### Energy and loss accounting

```julia
eb = energy_balance(sol)        # first law, all [W]
L  = loss_ledger(sol; T0)       # second law (exergy) — displays as a waterfall
losses_by_kind(L)               # [:aero => …, :pressure => …, …], largest first
ledger_rows(L)                  # (element, kind, W) rows, as printed
```

`energy_balance` fields: `Q_in` (heat sources with positive duty), `Q_out`
(radiators, and heat sources with negative duty), `Q_open` (net enthalpy
carried out by open streams — coolant loops, or the implicit cooler between a
dangling outlet and the seeded inlet of an open-loop model), `P_shaft`,
`imbalance` (Σ enthalpy change of adiabatic elements; ≈ 0), `residual`
(`Q_in − Q_out − Q_open − P_shaft`; ≈ 0), `by_element`.

`LossLedger` fields: `T0`, `Q_in`, `exergy_in` (`(name, W)` per heat source:
enthalpy rise minus `T0` × the entropy rise of the constant-pressure heating),
`P_shaft`, `losses::Vector{Loss}`, `closure` (exergy in − `P_shaft` − Σ losses;
≈ 0 by the Gouy–Stodola identity). `T0` defaults to the coldest available sink:
min of radiator `T_sink`, boundary-stream inlet temperatures, and the seed
inlet temperature.

`Loss(element, kind, exergy)` kinds: `:aero` (compressor/turbine
irreversibility), `:pressure` (total-pressure losses, incl. the pressure-loss
part of heat sources, radiators and heat exchangers), `:heat_transfer`
(finite-ΔT exchange), `:mixing`, `:heat_rejection` (exergy leaving with
radiated or rejected heat, or with open streams — netted per fluid),
`:imbalance` (numerical energy imbalance of an adiabatic element).

New element types join the accounting by defining `element_inlets(el)`,
`element_outlets(el)` (default: `el.inlet`/`el.outlet`), `element_role(el)`
(`:adiabatic`, `:work_in`, `:work_out`, `:heat_in`, `:heat_out`) and, if the
default single-kind entropy-generation entry is not right, `element_losses(el, T0)`.

## Unit conversion helpers

All exported, all with inverses: `R_to_K`/`K_to_R`, `psia_to_Pa`/`Pa_to_psia`,
`lbm_to_kg`/`kg_to_lbm`, `lbps_to_kgps`/`kgps_to_lbps`,
`btulbm_to_Jkg`/`Jkg_to_btulbm`, `btulbmR_to_JkgK`/`JkgK_to_btulbmR`,
`lbmft3_to_kgm3`/`kgm3_to_lbmft3`, `rpm_to_radps`/`radps_to_rpm`,
`hp_to_W`/`W_to_hp`.
