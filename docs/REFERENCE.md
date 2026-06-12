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
| `ConstantPropertyLiquid(; cp, rho, name="liquid")` | constant-property coolant | also `ConstantPropertyLiquid(path)` for function-style FPT files (`data/H2O.fpt`, `data/Oil.fpt`) |

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
| `viscosity(fl, T, P)` | μ [Pa·s] — `NobleGasMixture` only |
| `conductivity(fl, T, P)` | k [W/(m·K)] — `NobleGasMixture` only |
| `prandtl(fl, T, P)` | Pr [-] — `NobleGasMixture` only |

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
Compressor(name; PR=2.0, η_poly=0.87, η_type=:polytropic, map=nothing, mode=:design)
```
- `mode=:design` — `PR` fixed by the user.
- `mode=:off_design` — PR and η interpolated from `map` at `(Nc, Wc_map)`;
  `Wc_map` is a solver unknown with a flow-continuity residual.
- `η_type` — `:polytropic` (small-stage) or `:isentropic` (NPSS `effDes`).

### Turbine
```julia
Turbine(name; PR=2.0, η_poly=0.90, η_type=:polytropic, map=nothing,
        mode=:design, P_exit=101325.0)
```
- `mode=:design` — `PR` fixed.
- `mode=:pressure_closure` — PR computed so the exhaust hits `P_exit`
  (standard for closed-cycle design solves).
- `mode=:off_design` — as Compressor.

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
sol = solve!(net; tol=1e-6, maxiter=100, verbose=false)   # -> SolveResult
```

Mode is automatic: with no element unknowns, a back-edge Newton
(NonlinearSolve `NewtonRaphson`, ForwardDiff Jacobian) on the back-edge
`[Tt, Pt]` states; with unknowns (off-design `Wc_map`s, shaft `N`), everything
joins one normalized vector solved by `TrustRegion`. `tol` is on normalized
residuals. Non-convergence warns and returns `status = :failed`.

`SolveResult` fields: `status` (`:success`/`:failed`), `iterations`,
`residual_norm`, `net`. `sol["Comp"]` looks an element up by name.

```julia
net_power(sol)          # Σ turbine − Σ compressor power [W]
cycle_efficiency(sol)   # net power / Σ heat input
```

## Performance maps

```julia
m = PerformanceMap(Nc_axis, Wc_axis, PR_grid, eta_grid; bounds=:error)
query(m, Nc, Wc)               # -> (PR, η), bilinear
corrected_speed(N, Tt)          # N / √(Tt/288.15)
corrected_flow(W, Tt, Pt)       # W·√(Tt/288.15) / (Pt/101325)
scale_map(base; Nc_des, Wc_des, PR_des, eta_des,
          Nc_ref=center, Wc_ref=center, bounds=base.bounds)
```

`bounds` ∈ `:error` (DomainError outside the map), `:clamp`, `:warn`.
`scale_map` makes the base map pass through the design point: axes scale by
`des/ref`, PR scales as `1 + (PR−1)·s`, η multiplicatively. Pick `Nc_ref`/`Wc_ref`
off the grid nodes so the design point lands inside a smooth cell.

NPSS `.map` (NEO Table) files:

```julia
tables = read_npss_map("compressor.map")       # Dict{String, NPSSMapTable}
m = to_performance_map(tables; alpha=NaN, flow="TB_Wc", pr="TB_PR",
                       eff="TB_eff", nWc=25)   # Rline lines → rectangular map
```

Run the result through `scale_map` at your design point (this also makes the
lbm/s flow units in NPSS maps irrelevant).

## Output

```julia
summary(sol)                    # print NPSS-style station/component/cycle tables
sprint(summary, sol)            # same, captured as a String
stations(sol; branches=true)    # Vector of "Element.port" => FluidState, flow order
tsdiagram(sol)                  # Plots.jl recipe: T-s diagram with stations
mapplot(m); mapplot(comp)       # map speed lines; with operating point if element
```

`stations` walks the main loop in physical order (through HX hot sides via
back-edges), then bleed branches.

## Unit conversion helpers

All exported, all with inverses: `R_to_K`/`K_to_R`, `psia_to_Pa`/`Pa_to_psia`,
`lbm_to_kg`/`kg_to_lbm`, `lbps_to_kgps`/`kgps_to_lbps`,
`btulbm_to_Jkg`/`Jkg_to_btulbm`, `btulbmR_to_JkgK`/`JkgK_to_btulbmR`,
`lbmft3_to_kgm3`/`kgm3_to_lbmft3`, `rpm_to_radps`/`radps_to_rpm`,
`hp_to_W`/`W_to_hp`.
