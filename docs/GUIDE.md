# GasCycle.jl — User Guide

GasCycle models closed Brayton cycles (compressor → recuperator → heat source →
turbine, plus cold end) for design-point sizing and map-based off-design
analysis. It is a clean Julia reimplementation of the NPSS cycle-analysis
workflow, with two differentiators: exact design sensitivities via ForwardDiff,
and station-for-station validation against NPSS 3.3.

For every function and element option, see [REFERENCE.md](REFERENCE.md).
Runnable examples live in `examples/`.

## The mental model

A model is built in five steps:

1. **Pick a fluid** — a property backend implementing `h`, `s`, `cp`, `ρ`, …
   For He-Xe work use `HeXe(83.8)` (analytic, no table needed).
2. **Create elements** — `Compressor`, `Turbine`, `HeatExchanger`, etc. Each is
   a small mutable struct; constructor keywords set its parameters and mode.
3. **Wire a `FlowNetwork`** — `add!` the elements, `connect!` the flow chain,
   attach the shaft and recuperator pairing, and pin the entry state with
   `set_state!`.
4. **`solve!`** — the solver propagates states around the loop and Newton-iterates
   whatever is unknown (loop closure states, map operating points, shaft speed).
5. **Read results** — `summary(sol)` for NPSS-style tables, `net_power(sol)`,
   `cycle_efficiency(sol)`, `stations(sol)`, or any element's stored ports.

**Units are SI everywhere** — K, Pa, kg/s, W, J/kg. NPSS files and NASA reports
use English units; convert at the boundary with the exported helpers
(`R_to_K`, `psia_to_Pa`, `lbps_to_kgps`, …).

## A design-point model

```julia
using GasCycle

fluid = HeXe(83.8)                 # He-Xe at M = 83.8 kg/kmol

net   = FlowNetwork()
comp  = Compressor("Comp";  PR=2.5,    η_poly=0.88)
recup = HeatExchanger("Recup"; ε=0.92)
rx    = HeatSource("Reactor";  TtExit=1100.0, dPqP=0.02)
turb  = Turbine("Turb";     mode=:pressure_closure, P_exit=500e3, η_poly=0.90)

add!(net, comp, recup, rx, turb)
connect!(net, comp => recup => rx => turb => comp)   # serial loop
add_hx_pair!(net, recup; hot=turb)                   # turbine exhaust = recup hot side
set_state!(net, comp; Pt=500e3, Tt=400.0, W=10.0, fluid=fluid)

sol = solve!(net)
summary(sol)                       # station + component + cycle tables
net_power(sol), cycle_efficiency(sol)
```

What the solver does here: the recuperator makes the loop circular (its hot
inlet depends on the turbine outlet, which depends on its cold outlet). The
marked back-edge's `[Tt, Pt]` become Newton unknowns; the residual is "the
state I computed equals the state I seeded". Converges in 1–3 iterations.

Design-mode conventions:
- The **compressor PR** is a free design choice; the **turbine** usually runs
  `mode=:pressure_closure` with `P_exit` set so the loop pressure closes.
- The **shaft** in design mode has no residual — `power_balance(shaft)` simply
  reports the net power available to the generator.
- A `HeatSource` in `:fixed_TtExit` mode sets the turbine inlet temperature;
  `:fixed_Q` instead imposes the heat input.

## Closing the loop through the cold end

Instead of pinning the compressor inlet with `set_state!`, you can model heat
rejection so the inlet temperature responds to the operating point:

```julia
cool = HeatExchanger("Cooler"; UA=400.0)             # water-cooled ground test
# ...or... rad = Radiator("Rad"; mode=:fixed_TtExit, TtExit=300.0, T_sink=200.0)

connect_port!(net, recup, :hot_outlet, cool, :hot_inlet)
connect_port!(net, cool, :hot_outlet, comp, :inlet; back_edge=true)
set_boundary!(net, cool, :cold_inlet; Pt=600e3, Tt=290.0, W=1.0, fluid=water)
```

With `back_edge=true` the compressor inlet `[Tt, Pt]` joins the Newton
unknowns and `set_state!` only supplies mass flow, fluid, and the initial
guess. Note the loop pressure needs a free variable to close on — map-based
off-design has one naturally; a design solve with every PR fixed does not
(keep the pinned inlet there, or use `:pressure_closure`).

## Off-design analysis

Off-design replaces fixed design parameters with physical characteristics:
turbomachine performance comes from maps, shaft speed from the power balance,
and heat-exchanger effectiveness from flow-scaled UA. The recipe:

```julia
# 1. Solve the design point (as above), then scale maps through it.  cbase/tbase
#    are native maps — a CompressorMap/TurbineMap from compressor_map()/
#    turbine_map(), or a FunctionMap "script".
s1, s3 = comp.inlet[], turb.inlet[]
cmap = scale_map(cbase; Nc_des=corrected_speed(N_des, s1.Tt),
                 Wc_des=corrected_flow(s1.W, s1.Tt, s1.Pt),
                 PR_des=pressure_ratio(comp), eta_des=0.88)
tmap = scale_map(tbase; Np_des=corrected_speed(N_des, s3.Tt),
                 Wp_des=corrected_flow(s3.W, s3.Tt, s3.Pt),
                 PR_des=pressure_ratio(turb), eta_des=0.90)

# 2. Size the recuperator UA at the design point (ε now responds to flows):
size_UA!(recup)

# 3. Rebuild (or switch) the turbomachines and shaft into off-design mode.
#    Seed the turbine PR at its design expansion ratio:
comp  = Compressor("Comp"; map=cmap, mode=:off_design, η_poly=0.88)
turb  = Turbine("Turb";  PR=pressure_ratio(turb), map=tmap, mode=:off_design, η_poly=0.90)
shaft = Shaft("Shaft"; N=0.97*N_des, mode=:off_design, P_load=P_design)

# 4. Wire the same network shape and solve; then perturb and re-solve:
sol = solve!(net)            # must reproduce the design point
rx.TtExit = 0.9 * 1100.0     # throttle the reactor
sol = solve!(net)
```

How it works: maps are evaluated **forward in their native coordinates** (no
inversion). An off-design compressor owns its R-line coordinate `Rline`; an
off-design turbine owns its expansion ratio `PR` (the map's native input) — each
with a flow-continuity residual (the corrected flow the map passes must equal the
actual flow). An off-design shaft owns its speed `N` with the power-balance
residual (`P_load` is the generator extraction). These join the back-edge states
(which carry loop pressure closure) in one Newton vector, solved with a
trust-region method. An optional `reynolds=ReDesIndex(Re_des)` (or any
`ReynoldsModel`) drives the map's Reynolds-correction tables; the default applies
no correction. A fixed-speed (alternator-locked) study just leaves the shaft in
`:design` mode and sweeps a boundary condition — see
`examples/bru_tit_sweep_offdesign.jl`.

`solve!` picks the mode automatically: if any element declares independent
variables (or free variables are declared, below) it runs the coupled Newton,
otherwise the back-edge Newton.

## Targets and free variables

Most studies are specified by outputs, not inputs: "deliver 10 kW", "hold
design speed", "keep the compressor inlet at 400 K". Declare what the solver
may vary with `add_free!` and what it must achieve with `add_target!`:

```julia
# Design: size the loop mass flow for 10 kW net shaft power
add_free!(net, :W)                       # seed-state mass flow
add_target!(net, net_power, 10e3)
sol = solve!(net)

# Two pairs: flow and compressor PR for power and compressor exit temperature
add_free!(net, comp, :PR; lower=1.2, upper=4.0)
add_target!(net, n -> comp.outlet[].Tt, 520.0; name="T2")

# Off-design: reactor outlet temperature to hold design shaft speed
add_free!(net, reactor, :TtExit)
add_target!(net, n -> shaft.N, N_des)
```

The counts must balance (free variables + element unknowns = targets +
element residuals). The current values are the initial guesses; `ref` sets the
normalization of a free variable whose initial value is far from its typical
magnitude, and `scale` sets the normalization of a target near zero.
`clear_targets!(net)` removes all pairs. `solve!(net; verbose=true)` prints the
converged free variables and targets.

Targets compose with gradients: `ForwardDiff.derivative(η -> W_for_power(η), 0.88)`
differentiates the *solved* free variable, including through off-design map
solves and w.r.t. the target value itself. When differentiating, construct an
element whose field is a free variable in the gradient's number type
(`HeatSource("H"; TtExit = 1100.0 * one(eltype(p)))`) — the error message says
so if you forget.

## Where the energy goes

`energy_balance(sol)` is the first-law check (heat in = shaft power + heat
rejected, to solver tolerance). `loss_ledger(sol)` is the second-law view: of
the exergy the reactor puts into the gas, how much each irreversibility
destroys and how much the cold end throws away. It prints as a waterfall:

```
julia> loss_ledger(sol)
Loss ledger (exergy, dead state T0 = 300.0 K)
  heat input Q = 34.954 kW
  + exergy in  Heater                          24.661 kW   100.0 %
  − Comp  (aero)                                2.441 kW     9.9 %
  − Recup  (heat_transfer)                      1.765 kW    7.16 %
  − Recup  (pressure)                           0.589 kW    2.39 %
  …
  − Recup.hot_outlet  (heat_rejection)          4.267 kW    17.3 %
  = net shaft power                            13.119 kW   53.19 %
  closure residual: 5.82e-11 W
```

The closure residual is an identity, so anything beyond solver tolerance
points at a modeling error. The dead state `T0` defaults to the coldest sink
in the model (radiator sink, coolant inlet, or compressor inlet); pass
`T0 = …` to change it. Ledger values are differentiable, so a loss can be a
target or an optimization objective. `losswaterfall(L)` plots it.

## Bleed flows

Branches use explicit port wiring with `Splitter`/`Mixer`:

```julia
split = Splitter("Bleed"; fracs=[0.98, 0.02])
mix   = Mixer("Return")
connect_port!(net, split, :bleed_outlet, bearing_duct, :inlet)
connect_port!(net, bearing_duct, :outlet, mix, :bleed_inlet)
```

See `examples/bru_10kw.jl` for the BRU bearing-cooling bleed.

## Design sensitivities

All elements are parametric in their numeric type, so ForwardDiff propagates
through the whole cycle — including the Newton solves and real-gas property
inversions. Wrap model construction in a function of the parameters:

```julia
using ForwardDiff
∇W = ForwardDiff.gradient(p -> build_and_solve(p...) , [2.5, 0.92, 1100.0])
dW_dM = ForwardDiff.derivative(M -> power_with_fluid(HeXe(M)), 83.8)
```

See `examples/forwarddiff_sensitivity.jl`.

Because those gradients are exact, they plug straight into a gradient-based
optimizer. `examples/optimize_cycle.jl` maximizes net shaft power over
`[PR_comp, TIT]` (subject to a turbine material limit) by feeding the
ForwardDiff gradient — and the ForwardDiff *Hessian*, which also propagates
through the nested solve — to Optim.jl, then validates the optimum against a
brute-force grid sweep. Optim is not a GasCycle dependency; `] add Optim` into
your environment to run it.

## Plots

With `Plots.jl` loaded, two recipes are available: `tsdiagram(sol)` draws the
cycle on T-s axes with labeled stations, and `mapplot(comp)` draws a
turbomachine's map with its operating point.

## Pitfalls

- `cp` is **not exported** (it clashes with `Base.Filesystem.cp`); use
  `GasCycle.cp(fluid, T, P)` or `import GasCycle: cp`.
- NPSS `effDes` is an **isentropic** efficiency — pass `η_type=:isentropic`
  when replicating NPSS models (the default is `:polytropic`).
- `connect!(net, a => ... => a)` repeats the first element to close the loop;
  the closing edge is stripped — the loop is actually closed by `set_state!`
  or an explicit back-edge.
- Transport properties (`viscosity`, `conductivity`, `prandtl`) are available
  on `NobleGasMixture` (gases) and on `PolynomialLiquid` coolants whose file
  carries μ/k data; the FPT *table* reader doesn't load them.
