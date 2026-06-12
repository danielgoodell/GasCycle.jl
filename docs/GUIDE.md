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
# 1. Solve the design point (as above), then scale maps through it:
s1, s3 = comp.inlet[], turb.inlet[]
cmap = scale_map(cbase; Nc_des=corrected_speed(N_des, s1.Tt),
                 Wc_des=corrected_flow(s1.W, s1.Tt, s1.Pt),
                 PR_des=pressure_ratio(comp), eta_des=0.88)
tmap = scale_map(tbase; Nc_des=corrected_speed(N_des, s3.Tt),
                 Wc_des=corrected_flow(s3.W, s3.Tt, s3.Pt),
                 PR_des=pressure_ratio(turb), eta_des=0.90)

# 2. Size the recuperator UA at the design point (ε now responds to flows):
size_UA!(recup)

# 3. Rebuild (or switch) the turbomachines and shaft into off-design mode:
comp  = Compressor("Comp"; map=cmap, mode=:off_design, η_poly=0.88)
turb  = Turbine("Turb";  map=tmap, mode=:off_design, η_poly=0.90)
shaft = Shaft("Shaft"; N=0.97*N_des, mode=:off_design, P_load=P_design)

# 4. Wire the same network shape and solve; then perturb and re-solve:
sol = solve!(net)            # must reproduce the design point
rx.TtExit = 0.9 * 1100.0     # throttle the reactor
sol = solve!(net)
```

How it works: each off-design turbomachine owns one unknown (its map flow
coordinate `Wc_map`) with a flow-continuity residual; an off-design shaft owns
its speed `N` with the power-balance residual (`P_load` is the generator
extraction). These join the back-edge states in one Newton vector, solved with
a trust-region method. A fixed-speed (alternator-locked) study just leaves the
shaft in `:design` mode and sweeps a boundary condition — see
`examples/bru_tit_sweep_offdesign.jl`.

`solve!` picks the mode automatically: if any element declares independent
variables it runs the off-design Newton, otherwise the back-edge Newton.

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
- Transport properties (`viscosity`, `conductivity`, `prandtl`) are only
  available on the `NobleGasMixture` backend; FPT tables don't carry them.
