# GasCycle.jl

A Julia package for thermodynamic cycle analysis of closed Brayton cycles, with a focus on helium-xenon (He-Xe) working fluid mixtures for space nuclear power systems.

**Documentation:** [user guide](docs/GUIDE.md) · [API reference](docs/REFERENCE.md)

Designed as a clean, extensible reimplementation of the core NPSS cycle-analysis workflow — without the proprietary C++ infrastructure. Two key differentiators:

- **First-class automatic differentiation** — exact design sensitivities come for free, enabling gradient-based optimization over any cycle parameter, including the He-Xe mixture ratio itself.
- **Validated against NPSS** — running the same model on the same fluid table, GasCycle matches the NPSS 3.3 output station-for-station to print precision (see `validation/`).

## Features

- **Closed Brayton cycle modeling** — compressor, turbine, recuperator, reactor heat source, ducts, bleed flows (`Splitter`/`Mixer`), and a full cold end (radiator or coolant-loop cooler) so the loop closes physically
- **Design and off-design analysis** — design-point sizing, then map-based off-design with shaft power balance (`P_load` generator extraction), validated with constant-speed TIT sweeps down to the self-sustain threshold; heat-exchanger effectiveness responds to off-design flows via flow-scaled UA (`size_UA!` after the design solve)
- **Four fluid property backends** (see table below) — including a direct analytic noble-gas backend that needs no property table at all
- **Transport properties** — μ, k, Pr for all noble gases and binary mixtures, validated against experiment (groundwork for heat-exchanger sizing)
- **Exact design derivatives** via [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) — parametric element types (`Compressor{T<:Real}`, etc.) propagate Dual numbers through the full cycle, including through the Newton solves (implicit function theorem via nested Duals)
- **Newton solver** via [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl) — back-edge states (recuperator, cold-end loop closure) and off-design unknowns (map coordinates, shaft speed) solved together
- **Performance maps** — corrected speed/flow maps with design-point scaling, plus a reader for NPSS `.map` (NEO Table) files
- **NPSS-style output** — `summary(sol)` prints station/component/cycle tables in physical flow order; `stations(sol)` returns the data; Plots recipes for T-s diagrams and map operating points
- **Benchmark suite** — `benchmarks/runbenchmarks.jl`, including a three-backend comparison on an identical model

## Quick start

```julia
using GasCycle

fluid = HeXe(83.8)   # He-Xe at M = 83.8 kg/kmol — analytic, no table needed

net   = FlowNetwork()
comp  = Compressor("Comp";  PR=2.5,    η_poly=0.88)
recup = HeatExchanger("Recup"; ε=0.92)
rx    = HeatSource("Reactor";  TtExit=1100.0, dPqP=0.02)
turb  = Turbine("Turb";     mode=:pressure_closure, P_exit=500e3, η_poly=0.90)

add!(net, comp, recup, rx, turb)
connect!(net, comp => recup => rx => turb => comp)
add_hx_pair!(net, recup; hot=turb)
set_state!(net, comp; Pt=500e3, Tt=400.0, W=10.0, fluid=fluid)

sol = solve!(net)
println("Net power:  ", round(net_power(sol)/1000, digits=2), " kW")   # 234.58 kW
println("Efficiency: ", round(cycle_efficiency(sol)*100, digits=1), " %")  # 29.8 %
summary(sol)   # NPSS-style station + component + cycle tables
```

## Fluid property backends

| Backend | Source | Real-gas | Transport (μ, k, Pr) | AD |
|---|---|---|---|---|
| `NobleGasMixture` / `HeXe(M)` / `NobleGasFluid(gas)` | analytic — virial EOS + corresponding-states transport (Tournier/El-Genk AIAA 2006-4154) | yes | yes | full, incl. d/d(mixture ratio) |
| `FPTFluid("data/HeXe84.fpt")` | NPSS FPT table files (T-P grids; forward + inverse tables) | yes (whatever the table encodes) | no (the reader skips the μ/k/Pr tables present in HeXe84.fpt) | full |
| `IdealGasFluid` / `HeXeIdealGas(x_He)` | constant-cp monatomic ideal gas | no | no | full (closed form) |
| `ConstantPropertyLiquid` | constant cp/ρ coolants; reads function-style FPT files (`data/H2O.fpt`, `data/Oil.fpt`) | — | no | full |

The direct `NobleGasMixture` backend covers all five noble gases (He, Ne, Ar, Kr, Xe — `NobleGasFluid(ARGON)`, etc.) and their ten binary mixtures at any composition, with no table generation step. Thermodynamics come from the virial EOS (Eqs. 8–22 of the paper), transport from data-fitted pair correlations (Eqs. 2–7, 23–36); mixture Prandtl numbers match Taylor's 1988 He-Xe experiments to ≤1.1%. It is the convenient default; scalar calls are allocation-free and a full design solve costs ~1.4 ms vs ~0.3 ms with FPT tables (see `benchmarks/`).

Keep the FPT reader for cross-validation: reading the literal file the collaborator feeds NPSS is the cleanest apples-to-apples comparison available. The included `data/HeXe84.fpt` covers the BRU design composition (audited M ≈ 84.07 kg/kmol; provenance and a fidelity audit are in `data/README.md`).

## Design sensitivity with ForwardDiff

Because all element types are parametric in their numeric type, ForwardDiff Dual numbers propagate through the entire cycle — including through the Newton solve for back-edge convergence, and including through the real-gas property evaluations and their Newton inversions.

```julia
using ForwardDiff

function cycle_power(params)
    PR, ε, TIT = params
    net   = FlowNetwork()
    comp  = Compressor("C";  PR=PR,    η_poly=0.88)
    recup = HeatExchanger("R"; ε=ε)
    rx    = HeatSource("H";  TtExit=TIT, dPqP=0.02)
    turb  = Turbine("T";     mode=:pressure_closure, P_exit=500e3, η_poly=0.90)
    add!(net, comp, recup, rx, turb)
    connect!(net, comp => recup => rx => turb => comp)
    add_hx_pair!(net, recup; hot=turb)
    set_state!(net, comp; Pt=500e3, Tt=400.0, W=10.0, fluid=HeXe(83.8))
    net_power(solve!(net)) / 1000
end

∇W = ForwardDiff.gradient(cycle_power, [2.5, 0.90, 1100.0])
# ∂W_net/∂PR_comp  ≈  13.0  kW per unit PR
# ∂W_net/∂ε_recup  ≈   0.0  kW per unit ε  (correct: ε only affects η, not W)
# ∂W_net/∂TIT      ≈   0.68 kW/K
```

The mixture composition itself is differentiable — useful for working-fluid trade studies:

```julia
dW_dM = ForwardDiff.derivative(M -> cycle_power_with_fluid(HeXe(M)), 83.8)
# ≈ −2.75 kW per (kg/kmol) of mixture molecular weight
```

Derivatives match finite differences to better than 1e-5%.

## Validation

Two layers, documented in `validation/PLAN.md` and `validation/RESULTS.md`:

**NPSS cross-validation (primary).** GasCycle runs a replica of the collaborator's `BRU3.mdl` (`reference/`) on the identical `HeXe84.fpt` table, with NPSS's actual solver setup (tear-based, floating TIT and turbine PR). Result: every station temperature, pressure, and power matches the NPSS 3.3 output (`reference/HeXe.out`) within print precision — TIT to 0.025 °R with each code solving its own shaft balance. All historical discrepancies were attributed exactly (entropy-interpolation artifact — found in both codes, fixed here with pressure-detrended interpolation; isentropic `effDes` semantics; NPSS's workless-bypass bleed; kW-tagged HPX; Oil.fpt's cp).

**NASA 10.5 kW BRU (TN D-5815).** `examples/bru_10kw.jl` reproduces the Brayton Rotating Unit design point including the 2% bearing-cooling bleed: compressor outlet matches at print precision, no-bleed turbine outlet design value to 0.2 °R, ~10.6 kW estimated electrical vs 10.5 kW design.

The analytic `NobleGasMixture` backend is additionally validated against the El-Genk paper's spot values, Johnson NASA/CR-2006-214394 transport tables, Taylor's experimental Prandtl data, and `HeXe84.fpt` itself at cycle conditions (`test/test_noblegas.jl`, `validation/transport_notes.md`).

## Architecture

```
src/
├── units.jl                 # °R/K, psia/Pa, BTU, hp, … conversion helpers
├── thermo/
│   ├── FluidProperties.jl   # abstract interface: h, s, cp, ρ, γ, μ, k, Pr + inversions
│   ├── IdealGasFluid.jl     # constant-cp ideal gas (exact closed-form inversions)
│   ├── FPTFluid.jl          # NPSS FPT reader; forward + inverse tables, detrended s
│   ├── ConstantPropertyLiquid.jl  # coolants; reads function-style FPT files
│   └── NobleGasMixture.jl   # analytic virial EOS + transport, noble-gas binaries
├── core/
│   ├── FluidState.jl        # FluidState{T<:Real}: Pt, Tt, W, fluid
│   ├── Port.jl              # Port{T<:Real}: mutable ref to FluidState
│   └── Element.jl           # AbstractElement interface; polytropic stepping
├── elements/
│   ├── Compressor.jl        # compression; :design / :off_design (map) modes
│   ├── Turbine.jl           # expansion; :design / :pressure_closure / :off_design
│   ├── HeatExchanger.jl     # ε-NTU counter-flow; fixed-ε, :UA, or :scaled_UA
│   ├── HeatSource.jl        # reactor / heater; :fixed_Q or :fixed_TtExit
│   ├── Radiator.jl          # segmented σεA(T⁴−Tsink⁴); :fixed_area / :fixed_TtExit
│   ├── Duct.jl              # pressure loss only
│   ├── Splitter.jl          # mass-flow split by fixed fractions
│   ├── Mixer.jl             # mass and energy mixing of N inlet streams
│   └── Shaft.jl             # power balance; off-design speed + P_load extraction
├── maps/
│   ├── PerformanceMap.jl    # corrected speed/flow map with interpolation
│   ├── MapScaling.jl        # design-point scaling
│   └── NPSSMapReader.jl     # NPSS .map (NEO Table) → PerformanceMap
├── network/
│   └── FlowNetwork.jl       # directed port graph; flow plan; back-edges; boundaries
├── solver/
│   └── Solver.jl            # back-edge + off-design Newton (NonlinearSolve)
└── output/
    ├── Summary.jl           # summary(sol) station/component tables; stations(sol)
    └── PlotRecipes.jl       # tsdiagram, mapplot (Plots.jl recipes)
```

### How the solver works

The recuperator creates a circular dependency: the hot-side inlet depends on the turbine outlet, which depends on the cold-side outlet. Closing the loop through a cooler or radiator (`connect_port!(...; back_edge=true)`) adds another. GasCycle treats these as **back-edge Newton** unknowns:

1. The back-edge `[Tt, Pt]` states are explicit unknowns.
2. `one_pass!(net, z)` propagates the full cycle given back-edge seeds `z`.
3. The residual is `F(z) = [Tt_computed − z_Tt, Pt_computed − z_Pt]` (normalised).
4. NonlinearSolve.jl Newton with `AutoForwardDiff()` Jacobian converges in 1–3 iterations.

In off-design mode the map coordinates of each turbomachine (flow-continuity residual onto the map) and the shaft speed (power-balance residual) join the same Newton vector, with TrustRegion for robustness across map-cell kinks.

Because the Newton uses `AutoForwardDiff()` internally, an outer `ForwardDiff.gradient` call threads through correctly via nested (differently-tagged) Dual numbers — the back-edge sensitivities are handled by the implicit function theorem at no extra cost.

## Elements reference

| Element | Key parameters | Modes |
|---|---|---|
| `Compressor` | `PR`, `η_poly`, `map` | `:design`, `:off_design` (map) |
| `Turbine` | `PR`, `η_poly`, `P_exit`, `map` | `:design`, `:pressure_closure`, `:off_design` |
| `HeatExchanger` | `ε` or `UA`, `UA_exp`, `dPqP_hot`, `dPqP_cold` | fixed-ε, `:UA` (ε-NTU each pass), `:scaled_UA` (UA ∝ W^0.8 from design point via `size_UA!`) |
| `HeatSource` | `Q`, `TtExit`, `dPqP` | `:fixed_Q`, `:fixed_TtExit` |
| `Radiator` | `A`, `emissivity`, `T_sink`, `TtExit`, `N_seg` | `:fixed_TtExit` (sizing), `:fixed_area` (off-design) |
| `Duct` | `dPqP` | — |
| `Splitter` | `fracs` | fixed mass-fraction split |
| `Mixer` | `n_inlets` | mass + energy conservation |
| `Shaft` | `N`, `P_load` | `:design` (balance check), `:off_design` (N solved) |

`η_type=:isentropic` or `:polytropic` selects the efficiency semantics on `Compressor`/`Turbine` (NPSS `effDes` is isentropic).

## Running tests and benchmarks

```
julia --project -e 'using Pkg; Pkg.test()'    # ~530 tests
julia benchmarks/runbenchmarks.jl             # BenchmarkTools suite
```

Tests cover thermodynamics and transport (with literature oracles), individual elements, full-cycle design and off-design solves, the NPSS map reader, plot recipes, and ForwardDiff gradient correctness.

## Dependencies

| Package | Role |
|---|---|
| `NonlinearSolve.jl` | Newton solver for back-edge states and off-design unknowns |
| `ForwardDiff.jl` | Automatic differentiation through the full cycle |
| `Interpolations.jl` | B-spline interpolation of FPT property tables |
| `RecipesBase.jl` | Plots.jl recipes without a hard Plots dependency |

## Status

Design-point and map-based off-design analysis are complete and validated against NPSS to print precision. The direct noble-gas property backend (thermo + transport) replaces FPT tables for He-Xe work. Off-design heat-exchanger effectiveness follows from flow-scaled UA sized at the design point. Next on the roadmap (`ROADMAP.md`): inventory (charge-pressure) control with per-component volume bookkeeping, and transients — shaft dynamics first, then thermal capacitance.
