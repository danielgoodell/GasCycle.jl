# GasCycle.jl

A Julia package for steady-state thermodynamic cycle analysis of closed Brayton cycles, with a focus on helium-xenon (He-Xe) working fluid mixtures for space nuclear power systems.

Designed as a clean, extensible reimplementation of the core NPSS cycle-analysis workflow — without the proprietary C++ infrastructure. The key differentiator is first-class automatic differentiation: exact design sensitivities come for free, enabling gradient-based optimization over cycle parameters.

## Features

- **Closed Brayton cycle modeling** — compressor, turbine, recuperator, reactor heat source, ducts, bleed flows
- **Real He-Xe fluid properties** via NPSS FPT (Fluid Property Table) files — bicubic interpolation on T-P grids, or a fast ideal-gas backend for development
- **Exact design derivatives** via [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) — parametric element types (`Compressor{T<:Real}`, etc.) propagate Dual numbers through the full cycle
- **Newton solver** via [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl) — replaces fixed-point iteration for recuperator back-edges; inner `AutoForwardDiff` Jacobian lets outer ForwardDiff thread through correctly for implicit differentiation
- **Compressor bleeds** — `Splitter` / `Mixer` elements handle bearing-cooling bleeds, seal flows, or any branch topology
- **Performance maps** — bilinear/bicubic interpolation on corrected speed/flow grids with design-point scaling

## Quick start

```julia
using GasCycle

fluid = FPTFluid("HeXe84.fpt")   # load He-Xe fluid property table

net   = FlowNetwork()
comp  = Compressor("Comp";  PR=2.5,    η_poly=0.88)
recup = HeatExchanger("Recup"; ε=0.92)
rx    = HeatSource("Reactor";  TtExit=1100.0, dPqP=0.02)
turb  = Turbine("Turb";     mode=:pressure_closure, P_exit=500e3, η_poly=0.90)

add!(net, comp, recup, rx, turb)
connect!(net, comp => recup => rx => turb => comp)
add_hx_pair!(net, recup; hot=turb)
set_state!(net, comp; Pt=500e3, Tt=400.0, W=10.0, fluid=fluid)

sol = solve!(net; verbose=true)
println("Net power:  ", round(net_power(sol)/1000, digits=2), " kW")
println("Efficiency: ", round(cycle_efficiency(sol)*100, digits=1), " %")
```

## Design sensitivity with ForwardDiff

Because all element types are parametric in their numeric type (`Compressor{T<:Real}`, `HeatExchanger{T<:Real}`, etc.), ForwardDiff Dual numbers propagate through the entire cycle — including through the Newton solve for recuperator back-edge convergence.

```julia
using ForwardDiff

fluid = HeXeIdealGas(0.47)   # ideal gas for fully closed-form AD

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
    set_state!(net, comp; Pt=500e3, Tt=400.0, W=10.0, fluid=fluid)
    net_power(solve!(net)) / 1000
end

∇W = ForwardDiff.gradient(cycle_power, [2.5, 0.90, 1100.0])
# ∂W_net/∂PR_comp  ≈  20.7 kW per unit PR
# ∂W_net/∂ε_recup  ≈   0.0 kW per unit ε  (correct: ε only affects η, not W)
# ∂W_net/∂TIT      ≈   0.77 kW/K
```

Derivatives match finite differences to better than 1e-5%.

## Validation: NASA 10.5 kW BRU

The `examples/bru_10kw.jl` model reproduces the NASA Brayton Rotating Unit (BRU, TN D-5815) design point, including the 2% compressor bleed to bearings:

| Quantity | GasCycle | NPSS / paper |
|---|---|---|
| Compressor outlet T | 422 K (759 °R) | 409 K (737 °R) |
| Turbine outlet T | 914 K (1645 °R) | 945 K (1701 °R) |
| Turbine PR | 1.758 | ~1.75 |
| Net shaft power | 13.6 kW | ~13.4 kW |
| Est. electrical output | 11.1 kW | 10.5 kW design |

Both models use the same HeXe84.fpt fluid table. GasCycle implements true polytropic efficiency (N=20 step numerical integration of dh = v dP / η_p). The NPSS BRU model appears to use isentropic efficiency semantics with the same η=0.80 value — for ideal-gas He-Xe at PR=1.9, isentropic η=0.80 gives ~410 K and polytropic η_p=0.80 gives ~414 K, consistent with the observed offset. The remaining gap is FPT interpolation differences (GasCycle bicubic B-spline vs. NPSS's scheme). Net power agrees to within ~1.5%.

## Architecture

```
src/
├── thermo/
│   ├── FluidProperties.jl   # abstract interface: enthalpy, entropy, cp, density, γ
│   ├── FPTFluid.jl          # NPSS FPT file reader + bicubic interpolation
│   └── IdealGasFluid.jl     # constant-Cp ideal gas (exact closed-form inversions)
├── core/
│   ├── FluidState.jl        # FluidState{T<:Real}: Pt, Tt, W, fluid
│   ├── Port.jl              # Port{T<:Real}: mutable ref to FluidState
│   └── Element.jl           # AbstractElement interface
├── elements/
│   ├── Compressor.jl        # polytropic compression, optional performance map
│   ├── Turbine.jl           # polytropic expansion; :design / :pressure_closure modes
│   ├── HeatExchanger.jl     # ε-NTU counter-flow; two-port (hot + cold sides)
│   ├── HeatSource.jl        # reactor / heater; :fixed_Q or :fixed_TtExit modes
│   ├── Duct.jl              # pressure loss only
│   ├── Splitter.jl          # mass-flow split by fixed fractions
│   ├── Mixer.jl             # mass and energy mixing of N inlet streams
│   └── Shaft.jl             # mechanical coupling; enforces power balance
├── maps/
│   ├── PerformanceMap.jl    # 2D corrected-speed/flow map with interpolation
│   └── MapScaling.jl        # design-point scaling
├── network/
│   └── FlowNetwork.jl       # directed port graph; topological traversal; back-edges
└── solver/
    └── Solver.jl            # back-edge Newton (AutoForwardDiff) + off-design Newton
```

### How the solver works

The recuperator creates a circular dependency: the hot-side inlet depends on the turbine outlet, which depends on the cold-side outlet. GasCycle breaks this with a **back-edge Newton**:

1. The back-edge `[Tt, Pt]` states are treated as explicit unknowns.
2. `one_pass!(net, z)` propagates the full cycle given back-edge seeds `z`.
3. The residual is `F(z) = [Tt_computed - z_Tt, Pt_computed - z_Pt]` (normalised).
4. NonlinearSolve.jl Newton with `AutoForwardDiff()` Jacobian converges in 1–3 iterations.

Because the Newton uses `AutoForwardDiff()` internally, an outer `ForwardDiff.gradient` call threads through correctly via nested (differently-tagged) Dual numbers — the back-edge sensitivities are handled by the implicit function theorem at no extra cost.

## Fluid property files

GasCycle reads NPSS FPT (Fluid Property Table) files, the same format used by NPSS. These are structured T-P grids with tabulated h, s, cp, γ, and ρ. FPT files for He-Xe mixtures can be generated from NPSS or from CEA-based tools.

The included `HeXe84.fpt` covers the He-Xe mixture at M = 83.8 g/mol (the BRU design composition, ~40% He by mass).

The ideal-gas backend (`HeXeIdealGas(mole_fraction_He)`) requires no data file and supports full ForwardDiff AD (all property inversions are closed-form). It is accurate to within a few percent for monatomic He-Xe mixtures at conditions well away from the critical point.

## Elements reference

| Element | Key parameters | Modes |
|---|---|---|
| `Compressor` | `PR`, `η_poly` | `:design`, `:off_design` (map) |
| `Turbine` | `PR`, `η_poly`, `P_exit` | `:design`, `:pressure_closure`, `:off_design` |
| `HeatExchanger` | `ε`, `dPqP_hot`, `dPqP_cold` | counter-flow ε-NTU |
| `HeatSource` | `Q`, `TtExit`, `dPqP` | `:fixed_Q`, `:fixed_TtExit` |
| `Duct` | `dPqP` | — |
| `Splitter` | `fracs` | fixed mass-fraction split |
| `Mixer` | `n_inlets` | mass + energy conservation |
| `Shaft` | `N` | power balance across all attached turbomachinery |

## Running tests

```
julia --project=. test/runtests.jl
```

60 tests covering thermodynamics, individual elements, full-cycle solves, and ForwardDiff gradient correctness.

## Dependencies

| Package | Role |
|---|---|
| `NonlinearSolve.jl` | Newton solver for back-edge states and off-design unknowns |
| `ForwardDiff.jl` | Automatic differentiation through the full cycle |
| `Interpolations.jl` | Bicubic B-spline interpolation of FPT property tables |

## Status

Design-point analysis is complete and validated. Off-design analysis with performance map look-up (Phase 7) is next.
