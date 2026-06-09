# GasCycle Repository Review

## Assessment

This is a reasonable early Brayton-cycle prototype, but it is not yet close to
NPSS core functionality. It currently works best as a sequential design-point
cycle calculator with a special-case recuperator, fixed mass flow, fixed or
derived pressure ratios, and simple component models.

## Findings

- **High:** Off-design map residuals are ineffective. In
  `src/elements/Compressor.jl` and `src/elements/Turbine.jl`, `Wc_map` is set
  equal to actual corrected flow, then residuals compute `Wc_map - Wc_act`.
  This makes the residual identically zero, so off-design map solving is not
  actually constrained.
- **High:** The network is not a general NPSS-like graph. `FlowNetwork.one_pass!`
  resets each pass from `initial_state` and walks `flow_order`; `connections`
  are stored but not really used for graph propagation. This supports one serial
  loop plus a special heat-exchanger path, not arbitrary NPSS-style elements,
  ports, splitters, mixers, bleeds, or nested assemblies.
- **Medium:** Solver convergence is based only on outlet temperatures in
  fixed-point mode. It ignores pressure, flow, power balance, and residuals. A
  cycle could report convergence while pressure/flow closure or mechanical
  constraints are physically wrong.
- **Medium:** The code claims to be ForwardDiff-friendly, but `FluidState`
  forces `Float64`, and most thermo APIs dispatch on `Float64`. That blocks dual
  numbers and makes AD-based Newton solves difficult.
- **Medium:** FPT property lookup silently clamps outside table bounds. That can
  hide invalid operating points and produce plausible but wrong answers. FPT
  inverse lookups also fall back to bisection, making `test/test_thermo.jl` take
  about 67 seconds.
- **Medium:** Heat exchanger modeling is simplified for real-gas usage. It uses
  inlet `cp` and temperature deltas instead of enthalpy differences, so large
  temperature spans or non-constant-property fluids will be approximate.
- **Low:** `Pkg.test()` fails because there is no `test/runtests.jl`. The
  individual test files pass when run directly.
- **Low:** Example wording is misleading: `examples/recuperated_brayton.jl` says
  recuperation should raise efficiency above "Ideal Brayton eta", but the
  printed reference `1 - T_in/TIT` is not the ideal Brayton efficiency for the
  specified pressure ratio and is much higher than the actual recuperated
  result.

## Validation Run

- `julia --project=. test/test_elements.jl`: passed.
- `julia --project=. test/test_thermo.jl`: passed, but FPT section took about
  `1m07s`.
- `julia --project=. examples/simple_brayton.jl`: ran, efficiency `19.07%`.
- `julia --project=. examples/recuperated_brayton.jl`: ran, efficiency `29.46%`.
- `julia --project=. -e 'using Pkg; Pkg.test()'`: failed due missing
  `test/runtests.jl`.

## Suggestions

1. Add `test/runtests.jl` and make CI run the same command users will run.
2. Decide whether this is a design-point cycle calculator or an NPSS-like
   equation system. If NPSS-like is the goal, move toward explicit variables,
   residuals, balances, and graph connectivity rather than sequential
   propagation.
3. Redesign off-design turbomachinery maps around actual independent variables:
   shaft speed, corrected flow, map R-line/flow parameter, pressure ratio,
   efficiency, and shaft power balance.
4. Make `FluidState` and thermo methods numeric-type generic if AD or robust
   nonlinear solving is expected.
5. Replace silent FPT clamping with configurable behavior: `:error`, `:warn`, or
   `:clamp`.
6. Add end-to-end tests for simple closed loop, recuperated closed loop, failed
   convergence, map off-design, and pressure closure.
7. Consider using `NonlinearSolve.jl` or `NLsolve.jl` before implementing a
   custom Newton solver. NPSS-like systems need scaling, bounds, damping,
   diagnostics, and robust failure reporting.
