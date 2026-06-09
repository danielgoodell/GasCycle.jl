# GasCycle Repository Review

## Improvements Observed

- FPT parsing is much faster. `test/test_thermo.jl` FPT section dropped from
  about `1m07s` previously to about `1.5s`.
- `FluidState` and `Port` are now parametric over `T<:Real`, enabling
  Dual-number propagation for AD-sensitive paths.
- Ideal-gas thermo is now AD-compatible through generic method signatures and
  closed-form inversions.
- `ForwardDiff` is now tested against finite differences for net-power
  sensitivity.
- The network model is more general. It now uses explicit `PortEdge`s and
  supports splitter/mixer branch paths.
- Splitter and mixer elements were added and are tested.
- BRU validation is improved with a bleed branch and closer system-level
  comparison.

## Remaining Findings

- **High:** Off-design map residuals are still ineffective. `Wc_map` is set to
  `Wc_act`, then residuals compute `Wc_map - Wc_act`, so the residual is
  identically zero.
- **High:** Pressure-closure turbine reporting regressed. `PR_eff` is computed
  locally but `el.PR` is intentionally left unchanged, so examples print stale
  `2.0` even though the state update uses the correct pressure ratio. This shows
  up in `recuperated_brayton.jl` and `bru_tit_sweep.jl`.
- **Medium:** `Pkg.test()` still fails. `test/runtests.jl` imports `Pkg`, but
  `Pkg` is not declared for the isolated test environment. Direct test execution
  passes.
- **Medium:** `NonlinearSolve` is now a direct dependency, but the solver still
  uses the custom forward-difference Newton implementation. This adds large
  precompile cost without current benefit.
- **Medium:** Fixed-point convergence still only checks outlet temperatures. It
  does not check pressure, flow, back-edge residuals, power balance, or graph
  closure.
- **Medium:** `one_pass!` does not error if some elements remain unprocessed due
  to missing edges or cyclic dependencies. That can hide malformed networks.
- **Medium:** FPT AD is only partial. Forward property calls accept generic
  inputs, but `T_from_h` and `T_from_s` still use bisection and comments say
  derivatives are not propagated.
- **Low:** Heat exchanger remains `cp * ΔT` based rather than enthalpy based, so
  real-gas/large-temperature-span accuracy is still approximate.

## Validation

- `julia --project=. test/runtests.jl`: passed.
- `julia --project=. test/test_thermo.jl`: passed, FPT section about `1.5s`.
- `julia --project=. test/test_elements.jl`: passed, including ForwardDiff
  gradient test.
- `julia --project=. examples/recuperated_brayton.jl`: passed in about `9s`.
- `julia --project=. examples/forwarddiff_sensitivity.jl`: passed.
- `julia --project=. examples/bru_10kw.jl`: passed.
- `julia --project=. -e 'using Pkg; Pkg.test()'`: failed due `Pkg` import in
  test harness.

## Suggested Next Steps

1. Fix `Pkg.test()` by removing `Pkg.activate` from test files/runtests or
   declaring proper test dependencies.
2. Fix turbine pressure-closure reporting by storing a primal `PR` when safe, or
   adding a `pressure_ratio(turb)` accessor.
3. Fix off-design map residual formulation before investing further in map-based
   operation.
4. Either actually migrate the solver to `NonlinearSolve.jl`, or remove it from
   direct dependencies to avoid heavy precompile cost.
5. Add network validation: every non-seed inlet satisfied, every added element
   processed, no unintended cycles, back-edge residuals measurable.
6. Add end-to-end tests for the BRU bleed network and TIT sweep PR reporting.
