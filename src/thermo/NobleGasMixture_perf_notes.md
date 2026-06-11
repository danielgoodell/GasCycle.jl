# NobleGasMixture Performance Notes

Review date: 2026-06-10

## Current baseline

- `NobleGasMixture.jl` thermodynamic scalar calls are type-stable and allocation-free.
- He-Xe direct scalar thermodynamic calls are roughly `330 ns/property`; FPT table lookups are roughly `167 ns/property`, but allocate per lookup.
- Direct inverse calls are faster than current FPT fallback inversions in isolation:
  - `T_from_h`: about `1.8 us`
  - `T_from_s`: about `2.3 us`
  - FPT fallback inversions were about `4.3 us` and `5.5 us`
- A BRU-shaped primal `one_pass!` was faster with NobleGas than FPT, but the design-solver ForwardDiff back-edge Jacobian was slower with NobleGas.

## Main reminders

- ~~Fix `HeXe(M_molar)` for AD.~~ DONE 2026-06-10 (name built from `ForwardDiff.value`).
- Add shared helpers for inverse solves, such as `_h_cp(fp,T,P)` and `_s_cp(fp,T,P)`, so each Newton step computes virial coefficients and density once.
- Consider a bundled state API, such as `thermo_state(fp,T,P)`, returning `rho`, `h`, `s`, `cp`, and `gamma` together. Element code often asks for multiple properties at the same state.
- Precompute mixture constants in the `NobleGasMixture` constructor: `x2` and composition weights (transport pair constants are now resolved in the constructor, 2026-06-10).
- Keep the full 18-term LJ virial series for broad validity unless profiling shows `_B12_2` dominates. For ordinary He-Xe operating temperatures, fewer terms would likely be enough.

## Why inversion is the awkward part

The direct backend has no inverse table. Every `T_from_h` or `T_from_s` call solves a nonlinear inverse problem with Newton iterations.

The current implementation recomputes expensive shared work inside each Newton step:

- `T_from_h` calls `enthalpy(fp,T,P)` and then `cp(fp,T,P)`.
- `T_from_s` calls `entropy(fp,T,P)` and then `cp(fp,T,P)`.
- Each of those scalar property calls recomputes `_virial(fp,T)` and the molar-density root.

That means one Newton step pays for two full property evaluations even though both evaluations need the same virial coefficients and density at the same `(T,P)`.

This cost is multiplied in cycle solves because compressor/turbine outlet calculations, heat exchangers, and radiators all use inverse properties. It is multiplied again in the design solver's ForwardDiff Jacobian, where the same inversion logic runs with Dual numbers.

`T_from_s` also accepts `T_guess` but currently ignores it and always uses an ideal-gas entropy estimate. That estimate is usually a good start, but ignoring the caller's guess is surprising and should be revisited when the inverse helper is refactored.

## Transport reminder

DONE 2026-06-10: `viscosity`, `conductivity`, `prandtl` added to the interface and implemented for `NobleGasMixture` (pair constants resolved in the constructor; 0 allocations; μ 458 ns / k 495 ns / Pr 1.3 μs at He-Xe 900 K). Remaining micro-opportunity: `prandtl` evaluates the density root three times (once each in `cp`, `μ`, `k`); a bundled state API would share it.
