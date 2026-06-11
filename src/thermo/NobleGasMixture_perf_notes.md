# NobleGasMixture Performance Notes

Review date: 2026-06-10 (baseline re-measured 2026-06-11)

## Current baseline (2026-06-11, `SUITE["backends"]` in benchmarks/)

Identical recuperated-design model and gas (He-Xe M = 83.8) across the three
backends.  NOTE: the FPT numbers improved an order of magnitude vs the
2026-06-10 measurements (commit 35d1023 made inversions table-direct), so the
old "direct inversions beat FPT" conclusion is obsolete — NobleGas is now the
slow backend everywhere:

| benchmark (min time)        | idealgas | fpt     | noblegas | vs fpt |
|-----------------------------|----------|---------|----------|--------|
| scalar forward5 (5 props)   | 25 ns    | 373 ns  | 1.67 μs  | 4.5×   |
| T_from_h                    | 2.3 ns   | 74 ns   | 678 ns   | 9×     |
| T_from_s                    | 24 ns    | 93 ns   | 2.27 μs  | 24×    |
| h_from_s                    | 24 ns    | 97 ns   | 2.60 μs  | 27×    |
| AD dT/dh                    | 2.6 ns   | 84 ns   | 666 ns   | 8×     |
| AD dT/ds                    | 24 ns    | 96 ns   | 3.34 μs  | 35×    |
| solve! recuperated-design   | 171 μs   | 290 μs  | 2.13 ms  | 7.3×   |
| ForwardDiff cycle gradient  | 231 μs   | 386 μs  | 4.12 ms  | 10.7×  |

All scalar paths remain allocation-free; the gap is pure recomputation —
see "Why inversion is the awkward part" below.  The fix list in "Main
reminders" (shared `_h_cp`/`_s_cp` Newton helpers, honor `T_guess` in
`T_from_s`) directly targets the 24-35× rows.

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
