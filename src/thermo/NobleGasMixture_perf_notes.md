# NobleGasMixture Performance Notes

Review date: 2026-06-10 (baseline re-measured 2026-06-11)

## Current baseline (2026-06-11, `SUITE["backends"]` in benchmarks/)

Identical recuperated-design model and gas (He-Xe M = 83.8) across the three
backends.  NOTE: the FPT numbers improved an order of magnitude vs the
2026-06-10 measurements (commit 35d1023 made inversions table-direct), so the
old "direct inversions beat FPT" conclusion is obsolete — NobleGas is now the
slow backend everywhere:

| benchmark (min time)        | idealgas | fpt     | noblegas pre | post-opt | vs fpt |
|-----------------------------|----------|---------|--------------|----------|--------|
| scalar forward5 (5 props)   | 25 ns    | 373 ns  | 1.67 μs      | 1.63 μs  | 4.5×   |
| T_from_h                    | 2.3 ns   | 74 ns   | 678 ns       | 345 ns   | 4.4×   |
| T_from_s                    | 24 ns    | 93 ns   | 2.27 μs      | 359 ns   | 4.0×   |
| h_from_s                    | 24 ns    | 97 ns   | 2.60 μs      | 1.40 μs  | 15×    |
| AD dT/dh                    | 2.6 ns   | 84 ns   | 666 ns       | 370 ns   | 4.1×   |
| AD dT/ds                    | 24 ns    | 96 ns   | 3.34 μs      | 381 ns   | 3.7×   |
| solve! recuperated-design   | 171 μs   | 290 μs  | 2.13 ms      | 1.40 ms  | 4.9×   |
| ForwardDiff cycle gradient  | 231 μs   | 386 μs  | 4.12 ms      | 2.79 ms  | 7.4×   |

"post-opt" = after the 2026-06-11 inversion optimization (fused
`_h_cp`/`_s_cp` Newton helpers — one virial pass + one density root per
iteration — and `T_from_s` honoring the caller's `T_guess`; every element
call site passes one).  Inversions now cost ~1 fused evaluation from a
good guess, so the remaining solve/gradient gap tracks the scalar forward
cost: every property call pays its own ~330 ns virial pass + density
root.  The next lever is the bundled state API ("Main reminders") so
element code asking for several properties at one (T,P) pays once —
projected to bring solve! near ~500 μs.  `h_from_s` still uses the
generic no-guess fallback (closed-form start + full convergence); no
element calls it, so it stays unoptimized.

## Main reminders

- ~~Fix `HeXe(M_molar)` for AD.~~ DONE 2026-06-10 (name built from `ForwardDiff.value`).
- ~~Add shared helpers for inverse solves~~ DONE 2026-06-11: `_h_cp`/`_s_cp` fused helpers; scalar properties also refactored onto shared molar kernels (`_h_molar`, `_s_molar`, `_cp_molar`, `_cv_molar`), which de-triplicated the Eq. 13 cp body.
- Consider a bundled state API, such as `thermo_state(fp,T,P)`, returning `rho`, `h`, `s`, `cp`, and `gamma` together. Element code often asks for multiple properties at the same state. **Now the top remaining lever** — inversions are ~1 fused eval; the solve/gradient gap vs FPT is the per-property virial recomputation in forward calls.
- Precompute mixture constants in the `NobleGasMixture` constructor: `x2` and composition weights (transport pair constants are now resolved in the constructor, 2026-06-10).
- Keep the full 18-term LJ virial series for broad validity unless profiling shows `_B12_2` dominates. For ordinary He-Xe operating temperatures, fewer terms would likely be enough.

## Why inversion was the awkward part (RESOLVED 2026-06-11)

The direct backend has no inverse table; every `T_from_h`/`T_from_s` is a Newton solve. The original implementation paid two full property evaluations per Newton step (`enthalpy` then `cp`, each with its own `_virial` pass and density root) and `T_from_s` ignored the caller's `T_guess`, always restarting from the ideal-gas estimate.

Both fixed: Newton steps now use the fused `_h_cp`/`_s_cp` helpers (one virial pass, one density root per step), and `T_from_s` uses the caller's guess when given (`T_guess = nothing` falls back to the ideal-gas closed form, which is nearly exact at CBC conditions — real-gas departure ~1 K). The `T_from_s` step is capped at halving `T` per iteration so a poor caller guess cannot leave the log domain. Breaking on the pre-update residual applies the final dual correction before exit, which is what keeps the AD rows at parity with the primal rows.

## Transport reminder

DONE 2026-06-10: `viscosity`, `conductivity`, `prandtl` added to the interface and implemented for `NobleGasMixture` (pair constants resolved in the constructor; 0 allocations; μ 458 ns / k 495 ns / Pr 1.3 μs at He-Xe 900 K). Remaining micro-opportunity: `prandtl` evaluates the density root three times (once each in `cp`, `μ`, `k`); a bundled state API would share it.
