# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

GasCycle.jl — a Julia package for closed Brayton cycle analysis (He-Xe working fluid, space nuclear power), built as a clean reimplementation of the NPSS cycle-analysis workflow. Two pillars: full ForwardDiff automatic differentiation through the entire cycle (including Newton solves), and validation against NPSS 3.3 output to print precision (see `validation/RESULTS.md`).

## Commands

```bash
# Full test suite (~530 tests)
julia --project -e 'using Pkg; Pkg.test()'

# Single test file (each is standalone: does `using Test; using GasCycle` itself)
julia --project -e 'include("test/test_thermo.jl")'

# Benchmarks (separate project env in benchmarks/)
julia benchmarks/runbenchmarks.jl

# Validation scripts (NPSS cross-validation replicas)
julia --project validation/bru3_hexe_run.jl
```

CI (`.github/workflows/CI.yml`) deletes `Manifest.toml` and resolves fresh to check `[compat]` bounds — the committed Manifest is pinned to the dev machine's Julia version.

## Architecture

Include-order layering in `src/GasCycle.jl` (single module, no submodules): units → thermo → core → maps → elements → network → solver → output. Lower layers never reference higher ones.

- **`thermo/`** — `FluidProperties` abstract interface (`cp`, `enthalpy`, `entropy`, `density`, `gamma`, `viscosity`, `conductivity`, `prandtl` + inversions `T_from_h`, `T_from_s`, `h_from_s`). Elements call only this interface, never a concrete backend. Five backends: `NobleGasMixture`/`HeXe(M)` (analytic virial EOS + transport, Tournier/El-Genk AIAA 2006-4154 — the default for He-Xe work), `FPTFluid` (reads NPSS FPT table files — kept for apples-to-apples NPSS comparison), `IdealGasFluid`, `PolynomialLiquid` (coolants with polynomial T-dependence; reads `data/Water.fpt`, `data/DC200.fpt`, `data/WaterEG50.fpt`), `ConstantPropertyLiquid` (constant coolants; legacy NPSS-replication files).
- **`core/`** — `FluidState{T<:Real}` (Pt, Tt, W, fluid), `Port` (mutable ref to a state), `AbstractElement`.
- **`elements/`** — Compressor, Turbine, HeatExchanger, HeatSource, Radiator, Duct, Splitter, Mixer, Shaft. Each has design and (where applicable) off-design/map modes plus the off-design residual interface (`n_residuals`, `residuals`, `indep_vars`, `set_indep_vars!`).
- **`network/FlowNetwork.jl`** — directed port graph; computes a flow plan; loops (recuperator, cold-end closure) are declared as back-edges.
- **`solver/`** — `Solver.jl` includes `DesignSolve.jl`, `OffDesignSolve.jl`, `ResidualAssembly.jl`, etc. Back-edge `[Tt, Pt]` states are Newton unknowns; `one_pass!` propagates the cycle from seeds; NonlinearSolve.jl with `AutoForwardDiff()` closes the residual. Off-design adds map coordinates (`Wc_map` per turbomachine) and shaft speed to the same Newton vector (TrustRegion).

### AD is a hard constraint

Everything is parametric in the numeric type (`Compressor{T<:Real}`, `FluidState{T}`) so ForwardDiff Duals flow through the whole cycle, including nested through the inner Newton solves (implicit function theorem via differently-tagged Duals). When writing or modifying code: no `Float64` type annotations on values that derive from states or parameters, no non-differentiable branches on Dual values. Gradient correctness is tested against finite differences.

## Units

Internal computation is **SI everywhere** (K, Pa, J, kg). NPSS models, FPT files, and NASA reports use English units (°R, psia, BTU/lbm, lbm/s) — conversion helpers live in `src/units.jl`; use them rather than redefining factors. FPT files are always Rankine/psia/BTU and are converted to SI on load. `validation/` deliberately reports in °R/psia for direct comparison with NPSS prints.

## Gotchas

- `cp` is **not exported** (clashes with `Base.Filesystem.cp` on Julia ≥ 1.12) — use `import GasCycle: cp` or `GasCycle.cp(...)`.
- `FPTFluid`'s default `s_interp=:log_pressure` interpolates pressure-detrended entropy; `:linear` reproduces NPSS's linear-in-P scheme (a known artifact, documented in `validation/RESULTS.md`). Don't "fix" the `:linear` path — it exists for NPSS replication.
- NPSS `effDes` is **isentropic** efficiency; `η_type=:isentropic`/`:polytropic` selects semantics on Compressor/Turbine.
- `NobleGasMixture` scalar calls must stay allocation-free; performance history and the current optimization levers are in `src/thermo/NobleGasMixture_perf_notes.md`. Inversion call sites should pass a `T_guess`.

## Reference material

- `reference/` — the NPSS model being replicated (`BRU3.mdl`, `HeXe.out`) and the property/design papers (El-Genk AIAA 2006-4154, Johnson NASA/CR-2006-214394, NASA TN D-5815).
- `data/` — the FPT property files (`HeXe84.fpt`, `Oil.fpt`, `H2O.fpt`) with provenance and audit notes in `data/README.md`. Oil.fpt deliberately does not match real Dow Corning 200, and HeXe84.fpt's Cp/gam/Pr tables have a known wrong-signed real-gas departure — don't "fix" these files; they replicate the NPSS inputs byte-for-byte.
- `validation/PLAN.md` + `validation/RESULTS.md` — the NPSS cross-validation ladder and its findings ledger. New discrepancies against NPSS should be logged there, not just fixed silently.
- `ROADMAP.md` — ordered next steps (inventory control, transients, additional working fluids).
- `docs/GUIDE.md` + `docs/REFERENCE.md` — user-facing guide and API reference. The reference lists every exported function and element option — keep it in sync when changing public API.
