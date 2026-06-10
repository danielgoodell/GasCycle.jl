# NPSS Cross-Validation Results

Running log. See PLAN.md for the ladder definition. English units (°R, psia)
for direct comparison with NPSS prints.

## 2026-06-10 — Rung 0: entropy interpolation artifact found and fixed

Probing FPT property lookups at the BRU station states (`HeXe84.fpt`):
cp, γ, ρ all within ~1 % of ideal monatomic M = 83.8, but the s-table's
pressure derivative was wrong:

    ds/dlnP = −113.23 J/(kg·K)   measured at 300 K, 165 kPa (mid-cell)
            = −99.22             required (−R; real-gas correction ≪ 1 % here)

Root cause (confirmed to 4 significant figures): bilinear interpolation
**linear in P** of s ∝ −R·ln P across the coarse Pt grid (nodes at 101.6
and 198.2 kPa bracket the BRU compressor inlet at 163.4 kPa — worst case).
Predicted chord slope −R·P·ln(P₂/P₁)/(P₂−P₁) = −113.23 — identical to the
measured value. This single artifact inflated the implied isentropic
exponent to γ_eff ≈ 1.74 (> 5/3, unphysical) and produced **+14 °R at the
compressor outlet** — most of the historical "fluid model difference".

**Fix:** `FPTFluid(...; s_interp=:log_pressure)` (now the default)
interpolates the detrended σ = s + R·ln(Pt/P_ref) and restores the log term
analytically; R is fitted from the table itself (median −Δs/ΔlnP over all
node pairs: 99.27 vs 99.22 from M = 83.8). Exact at table nodes; exact in P
for ideal gas. `s_interp=:linear` keeps the legacy behavior for
apples-to-apples runs against interpolators that are linear in P (NPSS's
scheme: to be confirmed at rung 0 with the run listing).

After the fix the isentropic step at the compressor inlet implies
γ_eff = 1.666 vs table γ = 1.667. Regression tests in `test_thermo.jl`.

## 2026-06-10 — Rung 1: isolation results (GasCycle side; NPSS couts pending)

`julia validation/bru3_isolation.jl`, HeXe84.fpt, effDes (isentropic):

| Element (BRU3.mdl isolation inputs)        | GasCycle | NPSS      | anchor |
|--------------------------------------------|----------|-----------|--------|
| Compressor out (540 °R/23.7 psia, PR 1.9)  | 737.16 °R| PENDING   | ≈737 °R (in-loop) |
| Turbine out (2060 °R/43.2 psia, PR 1.7497) | 1700.83 °R| PENDING  | 1701 °R (tear value) |
| Recup cold out (737 + 1701 °R, ε 0.95)     | 1652.74 °R| PENDING  | ~1651 °R (.mdl note) |
| Recup hot out                              | 793.51 °R| PENDING   | 786 °R (in-loop) |
| HeatSinkHx gas out (786 °R + 527 °R oil)   | 540.92 °R| PENDING   | ≈540 °R (loop closure) |

Sub-°R agreement with every isolation-consistent anchor. Factor sizes at
the compressor outlet (°R): s-interp artifact +14.3 (fixed), semantics
isen−poly −6.3, HeXe84-vs-ideal −0.4.

Key inference: the BRU3.mdl 1701 °R HotStart tear value **is** the no-bleed
isolated turbine output (we get 1700.83). The .mdl pins station 6 there
even though its own 2 % bleed (cooled to 559 °R, reinjected at the turbine
inlet) lowers the effective turbine inlet to ~2030 °R and hence the outlet
to ~1673 °R — exactly the .mdl header's open issue ("look at Interstage
Bleed Port").

## 2026-06-10 — Rung 2 preview: full loop (examples/bru_10kw.jl, isentropic)

| Quantity        | GasCycle | NPSS/anchor | status |
|-----------------|----------|-------------|--------|
| Comp outlet     | 737.0 °R | ≈737 °R     | matches at print precision |
| Turb outlet     | 1673 °R  | 1701 °R     | Δ = bleed premix vs no-bleed tear (see above) |
| Recup cold out  | 1627 °R  | ~1651 °R    | cascades from turb outlet Δ |
| Heater input    | 34.95 kW | ~33.1 kW    | cascades from recup cold out Δ |
| Net shaft power | 13.12 kW | "13.4" (HPX) | HPX hp-vs-kW ambiguity pending |
| Est. electrical | 10.6 kW  | 10.5 kW design | — |

Everything left on the list traces to (a) the bleed/tear treatment — needs
the NPSS run listing to see what NPSS actually converges to, and ideally a
no-bleed run of both models, and (b) the CEAT-vs-HeXe84 table at hot-end
conditions — bounded small by the turbine isolation row but only directly
checkable with the artifacts in PLAN.md.
