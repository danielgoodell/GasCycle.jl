"""
Replication of the actual NPSS output recorded in
reference/BRU-ModelOutput-ParameterBookKeeping.xlsx ("NPSS" column), and
decomposition of its flagged Test-vs-NPSS temperature discrepancies.

That spreadsheet is the first true NPSS output listing we have, and it
recontextualizes the .mdl comments: the 737 / 1701 / 786 °R "anchors" are
the **Test Data** column (NASA TN D-5815 design/test values), not NPSS
outputs.  The recorded NPSS run differs from the .mdl as-saved:

  - TIT pinned at 2060 °R (station 5 exact), no shaft-balance involvement
  - heater pressure drop = 0  (st4→st5: 44.535 → 44.535; .mdl says 0.027)
    ⇒ turbine PR = 44.535/24.356 = 1.8285, not ~1.75
  - bleed flowing: 0.03 lb/s (2.27 %) at 540 °R (not heated to 559 °R)
  - HeatSinkHx gas-side dPqP = 0.005 (not the .mdl's 0.0173)
  - the HotStart tear destroys 0.04 lb/s (st6 W=1.32, st7 W=1.28)
    and pins the recup hot inlet at the converged 1642.05 °R

Because TIT and the tear are pinned, every element is explicit — no loop
solve.  Each element below is fed the exact NPSS-column inlet state, so
disagreement isolates element/property formulation, never propagation.
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using GasCycle
import GasCycle: cp
using Printf

# CEAT-equivalent fluid (see RESULTS.md: CEA noble-gas ≡ ideal monatomic)
M_He, M_Xe = 4.002602, 131.293
w_He  = 0.0181
fluid = IdealGasFluid(M_molar = 1 / (w_He / M_He + (1 - w_He) / M_Xe))

toR(K)  = K_to_R(K)
st(Tt_R, Pt_psia, W_lbps) =
    FluidState(psia_to_Pa(Pt_psia), R_to_K(Tt_R), lbps_to_kgps(W_lbps), fluid)

rule() = println("─"^74)

# ── Station 1: compressor (540/23.7 → PR 1.9), NPSS column: 741.7 °R ─────────
println("\n[1] Compressor: 540 °R / 23.7 psia, PR = 1.9   (NPSS output: 741.70 °R)")
rule()
for (lbl, sem, η) in (("effDes 0.80 isentropic", :isentropic, 0.80),
                      ("0.80 polytropic",        :polytropic, 0.80),
                      ("0.784 isentropic",       :isentropic, 0.784),
                      ("0.809 polytropic",       :polytropic, 0.809))
    c = Compressor("C"; PR = 1.9, η_poly = η, η_type = sem)
    out = compute!(c, Port(st(540.0, 23.7, 1.32)))[]
    @printf("  %-26s -> %8.2f °R   (Δ vs NPSS %+6.2f)\n", lbl, toR(out.Tt),
            toR(out.Tt) - 741.70)
end

# ── Station 6: turbine with bleed, NPSS column: 1642.07 °R ───────────────────
# Inlet: 1.29 lb/s at 2060 °R mixed with 0.03 lb/s bleed at 540 °R (Pfract=0
# ⇒ bleed expands through the full turbine).  PR = 44.535/24.356.
PR_t = 44.535 / 24.356
W_mix  = 1.32
Tt_mix_R = (1.29 * 2060.0 + 0.03 * 540.0) / 1.32      # 2025.45 °R (ideal gas: cp cancels)
println("\n[2] Turbine: mixed inlet $(round(Tt_mix_R, digits=2)) °R / 44.535 psia, PR = $(round(PR_t, digits=4))   (NPSS output: 1642.07 °R)")
rule()
for (lbl, sem, η) in (("effDes 0.87 isentropic", :isentropic, 0.87),
                      ("0.87 polytropic",        :polytropic, 0.87))
    t = Turbine("T"; PR = PR_t, η_poly = η, η_type = sem)
    out = compute!(t, Port(st(Tt_mix_R, 44.535, W_mix)))[]
    @printf("  %-26s -> %8.2f °R   (Δ vs NPSS %+6.2f)\n", lbl, toR(out.Tt),
            toR(out.Tt) - 1642.07)
end

# ── Stations 4 & 8: recuperator at the exact NPSS port states ────────────────
# cold: 741.7 °R / 45.03 psia / 1.29   hot: 1642.05 °R / 24.355 psia / 1.28
println("\n[3] Recuperator: cold 741.7/45.03/1.29, hot 1642.05/24.355/1.28, ε = 0.95")
rule()
hx = HeatExchanger("Recup"; ε = 0.95, dPqP_cold = 0.011, dPqP_hot = 0.022)
hx.cold_inlet = Port(st(741.70, 45.03, 1.29))
hx.hot_inlet  = Port(st(1642.05, 24.355, 1.28))
compute_hx!(hx)
@printf("  cold out  %8.2f °R   (NPSS 1580.02, Δ %+6.2f)\n",
        toR(hx.cold_outlet[].Tt), toR(hx.cold_outlet[].Tt) - 1580.02)
@printf("  hot  out  %8.2f °R   (NPSS  794.77, Δ %+6.2f)\n",
        toR(hx.hot_outlet[].Tt), toR(hx.hot_outlet[].Tt) - 794.77)
ε_impl_hot  = (1642.05 - 794.77) / (1642.05 - 741.70)
ε_impl_cold = (1580.02 - 741.70) / (1642.05 - 741.70) * (1.29 / 1.28)
@printf("  implied NPSS ε (C_min = hot basis):  hot side %.4f, cold side %.4f\n",
        ε_impl_hot, ε_impl_cold)

# ── Stations 9 & 12: HeatSinkHx, gas out and oil out ─────────────────────────
# gas: 794.77 °R / 1.28 lb/s.  oil in: 525.54 °R / 0.14 lb/s.
println("\n[4] HeatSinkHx: gas 794.77/1.28, oil 525.54/0.14, ε = 0.946")
rule()
cp_gas_btu = JkgK_to_btulbmR(cp(fluid, 400.0, 1e5))
Q_gas      = 1.28 * cp_gas_btu * (794.77 - 540.0)            # BTU/s (NPSS gas out = 540)
T9_pred    = 794.77 - 0.946 * (794.77 - 525.54)
@printf("  ε-NTU gas out (gas C_min):  %7.2f °R  (NPSS 540.00, Δ %+5.2f)\n",
        T9_pred, T9_pred - 540.0)
cp_oil_implied = Q_gas / (0.14 * (731.56 - 525.54))
@printf("  oil cp implied by NPSS oil out 731.56 °R:    %6.4f BTU/(lbm·R) = %5.3f kJ/(kg·K)\n",
        cp_oil_implied, btulbmR_to_JkgK(cp_oil_implied) / 1e3)
cp_oil_test = (1.32 * cp_gas_btu * (786.0 - 540.0)) / (0.14 * (746.0 - 527.0))
@printf("  oil cp implied by TEST column (746 °R out):  %6.4f BTU/(lbm·R) = %5.3f kJ/(kg·K)\n",
        cp_oil_test, btulbmR_to_JkgK(cp_oil_test) / 1e3)
println("  candidates on the sheet: 'Current Value' 1.8, Dow200 'correct' 0.884,")
println("  Oil.fpt 0.8  [BTU/(lbm·R)] — none matches the implied values above.")

# ── Decomposition of the flagged Test-vs-NPSS diffs ──────────────────────────
println("\n[5] Decomposition of the spreadsheet's flagged discrepancies")
rule()
T_nobleed_175  = 2060.0 * (1 - 0.87 * (1 - (1 / 1.75)^0.4))
T_nobleed_1828 = 2060.0 * (1 - 0.87 * (1 - (1 / PR_t)^0.4))
T_mix_1828     = Tt_mix_R * (1 - 0.87 * (1 - (1 / PR_t)^0.4))
@printf("""  Station 6/7 (turbine out), test 1701 vs NPSS 1642.07 (Δ ≈ 59 °R):
      missing heater ΔP ⇒ PR 1.8285 not 1.75:  %+6.1f °R
      2.27%% bleed dilution at 540 °R:           %+6.1f °R
      residual (semantics/CEA/η):                %+6.1f °R
""",
        T_nobleed_1828 - T_nobleed_175,
        T_mix_1828 - T_nobleed_1828,
        1642.07 - T_mix_1828)
println("  Station 4 (recup cold out), test 1653 vs NPSS 1580 (Δ ≈ 73 °R):")
println("      ≈ ε × (station-7 Δ) = 0.95 × 59 ≈ 56 °R, plus the comp-outlet Δ and")
println("      the unequal recup flows (1.29/1.28 vs test 1.267/1.28) — i.e. fully")
println("      downstream of the station-6 causes, not an independent recup issue.")
println("  Station 12 (oil out), test 746 vs NPSS 731.6 (Δ ≈ 14 °R):")
println("      oil-side cp bookkeeping (see [4]) — the sheet itself flags oil Cp.")
