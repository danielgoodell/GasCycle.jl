"""
Rung-1 of the NPSS cross-validation ladder (see PLAN.md): single elements in
isolation, mirroring the diagnostic section at the end of reference/BRU3.mdl
line for line.  All inputs are the exact NPSS values (English units in, so
the table prints °R / psia for digit-by-digit comparison with NPSS couts).

NPSS output columns are PENDING until the collaborator provides the run
listing (PLAN.md §"NPSS-side artifacts needed").  In-loop anchor values from
BRU3.mdl comments are shown where they exist — they are loop states, not
isolation outputs, so treat them as ~anchors, not oracles.

Also computes the two rung-3 attribution factors available today:
  - efficiency semantics: effDes (isentropic) vs polytropic at same η
  - fluid backend: HeXe84.fpt vs ideal monatomic gas M = 83.8
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using GasCycle
using Printf

const root = joinpath(@__DIR__, "..")
fpt   = FPTFluid(joinpath(root, "HeXe84.fpt"))
ideal = IdealGasFluid(M_molar = 83.8)
oil   = ConstantPropertyLiquid(joinpath(root, "Oil.fpt"))   # BRU3 "Oil" (Dow 200)

toR(T_K)     = K_to_R(T_K)
topsia(P_Pa) = Pa_to_psia(P_Pa)
st(Tt_R, Pt_psia, W_lbps, fl) =
    FluidState(psia_to_Pa(Pt_psia), R_to_K(Tt_R), lbps_to_kgps(W_lbps), fl)

rule() = println("─"^78)

# ── 1. Compressor in isolation ────────────────────────────────────────────────
# Not in the .mdl diagnostic section, but station 1 (737 °R / 45.03 psia) is
# the strongest in-loop anchor we have.  NPSS: PRdes=1.9, effDes=0.80.
println("\n[1] Compressor isolation — inlet 540 °R / 23.7 psia / 1.32 lb/s, PR=1.9")
rule()
@printf("%-34s %12s %12s\n", "variant", "Tt_out [°R]", "Pt_out [psia]")
rule()
comp_out = Dict{String,Float64}()
for (fl_name, fl) in ("HeXe84.fpt" => fpt, "ideal M=83.8" => ideal),
    (sem_name, sem) in ("effDes (isentropic)" => :isentropic,
                        "polytropic"          => :polytropic)
    c = Compressor("C"; PR = 1.9, η_poly = 0.80, η_type = sem)
    out = compute!(c, Port(st(540.0, 23.7, 1.32, fl)))[]
    comp_out["$fl_name/$sem_name"] = toR(out.Tt)
    @printf("%-34s %12.2f %12.3f\n", "$fl_name, $sem_name", toR(out.Tt), topsia(out.Pt))
end
@printf("%-34s %12s %12s\n", "NPSS (CEAT, effDes)  [in-loop]", "≈737", "45.03")
rule()
@printf("Δ semantics (FPT, isen−poly):   %+7.2f °R\n",
        comp_out["HeXe84.fpt/effDes (isentropic)"] - comp_out["HeXe84.fpt/polytropic"])
@printf("Δ fluid    (isen, FPT−ideal):   %+7.2f °R\n",
        comp_out["HeXe84.fpt/effDes (isentropic)"] - comp_out["ideal M=83.8/effDes (isentropic)"])
@printf("Δ vs NPSS anchor (FPT, isen):   %+7.2f °R   ← CEAT-vs-HeXe84 table diff (+ loop effects)\n",
        comp_out["HeXe84.fpt/effDes (isentropic)"] - 737.0)

# ── 2. Turbine in isolation (BRU3.mdl lines ~288-294) ─────────────────────────
# Turb.Fl_I.setTotalTP(2060., 43.2); W = 1.31; PRdes = 43.2/24.69; effDes=0.87
PR_t = 43.2 / 24.69
println("\n[2] Turbine isolation — inlet 2060 °R / 43.2 psia / 1.31 lb/s, PR=$(round(PR_t, digits=5))")
rule()
@printf("%-34s %12s %12s\n", "variant", "Tt_out [°R]", "Pt_out [psia]")
rule()
turb_out = Dict{String,Float64}()
for (fl_name, fl) in ("HeXe84.fpt" => fpt, "ideal M=83.8" => ideal),
    (sem_name, sem) in ("effDes (isentropic)" => :isentropic,
                        "polytropic"          => :polytropic)
    t = Turbine("T"; PR = PR_t, η_poly = 0.87, η_type = sem)
    out = compute!(t, Port(st(2060.0, 43.2, 1.31, fl)))[]
    turb_out["$fl_name/$sem_name"] = toR(out.Tt)
    @printf("%-34s %12.2f %12.3f\n", "$fl_name, $sem_name", toR(out.Tt), topsia(out.Pt))
end
@printf("%-34s %12s %12s\n", "NPSS isolation cout", "PENDING", "PENDING")
rule()
@printf("Δ semantics (FPT, isen−poly):   %+7.2f °R\n",
        turb_out["HeXe84.fpt/effDes (isentropic)"] - turb_out["HeXe84.fpt/polytropic"])
@printf("Δ fluid    (isen, FPT−ideal):   %+7.2f °R\n",
        turb_out["HeXe84.fpt/effDes (isentropic)"] - turb_out["ideal M=83.8/effDes (isentropic)"])
println("(in-loop station 6 = 1701 °R is the HotStart tear value, not a turbine output)")

# ── 3. Recuperator in isolation (BRU3.mdl lines ~296-302) ─────────────────────
# Fl_I1 (cold/high-P) 737 °R / 45.03 psia / 1.267 lb/s
# Fl_I2 (hot /low-P) 1701 °R / 24.68 psia / 1.28 lb/s
# effect = 0.95, dPqP1 = 0.011 (cold), dPqP2 = 0.022 (hot)
println("\n[3] Recuperator isolation — cold 737 °R/45.03 psia/1.267, hot 1701 °R/24.68 psia/1.28, ε=0.95")
rule()
@printf("%-22s %16s %16s\n", "fluid", "cold_out [°R]", "hot_out [°R]")
rule()
for (fl_name, fl) in ("HeXe84.fpt" => fpt, "ideal M=83.8" => ideal)
    hx = HeatExchanger("Recup"; ε = 0.95, dPqP_cold = 0.011, dPqP_hot = 0.022)
    hx.cold_inlet = Port(st(737.0, 45.03, 1.267, fl))
    hx.hot_inlet  = Port(st(1701.0, 24.68, 1.28, fl))
    compute_hx!(hx)
    @printf("%-22s %16.2f %16.2f\n", fl_name,
            toR(hx.cold_outlet[].Tt), toR(hx.hot_outlet[].Tt))
end
@printf("%-22s %16s %16s\n", "NPSS isolation cout", "PENDING", "PENDING")
println("(in-loop station 8 anchor: hot_out ≈ 786 °R)")

# ── 4. HeatSinkHx in isolation (BRU3.mdl lines ~304-311) ──────────────────────
# Fl_I1 (gas) 786 °R / 1.267 lb/s at the prior-run pressure (≈ station-8
#   pressure: 24.68 psia × (1 − 0.022) ≈ 24.14 psia)
# Fl_I2 (oil) 527 °R / 1.28 lb/s (sic — the .mdl reuses 1.28; the design
#   oil flow is 0.14 lb/s) at ≈ pump outlet 90.2 psia
# effect = 0.946, dPqP1 = (14.45−14.2)/14.45 ≈ 0.01730, dPqP2 = 0.005
println("\n[4] HeatSinkHx isolation — gas 786 °R/1.267 lb/s, oil(Dow200) 527 °R/1.28 lb/s, ε=0.946")
rule()
@printf("%-22s %16s %16s\n", "gas fluid", "gas_out [°R]", "oil_out [°R]")
rule()
P_gas = 24.68 * (1 - 0.022)
for (fl_name, fl) in ("HeXe84.fpt" => fpt, "ideal M=83.8" => ideal)
    hx = HeatExchanger("HeatSinkHx"; ε = 0.946,
                       dPqP_hot = (14.45 - 14.2) / 14.45, dPqP_cold = 0.005)
    hx.hot_inlet  = Port(st(786.0, P_gas, 1.267, fl))
    hx.cold_inlet = Port(st(527.0, 90.2, 1.28, oil))
    compute_hx!(hx)
    @printf("%-22s %16.2f %16.2f\n", fl_name,
            toR(hx.hot_outlet[].Tt), toR(hx.cold_outlet[].Tt))
end
@printf("%-22s %16s %16s\n", "NPSS isolation cout", "PENDING", "PENDING")
println("(loop closure check: gas_out should land near 540 °R = comp inlet)")

println("\nDone.  Fill PENDING columns from the NPSS run listing / isolation couts;")
println("rung-0 (property lookup) comparison comes first if any rung-1 row disagrees.")
