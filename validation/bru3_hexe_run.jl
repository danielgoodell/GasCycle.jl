"""
Rung 0/1/2 against reference/HeXe.out — the authoritative apples-to-apples
target: NPSS 3.3 running BRU3.mdl with the SAME HeXe84.fpt table GasCycle
reads (run by dgoodell, 06/10/26, converged in 4 iterations).

Configuration recovered from the listing:
  - machines isentropic (eff = 0.8000 / 0.8700 exact; efPoly derived)
  - 2 % compressor-exit bleed, forced to 559 °R, reinjected at turbine inlet
  - heater dPqP = 0.027 present (43.332/44.535); sink gas dPqP = 0.0173
  - turbine PR floated for pressure closure → 1.757, exit 24.660 psia
  - TIT floated for shaft balance vs HPX → 2024.04 °R
  - HPX displays 18.00 (hp) = 13.42 kW: the .mdl assigns 13.4198 with an
    explicit "kW" tag and NPSS stores hp internally — units question closed
  - oil loop: pump exit 522.88 °R into the sink HX; implied oil cp from the
    listing's energy balance = 0.7999 BTU/(lbm·R) = Oil.fpt's 0.8

GasCycle is run in both entropy-interpolation modes:
  :linear        — NPSS-compat (NPSS's FPT lookup is linear in P)
  :log_pressure  — physically correct (GasCycle default)
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using GasCycle
import GasCycle: cp
using Printf

const root = joinpath(@__DIR__, "..")
oil = ConstantPropertyLiquid(joinpath(root, "data", "Oil.fpt"))

T0, P0 = R_to_K(540.0), psia_to_Pa(23.7)
W      = lbps_to_kgps(1.32)
P_load = 13.4198e3                      # W; .mdl HPX expression, kW-tagged
P_exit = P0 / ((1 - 0.022) * (1 - (14.45 - 14.2) / 14.45))

toR(K) = K_to_R(K)
rule() = println("─"^78)

# ── Rung 0: property lookups vs the listing's ht / s / gamt columns ──────────
println("\n[0] Property lookups at listing states (HeXe84.fpt; NPSS prints 2-3 decimals)")
rule()
@printf("%-26s %-6s %9s %9s %9s %9s %9s\n",
        "station (Tt °R, Pt psia)", "mode", "ht GC", "ht NPSS", "s GC", "s NPSS", "γ GC")
for (lbl, Tt_R, Pt_psia, ht_n, s_n, gam_n) in
    (("st0 (540.00, 23.700)",  540.00, 23.700,  61.54, 0.282, 1.66743),
     ("st1 (751.47, 45.030)",  751.47, 45.030,  74.03, 0.286, 1.66880),
     ("st5 (2024.04, 43.332)", 2024.04, 43.332, 149.31, 0.345, 1.66778),
     ("st7 (1624.41, 24.660)", 1624.41, 24.660, 125.70, 0.346, 1.66746))
    for (mode, fl) in (("linear", FPTFluid(joinpath(root, "data", "HeXe84.fpt"); s_interp = :linear)),
                       ("log_p",  FPTFluid(joinpath(root, "data", "HeXe84.fpt"))))
        T, P = R_to_K(Tt_R), psia_to_Pa(Pt_psia)
        @printf("%-26s %-6s %9.2f %9.2f %9.4f %9.3f %9.5f\n", lbl, mode,
                Jkg_to_btulbm(enthalpy(fl, T, P)), ht_n,
                JkgK_to_btulbmR(entropy(fl, T, P)), s_n, gamma(fl, T, P))
        lbl = ""; ht_n = NaN; s_n = NaN
    end
end
println("(γ NPSS column: 1.66743 / 1.66880 / 1.66778 / 1.66746)")

# ── Rungs 1-2: full loop in both entropy modes ───────────────────────────────
function build_and_solve(fluid; TIT0 = R_to_K(2024.0))
    comp   = Compressor("Comp"; PR = 1.9, η_poly = 0.80, η_type = :isentropic)
    bsplit = Splitter("BldSplit"; fracs = [0.98, 0.02])
    recup  = HeatExchanger("Recup"; ε = 0.95, dPqP_cold = 0.011, dPqP_hot = 0.022)
    bcool  = HeatSource("BldCool"; TtExit = R_to_K(559.0), mode = :fixed_TtExit, dPqP = 0.0)
    heater = HeatSource("Heater"; TtExit = TIT0, dPqP = 0.027)
    bmix   = Mixer("BldMix")
    turb   = Turbine("Turb"; mode = :pressure_closure, P_exit = P_exit,
                     η_poly = 0.87, η_type = :isentropic)
    sink   = HeatExchanger("HeatSinkHx"; ε = 0.946,
                           dPqP_hot = (14.45 - 14.2) / 14.45, dPqP_cold = 0.005)
    shaft  = Shaft("Shaft"; N = 36_000.0)

    net = FlowNetwork()
    add!(net, comp, bsplit, recup, bcool, heater, bmix, turb, sink)
    connect!(net, comp => bsplit => recup => heater => bmix => turb)
    connect_port!(net, bsplit, :bleed_outlet, bcool, :inlet)
    connect_port!(net, bcool, :outlet, bmix, :bleed_inlet)
    add_shaft!(net, shaft; drives = comp, driven_by = turb)
    add_hx_pair!(net, recup; hot = turb)
    connect_port!(net, recup, :hot_outlet, sink, :hot_inlet)
    set_state!(net, comp; Pt = P0, Tt = T0, W = W, fluid = fluid)
    set_boundary!(net, sink, :cold_inlet;
                  Pt = psia_to_Pa(90.653), Tt = R_to_K(522.88),
                  W = lbps_to_kgps(0.14), fluid = oil)

    # TIT from shaft balance (NPSS: ind_HeatSource ↔ HPX power balance)
    lo, hi = R_to_K(1700.0), R_to_K(2400.0)
    local sol
    for _ in 1:60
        heater.TtExit = 0.5 * (lo + hi)
        sol = solve!(net)
        sol.status == :success || error("solve failed at TIT=$(toR(heater.TtExit))")
        net_power(sol) < P_load ? (lo = heater.TtExit) : (hi = heater.TtExit)
    end
    (; sol, comp, recup, heater, bmix, turb, sink, shaft)
end

npss = [  # (label, value °R or psia or kW) from HeXe.out
    ("TIT = st5 Tt [°R]",          2024.04),
    ("st5 Pt [psia]",                43.332),
    ("st1 comp out Tt [°R]",        751.47),
    ("turb inlet (mixed) PR",        1.757),
    ("st6 turb out Tt [°R]",       1624.41),
    ("st6 turb out Pt [psia]",       24.660),
    ("st4 recup cold out Tt [°R]", 1580.77),
    ("st8 recup hot out Tt [°R]",   810.74),
    ("st8 Pt [psia]",                24.117),
    ("st9 sink gas out Tt [°R]",    540.00),
    ("comp power [kW]",               17.1),
    ("turb power [kW]",               30.5),
    ("net = HPX [kW]",               13.42),
]

function gascycle_row(m)
    Dict(
        "TIT = st5 Tt [°R]"          => toR(m.heater.outlet[].Tt),
        "st5 Pt [psia]"              => Pa_to_psia(m.heater.outlet[].Pt),
        "st1 comp out Tt [°R]"       => toR(m.comp.outlet[].Tt),
        "turb inlet (mixed) PR"      => pressure_ratio(m.turb),
        "st6 turb out Tt [°R]"       => toR(m.turb.outlet[].Tt),
        "st6 turb out Pt [psia]"     => Pa_to_psia(m.turb.outlet[].Pt),
        "st4 recup cold out Tt [°R]" => toR(m.recup.cold_outlet[].Tt),
        "st8 recup hot out Tt [°R]"  => toR(m.recup.hot_outlet[].Tt),
        "st8 Pt [psia]"              => Pa_to_psia(m.recup.hot_outlet[].Pt),
        "st9 sink gas out Tt [°R]"   => toR(m.sink.hot_outlet[].Tt),
        "comp power [kW]"            => specific_work(m.comp) * m.comp.inlet[].W / 1e3,
        "turb power [kW]"            => specific_work(m.turb) * m.turb.inlet[].W / 1e3,
        "net = HPX [kW]"             => net_power(m.sol) / 1e3,
    )
end

m_lin = build_and_solve(FPTFluid(joinpath(root, "data", "HeXe84.fpt"); s_interp = :linear))
m_log = build_and_solve(FPTFluid(joinpath(root, "data", "HeXe84.fpt")))
g_lin, g_log = gascycle_row(m_lin), gascycle_row(m_log)

println("\n[1] Full-loop comparison vs HeXe.out (NPSS-compat :linear and physical :log_pressure)")
rule()
@printf("%-28s %10s %10s %8s %10s\n", "quantity", "NPSS", "GC linear", "Δ", "GC log_p")
rule()
for (lbl, v) in npss
    d = g_lin[lbl] - v
    @printf("%-28s %10.3f %10.3f %8.3f %10.3f\n", lbl, v, g_lin[lbl], d, g_log[lbl])
end
rule()
println("""
Notes:
  - NPSS station temps print 2 decimals; powers 1 decimal (±0.05 kW).
  - st9 uses the listing's oil inlet (522.88 °R) as a fixed boundary; in
    NPSS the oil temperature floats to force exactly 540.00, so the GC st9
    value measures HX-formulation agreement, not closure.""")

# ── [2] NPSS-faithful bleed bookkeeping (bleed = workless bypass) ─────────────
# The listing's BLEEDS section (dhb/dh = 0, dPb/dP = 1) plus the power and
# outlet-temperature forensics show NPSS's bleed does no work in either
# machine: comp pwr is main-flow-only Δh (17.05 ≈ "17.1"), and the turbine
# expands the main flow alone, with the 559 °R bleed diluting it at the
# EXIT (main-only outlet 1646.4 °R → mixed 1624.6 ≈ 1624.41).  Replicate by
# splitting before the compressor and mixing after the turbine, with the
# bleed "freely" repressurized to 45.03 psia exactly as NPSS grants it.
println("\n[2] NPSS-faithful topology: bleed bypasses both machines (s_interp=:linear)")
rule()
function build_and_solve_npss_faithful(fluid)
    bsplit = Splitter("BldSplit"; fracs = [0.98, 0.02])
    comp   = Compressor("Comp"; PR = 1.9, η_poly = 0.80, η_type = :isentropic)
    recup  = HeatExchanger("Recup"; ε = 0.95, dPqP_cold = 0.011, dPqP_hot = 0.022)
    # dhb/dh = 0 and the 559 °R bearing forcing; dPb/dP = 1 grants the bleed
    # compressor-exit pressure for free: dPqP = 1 − 45.03/23.7
    bcool  = HeatSource("BldCool"; TtExit = R_to_K(559.0), mode = :fixed_TtExit,
                        dPqP = 1 - 45.03 / 23.7)
    heater = HeatSource("Heater"; TtExit = R_to_K(2024.0), dPqP = 0.027)
    turb   = Turbine("Turb"; mode = :pressure_closure, P_exit = P_exit,
                     η_poly = 0.87, η_type = :isentropic)
    bmix   = Mixer("BldMix")
    sink   = HeatExchanger("HeatSinkHx"; ε = 0.946,
                           dPqP_hot = (14.45 - 14.2) / 14.45, dPqP_cold = 0.005)
    shaft  = Shaft("Shaft"; N = 36_000.0)

    net = FlowNetwork()
    add!(net, bsplit, comp, recup, bcool, heater, turb, bmix, sink)
    connect!(net, bsplit => comp => recup => heater => turb => bmix)
    connect_port!(net, bsplit, :bleed_outlet, bcool, :inlet)
    connect_port!(net, bcool, :outlet, bmix, :bleed_inlet)
    add_shaft!(net, shaft; drives = comp, driven_by = turb)
    add_hx_pair!(net, recup; hot = bmix)
    connect_port!(net, recup, :hot_outlet, sink, :hot_inlet)
    set_state!(net, bsplit; Pt = P0, Tt = T0, W = W, fluid = fluid)
    set_boundary!(net, sink, :cold_inlet;
                  Pt = psia_to_Pa(90.653), Tt = R_to_K(522.88),
                  W = lbps_to_kgps(0.14), fluid = oil)

    lo, hi = R_to_K(1700.0), R_to_K(2400.0)
    local sol
    for _ in 1:60
        heater.TtExit = 0.5 * (lo + hi)
        sol = solve!(net)
        sol.status == :success || error("solve failed at TIT=$(toR(heater.TtExit))")
        net_power(sol) < P_load ? (lo = heater.TtExit) : (hi = heater.TtExit)
    end
    (; sol, comp, recup, heater, turb, bmix, sink, shaft)
end

f = build_and_solve_npss_faithful(FPTFluid(joinpath(root, "data", "HeXe84.fpt"); s_interp = :linear))
@printf("%-44s %10s %10s %8s\n", "quantity", "GasCycle", "NPSS", "Δ")
rule()
for (lbl, gc, np) in (
    ("TIT = st5 Tt [°R]  (both solved)", toR(f.heater.outlet[].Tt), 2024.04),
    ("st1 comp out Tt [°R]",             toR(f.comp.outlet[].Tt),    751.47),
    ("st6 turb out, after mix [°R]",     toR(f.bmix.outlet[].Tt),   1624.41),
    ("st6 Pt [psia]",                    Pa_to_psia(f.bmix.outlet[].Pt), 24.660),
    ("turbine PR",                       pressure_ratio(f.turb),      1.757),
    ("st4 recup cold out Tt [°R]",       toR(f.recup.cold_outlet[].Tt), 1580.77),
    ("st8 recup hot out Tt [°R]",        toR(f.recup.hot_outlet[].Tt),  810.74),
    ("comp pwr [kW]",  specific_work(f.comp) * f.comp.inlet[].W / 1e3,   17.1),
    ("turb pwr [kW]",  specific_work(f.turb) * f.turb.inlet[].W / 1e3,   30.5),
    ("net = HPX [kW]", net_power(f.sol) / 1e3,                          13.42))
    @printf("%-44s %10.3f %10.3f %8.3f\n", lbl, gc, np, gc - np)
end

# Effective sink-HX ε implied by the listing (one remaining open question)
ε_sink_npss = (810.74 - 540.0) / (810.74 - 522.88)
@printf("\nsink HX: NPSS behaves as ε = %.4f (effect set 0.946; recup matches its 0.95 exactly)\n",
        ε_sink_npss)
