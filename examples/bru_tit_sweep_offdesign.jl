"""
BRU Turbine Inlet Temperature Sweep — OFF-DESIGN (map-based) version.

Unlike bru_tit_sweep.jl (which re-solves the design problem at each TIT with a
fixed compressor PR), this sweep runs the cycle in true off-design mode:
compressor and turbine operating points come from performance maps scaled
through the design point, and the loop back-edge is solved simultaneously
with the map unknowns.

The shaft is alternator-locked at the design 36 000 rpm (the BRU alternator
is synchronous), so net shaft power varies with TIT.

Modeling assumptions:
  - Compressor inlet held at fixed Pt/Tt (perfect heat sink + inventory
    control).  Loop pressure closure is therefore not re-enforced off-design;
    the turbine exhaust pressure follows from its map PR.
  - Synthetic map shapes (smooth analytic surfaces) scaled through the BRU
    design point — real BRU maps can be dropped in via PerformanceMap.
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using GasCycle

# ── Fluid ─────────────────────────────────────────────────────────────────────
fpt_path = joinpath(@__DIR__, "..", "HeXe84.fpt")
fluid    = FPTFluid(fpt_path)

# ── Unit helpers ──────────────────────────────────────────────────────────────
R_to_K(T_R)     = T_R * (5/9)
psia_to_Pa(P)   = P * 6894.757
lbps_to_kgps(W) = W * 0.453592
K_to_R(T_K)     = T_K * (9/5)

# ── Fixed design parameters (same as bru_tit_sweep.jl) ────────────────────────
T_comp_in = R_to_K(540.0)
P_comp_in = psia_to_Pa(23.7)
W_flow    = lbps_to_kgps(1.32)

PR_comp          = 1.9
η_comp           = 0.80
ε_recup          = 0.95
dPqP_recup_cold  = 0.011
dPqP_recup_hot   = 0.022
dPqP_heatsrc     = 0.027
dPqP_downstream  = 0.017
η_turb           = 0.87
TIT_des_R        = 2060.0
N_des            = 36_000.0

P_turb_exit = P_comp_in / ((1 - dPqP_recup_hot) * (1 - dPqP_downstream))

η_alt           = 0.92
W_parasitic_kW  = 0.65 + 0.236 + 0.272 + 0.064 + 0.1 + 0.05 + 0.2  # 1.572 kW

function build_loop(; comp_kw=(;), turb_kw=(;), shaft_kw=(;))
    comp   = Compressor("Comp"; η_poly=η_comp, comp_kw...)
    recup  = HeatExchanger("Recup"; ε=ε_recup,
                           dPqP_hot=dPqP_recup_hot, dPqP_cold=dPqP_recup_cold)
    heater = HeatSource("Heater"; TtExit=R_to_K(TIT_des_R), dPqP=dPqP_heatsrc)
    turb   = Turbine("Turb"; η_poly=η_turb, turb_kw...)
    shaft  = Shaft("Shaft"; shaft_kw...)

    net = FlowNetwork()
    add!(net, comp, recup, heater, turb)
    connect!(net, comp => recup => heater => turb => comp)
    add_shaft!(net, shaft; drives=comp, driven_by=turb)
    add_hx_pair!(net, recup; hot=turb)
    set_state!(net, comp; Pt=P_comp_in, Tt=T_comp_in, W=W_flow, fluid=fluid)
    (net, comp, recup, heater, turb, shaft)
end

# ── Stage 1: design solve (fixed PR, pressure closure) ────────────────────────
net_d, comp_d, recup_d, heater_d, turb_d, shaft_d =
    build_loop(comp_kw=(PR=PR_comp,),
               turb_kw=(mode=:pressure_closure, P_exit=P_turb_exit),
               shaft_kw=(N=N_des,))
sol_d = solve!(net_d; verbose=true)
sol_d.status == :success || error("design solve failed")

P_net_des = net_power(sol_d)
PR_t_des  = pressure_ratio(turb_d)
s_t       = turb_d.inlet[]
Nc_c, Wc_c = corrected_speed(N_des, T_comp_in), corrected_flow(W_flow, T_comp_in, P_comp_in)
Nc_t, Wc_t = corrected_speed(N_des, s_t.Tt),    corrected_flow(W_flow, s_t.Tt, s_t.Pt)

println("\nDesign point: W_shaft = $(round(P_net_des/1000, digits=2)) kW, " *
        "PR_turb = $(round(PR_t_des, digits=4))")

# ── Stage 2: scale synthetic maps through the design point ────────────────────
Nc_ax = collect(0.5:0.05:1.3)
Wc_ax = collect(0.4:0.05:1.4)
cbase = PerformanceMap(Nc_ax, Wc_ax,
    [1.0 + 1.5 * n^2 * (1.3 - 0.5 * w) for n in Nc_ax, w in Wc_ax],
    [0.83 - 0.3 * (w - n)^2            for n in Nc_ax, w in Wc_ax])
tbase = PerformanceMap(Nc_ax, Wc_ax,
    [1.0 + 2.0 * w * sqrt(n)           for n in Nc_ax, w in Wc_ax],
    [0.88 - 0.2 * (w - n)^2            for n in Nc_ax, w in Wc_ax])

cmap = scale_map(cbase; Nc_des=Nc_c, Wc_des=Wc_c, PR_des=PR_comp, eta_des=η_comp,
                 Nc_ref=0.93, Wc_ref=0.87)
tmap = scale_map(tbase; Nc_des=Nc_t, Wc_des=Wc_t, PR_des=PR_t_des, eta_des=η_turb,
                 Nc_ref=0.93, Wc_ref=0.87)

# ── Stage 3: off-design TIT sweep at constant shaft speed ─────────────────────
net, comp, recup, heater, turb, shaft =
    build_loop(comp_kw=(map=cmap, mode=:off_design),
               turb_kw=(map=tmap, mode=:off_design),
               shaft_kw=(N=N_des,))   # alternator-locked: N fixed, no shaft residual

TIT_range_R = range(0.6 * TIT_des_R, TIT_des_R, length=21)

println("\nTIT_R, TIT_K, W_shaft_kW, W_elec_kW, Q_heater_kW, Q_recup_kW, " *
        "eta_cycle_pct, PR_comp, PR_turb, Wc_turb")

results = []
for TIT_R in reverse(collect(TIT_range_R))   # sweep down from design
    heater.TtExit = R_to_K(TIT_R)
    sol = solve!(net)
    if sol.status != :success
        @warn "Did not converge at TIT = $TIT_R °R"
        continue
    end

    W_shaft = net_power(sol) / 1000
    W_elec  = max(0.0, (W_shaft - W_parasitic_kW) * η_alt)
    Q_heat  = heater.Q / 1000
    Q_rec   = Q_transferred(recup) / 1000
    η_cyc   = cycle_efficiency(sol) * 100

    println("$(round(TIT_R,digits=1)), $(round(R_to_K(TIT_R),digits=1)), " *
            "$(round(W_shaft,digits=3)), $(round(W_elec,digits=3)), " *
            "$(round(Q_heat,digits=3)), $(round(Q_rec,digits=3)), " *
            "$(round(η_cyc,digits=2)), $(round(pressure_ratio(comp),digits=4)), " *
            "$(round(pressure_ratio(turb),digits=4)), $(round(turb.Wc_map,digits=4))")
    push!(results, (TIT_R, W_shaft, W_elec))
end

W_dev = results[1]
println("\nAt design TIT (2060 °R): W_shaft = $(round(W_dev[2],digits=2)) kW " *
        "(design solve gave $(round(P_net_des/1000,digits=2)) kW)")
println("Lowest TIT ($(round(results[end][1],digits=0)) °R): " *
        "W_shaft = $(round(results[end][2],digits=2)) kW")
