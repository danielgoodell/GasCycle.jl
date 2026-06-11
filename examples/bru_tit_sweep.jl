"""
BRU Turbine Inlet Temperature Sweep
Varies TIT from 1460 °R to 2060 °R (design point) at fixed compressor
pressure ratio, mass flow, and all other cycle parameters.
Prints CSV-style rows and returns the same data as named tuples.
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using GasCycle

function run_bru_tit_sweep(; TIT_range_R = range(1460.0, 2060.0, length = 25),
                           print_rows::Bool = true)
    # ── Fluid ─────────────────────────────────────────────────────────────────
    fpt_path = joinpath(@__DIR__, "..", "HeXe84.fpt")
    fluid    = FPTFluid(fpt_path)

    # Unit helpers (R_to_K, K_to_R, psia_to_Pa, lbps_to_kgps, …) are exported by GasCycle.

    # ── Fixed design parameters ───────────────────────────────────────────────
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

    P_turb_exit = P_comp_in / ((1 - dPqP_recup_hot) * (1 - dPqP_downstream))

    # Alternator + parasitic loss estimates (from BRU3.mdl HPX breakdown)
    η_alt           = 0.92
    W_parasitic_kW  = 0.65 + 0.236 + 0.272 + 0.064 + 0.1 + 0.05 + 0.2  # 1.572 kW

    # ── Build network (once) ──────────────────────────────────────────────────
    net    = FlowNetwork()
    comp   = Compressor("Comp";  PR=PR_comp,  η_poly=η_comp)
    recup  = HeatExchanger("Recup"; ε=ε_recup,
                                     dPqP_hot=dPqP_recup_hot,
                                     dPqP_cold=dPqP_recup_cold)
    heater = HeatSource("Heater"; TtExit=R_to_K(2060.0), dPqP=dPqP_heatsrc)
    turb   = Turbine("Turb"; mode=:pressure_closure, P_exit=P_turb_exit, η_poly=η_turb)
    shaft  = Shaft("Shaft"; N=36000.0)

    add!(net, comp, recup, heater, turb)
    connect!(net, comp => recup => heater => turb => comp)
    add_shaft!(net, shaft; drives=comp, driven_by=turb)
    add_hx_pair!(net, recup; hot=turb, cold=comp)
    set_state!(net, comp; Pt=P_comp_in, Tt=T_comp_in, W=W_flow, fluid=fluid)

    print_rows && println("TIT_R, TIT_K, W_shaft_kW, W_elec_kW, Q_heater_kW, Q_recup_kW, eta_cycle_pct, T_turb_out_R, T_recup_hot_out_R, PR_turb")

    results = NamedTuple[]
    for TIT_R in TIT_range_R
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

        T_tout  = K_to_R(turb.outlet[].Tt)
        T_rout  = K_to_R(recup.hot_outlet[].Tt)
        PR_t    = pressure_ratio(turb)

        row = (TIT_R = TIT_R,
               TIT_K = R_to_K(TIT_R),
               W_shaft_kW = W_shaft,
               W_elec_kW = W_elec,
               Q_heater_kW = Q_heat,
               Q_recup_kW = Q_rec,
               eta_cycle_pct = η_cyc,
               T_turb_out_R = T_tout,
               T_recup_hot_out_R = T_rout,
               PR_turb = PR_t)

        if print_rows
            println("$(round(row.TIT_R,digits=1)), $(round(row.TIT_K,digits=1)), " *
                    "$(round(row.W_shaft_kW,digits=3)), $(round(row.W_elec_kW,digits=3)), " *
                    "$(round(row.Q_heater_kW,digits=3)), $(round(row.Q_recup_kW,digits=3)), " *
                    "$(round(row.eta_cycle_pct,digits=2)), $(round(row.T_turb_out_R,digits=1)), " *
                    "$(round(row.T_recup_hot_out_R,digits=1)), $(round(row.PR_turb,digits=4))")
        end
        push!(results, row)
    end
    results
end

if abspath(PROGRAM_FILE) == @__FILE__
    results = run_bru_tit_sweep()
    if !isempty(results)
        dp = results[end]
        println("\nDesign point (2060 °R):  W_shaft = $(round(dp.W_shaft_kW,digits=2)) kW,  W_elec = $(round(dp.W_elec_kW,digits=2)) kW,  η = $(round(dp.eta_cycle_pct,digits=1)) %")
        println("Design goal: 10.5 kW net electrical")
    end
end
