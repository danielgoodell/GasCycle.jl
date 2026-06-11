using Pkg

Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=dirname(@__DIR__)))
Pkg.instantiate()

using BenchmarkTools
using ForwardDiff
using GasCycle

const SUITE = BenchmarkGroup()

function fpt_forward5(fluid, T, P)
    GasCycle.cp(fluid, T, P),
    enthalpy(fluid, T, P),
    entropy(fluid, T, P),
    density(fluid, T, P),
    gamma(fluid, T, P)
end

fpt_dTdh(fluid, h, P, T_guess) =
    ForwardDiff.derivative(hx -> T_from_h(fluid, hx, P; T_guess), h)

fpt_dTds(fluid, s, P, T_guess) =
    ForwardDiff.derivative(sx -> T_from_s(fluid, sx, P; T_guess), s)

function simple_brayton_net()
    fluid = HeXeIdealGas(0.47)
    T_in, P_in, W_flow = 400.0, 500e3, 10.0

    comp = Compressor("Comp"; PR=2.5, η_poly=0.87)
    heat = HeatSource("Reactor"; TtExit=1100.0, dPqP=0.02)
    turb = Turbine("Turb"; mode=:pressure_closure, P_exit=P_in, η_poly=0.90)
    shaft = Shaft("Main"; N=15_000.0)

    net = FlowNetwork()
    add!(net, comp, heat, turb)
    connect!(net, comp => heat => turb => comp)
    add_shaft!(net, shaft; drives=comp, driven_by=turb)
    set_state!(net, comp; Pt=P_in, Tt=T_in, W=W_flow, fluid=fluid)
    solve!(net)
    net
end

function recuperated_design_net()
    fluid = HeXeIdealGas(0.47)
    T_in, P_in, W_flow = 400.0, 500e3, 10.0

    comp = Compressor("Comp"; PR=2.5, η_poly=0.87)
    recup = HeatExchanger("Recup"; ε=0.92)
    heat = HeatSource("Reactor"; TtExit=1100.0, dPqP=0.02)
    turb = Turbine("Turb"; mode=:pressure_closure, P_exit=P_in, η_poly=0.90)
    shaft = Shaft("Main"; N=15_000.0)

    net = FlowNetwork()
    add!(net, comp, recup, heat, turb)
    connect!(net, comp => recup => heat => turb => comp)
    add_shaft!(net, shaft; drives=comp, driven_by=turb)
    add_hx_pair!(net, recup; hot=turb)
    set_state!(net, comp; Pt=P_in, Tt=T_in, W=W_flow, fluid=fluid)
    solve!(net)
    net
end

function fpt_recuperated_design_net(fluid)
    T_in, P_in, W_flow = 400.0, 500e3, 10.0

    comp = Compressor("Comp"; PR=2.5, η_poly=0.87)
    recup = HeatExchanger("Recup"; ε=0.92)
    heat = HeatSource("Reactor"; TtExit=1100.0, dPqP=0.02)
    turb = Turbine("Turb"; mode=:pressure_closure, P_exit=P_in, η_poly=0.90)
    shaft = Shaft("Main"; N=15_000.0)

    net = FlowNetwork()
    add!(net, comp, recup, heat, turb)
    connect!(net, comp => recup => heat => turb => comp)
    add_shaft!(net, shaft; drives=comp, driven_by=turb)
    add_hx_pair!(net, recup; hot=turb)
    set_state!(net, comp; Pt=P_in, Tt=T_in, W=W_flow, fluid=fluid)
    solve!(net)
    net
end

function synthetic_base_maps()
    Nc_ax = collect(0.5:0.05:1.3)
    Wc_ax = collect(0.4:0.05:1.4)
    cbase = PerformanceMap(Nc_ax, Wc_ax,
        [1.0 + 1.5 * n^2 * (1.3 - 0.5 * w) for n in Nc_ax, w in Wc_ax],
        [0.83 - 0.3 * (w - n)^2 for n in Nc_ax, w in Wc_ax])
    tbase = PerformanceMap(Nc_ax, Wc_ax,
        [1.0 + 2.0 * w * sqrt(n) for n in Nc_ax, w in Wc_ax],
        [0.88 - 0.2 * (w - n)^2 for n in Nc_ax, w in Wc_ax])
    (cbase, tbase)
end

function offdesign_map_net()
    fluid = IdealGasFluid(M_molar=83.8)

    P1, T1, W = 2.0e5, 400.0, 1.0
    PR_c, η_c = 2.0, 0.85
    η_t, TIT, N_des = 0.90, 1100.0, 40_000.0

    comp_d = Compressor("Comp"; PR=PR_c, η_poly=η_c)
    heat_d = HeatSource("Reactor"; TtExit=TIT, dPqP=0.02)
    p2 = compute!(comp_d, Port(FluidState(P1, T1, W, fluid)))
    p3 = compute!(heat_d, p2)
    w_comp = specific_work(comp_d)

    function turb_work(PR)
        turb = Turbine("Turb"; PR=PR, η_poly=η_t)
        compute!(turb, Port(p3[]))
        specific_work(turb)
    end

    lo, hi = 1.2, 4.0
    for _ in 1:60
        mid = 0.5 * (lo + hi)
        turb_work(mid) < w_comp ? (lo = mid) : (hi = mid)
    end
    PR_t = 0.5 * (lo + hi)

    Tt3, Pt3 = p3[].Tt, p3[].Pt
    Nc_c, Wc_c = corrected_speed(N_des, T1), corrected_flow(W, T1, P1)
    Nc_t, Wc_t = corrected_speed(N_des, Tt3), corrected_flow(W, Tt3, Pt3)

    cbase, tbase = synthetic_base_maps()
    cmap = scale_map(cbase; Nc_des=Nc_c, Wc_des=Wc_c, PR_des=PR_c,
                     eta_des=η_c, Nc_ref=0.93, Wc_ref=0.87)
    tmap = scale_map(tbase; Nc_des=Nc_t, Wc_des=Wc_t, PR_des=PR_t,
                     eta_des=η_t, Nc_ref=0.93, Wc_ref=0.87)

    comp = Compressor("Comp"; η_poly=η_c, map=cmap, mode=:off_design)
    heat = HeatSource("Reactor"; TtExit=TIT, dPqP=0.02)
    turb = Turbine("Turb"; η_poly=η_t, map=tmap, mode=:off_design)
    shaft = Shaft("Main"; N=0.95 * N_des, mode=:off_design)

    net = FlowNetwork()
    add!(net, comp, heat, turb)
    connect!(net, comp => heat => turb)
    add_shaft!(net, shaft; drives=comp, driven_by=turb)
    set_state!(net, comp; Pt=P1, Tt=T1, W=W, fluid=fluid)
    solve!(net)
    net
end

function cycle_power(x::AbstractVector)
    PR_comp, ε_recup = x[1], x[2]
    fluid = HeXeIdealGas(0.47)
    T0, P0, W = 400.0, 500e3, 10.0

    comp = Compressor("Comp"; PR=PR_comp, η_poly=0.88)
    recup = HeatExchanger("Recup"; ε=ε_recup, dPqP_hot=0.01, dPqP_cold=0.01)
    heat = HeatSource("Heater"; TtExit=1100.0, dPqP=0.02)
    p_exit = P0 / ((1 - 0.01) * (1 - 0.01))
    turb = Turbine("Turb"; mode=:pressure_closure, P_exit=p_exit, η_poly=0.90)

    net = FlowNetwork()
    add!(net, comp, recup, heat, turb)
    connect!(net, comp => recup => heat => turb => comp)
    add_hx_pair!(net, recup; hot=turb)
    set_state!(net, comp; Pt=P0, Tt=T0, W=W, fluid=fluid)
    sol = solve!(net; maxiter=200)
    net_power(sol)
end

function fpt_cycle_power(x::AbstractVector, fluid)
    PR_comp, T_reactor_exit = x[1], x[2]
    T0, P0, W = 400.0, 500e3, 10.0

    comp = Compressor("Comp"; PR=PR_comp, η_poly=0.88)
    recup = HeatExchanger("Recup"; ε=0.90, dPqP_hot=0.01, dPqP_cold=0.01)
    heat = HeatSource("Heater"; TtExit=T_reactor_exit, dPqP=0.02)
    p_exit = P0 / ((1 - 0.01) * (1 - 0.01))
    turb = Turbine("Turb"; mode=:pressure_closure, P_exit=p_exit, η_poly=0.90)

    net = FlowNetwork()
    add!(net, comp, recup, heat, turb)
    connect!(net, comp => recup => heat => turb => comp)
    add_hx_pair!(net, recup; hot=turb)
    set_state!(net, comp; Pt=P0, Tt=T0, W=W, fluid=fluid)
    sol = solve!(net; maxiter=200)
    net_power(sol)
end

function build_suite!()
    simple = simple_brayton_net()
    recup = recuperated_design_net()
    offdesign = offdesign_map_net()
    x0 = [2.5, 0.90]
    fpt = FPTFluid(joinpath(dirname(@__DIR__), "HeXe84.fpt"))
    fpt_recup = fpt_recuperated_design_net(fpt)
    T_fpt, P_fpt = 900.0, 2.0e5
    h_fpt = enthalpy(fpt, T_fpt, P_fpt)
    s_fpt = entropy(fpt, T_fpt, P_fpt)
    x0_fpt = [2.5, 1100.0]

    SUITE["solve"]["simple-brayton"] = @benchmarkable solve!($simple)
    SUITE["solve"]["recuperated-design"] = @benchmarkable solve!($recup)
    SUITE["solve"]["offdesign-map"] = @benchmarkable solve!($offdesign)
    SUITE["sensitivity"]["forwarddiff-gradient"] =
        @benchmarkable ForwardDiff.gradient($cycle_power, $x0)
    SUITE["fpt"]["scalar"]["forward5"] =
        @benchmarkable fpt_forward5($fpt, $T_fpt, $P_fpt)
    SUITE["fpt"]["scalar"]["T-from-h"] =
        @benchmarkable T_from_h($fpt, $h_fpt, $P_fpt; T_guess=$T_fpt)
    SUITE["fpt"]["scalar"]["T-from-s"] =
        @benchmarkable T_from_s($fpt, $s_fpt, $P_fpt; T_guess=$T_fpt)
    SUITE["fpt"]["scalar"]["h-from-s"] =
        @benchmarkable h_from_s($fpt, $s_fpt, $P_fpt)
    SUITE["fpt"]["scalar-ad"]["dTdh"] =
        @benchmarkable fpt_dTdh($fpt, $h_fpt, $P_fpt, $T_fpt)
    SUITE["fpt"]["scalar-ad"]["dTds"] =
        @benchmarkable fpt_dTds($fpt, $s_fpt, $P_fpt, $T_fpt)
    SUITE["fpt"]["solve"]["recuperated-design"] =
        @benchmarkable solve!($fpt_recup)
    SUITE["fpt"]["sensitivity"]["forwarddiff-gradient"] =
        @benchmarkable ForwardDiff.gradient(x -> fpt_cycle_power(x, $fpt), $x0_fpt)

    hexe = HeXe(83.8)
    SUITE["noblegas"]["transport"]["viscosity"] =
        @benchmarkable viscosity($hexe, $T_fpt, $P_fpt)
    SUITE["noblegas"]["transport"]["conductivity"] =
        @benchmarkable conductivity($hexe, $T_fpt, $P_fpt)
    SUITE["noblegas"]["transport"]["prandtl"] =
        @benchmarkable prandtl($hexe, $T_fpt, $P_fpt)

    SUITE
end

function main()
    suite = build_suite!()
    tune!(suite)
    results = run(suite; verbose=true)
    display(results)
end

main()
