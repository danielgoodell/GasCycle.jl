using Test
using GasCycle

"""
Off-design map-based operation (NPSS-style formulation).

Unknowns per off-design turbomachine: the map flow coordinate Wc_map.
Unknown per off-design shaft: speed N.
Residuals: Wc_map - Wc_actual per machine, shaft power balance.

Gold-standard check: an off-design solve at design boundary conditions, with
maps scaled through the design point, must recover the design shaft speed,
pressure ratios, and efficiencies.
"""

# Synthetic base maps shared by the off-design testsets: smooth analytic
# surfaces with physically-signed slopes (compressor PR falls with flow,
# turbine PR rises with flow; both rise with speed).
function synthetic_base_maps()
    Nc_ax = collect(0.5:0.05:1.3)
    Wc_ax = collect(0.4:0.05:1.4)
    cbase = PerformanceMap(Nc_ax, Wc_ax,
        [1.0 + 1.5 * n^2 * (1.3 - 0.5 * w) for n in Nc_ax, w in Wc_ax],
        [0.83 - 0.3 * (w - n)^2            for n in Nc_ax, w in Wc_ax])
    tbase = PerformanceMap(Nc_ax, Wc_ax,
        [1.0 + 2.0 * w * sqrt(n)           for n in Nc_ax, w in Wc_ax],
        [0.88 - 0.2 * (w - n)^2            for n in Nc_ax, w in Wc_ax])
    (cbase, tbase)
end

@testset "Off-design map operation" begin
    fluid = IdealGasFluid(M_molar = 83.8)   # HeXe84

    P1, T1, W = 2.0e5, 400.0, 1.0
    PR_c, η_c = 2.0, 0.85
    η_t       = 0.90
    TIT       = 1100.0
    N_des     = 40_000.0

    # ── Design point: fixed-PR chain; bisect turbine PR for shaft balance ────
    comp_d = Compressor("Comp"; PR=PR_c, η_poly=η_c)
    heat_d = HeatSource("Reactor"; TtExit=TIT, dPqP=0.02)

    p2 = compute!(comp_d, Port(FluidState(P1, T1, W, fluid)))
    p3 = compute!(heat_d, p2)
    w_comp = specific_work(comp_d)

    function turb_work(PR)
        t = Turbine("T"; PR=PR, η_poly=η_t)
        compute!(t, Port(p3[]))
        specific_work(t)
    end
    lo, hi = 1.2, 4.0
    for _ in 1:60
        mid = 0.5 * (lo + hi)
        turb_work(mid) < w_comp ? (lo = mid) : (hi = mid)
    end
    PR_t = 0.5 * (lo + hi)
    @test turb_work(PR_t) ≈ w_comp rtol=1e-8

    # Design-point corrected conditions for map scaling
    Tt3, Pt3 = p3[].Tt, p3[].Pt
    Nc_c, Wc_c = corrected_speed(N_des, T1),  corrected_flow(W, T1, P1)
    Nc_t, Wc_t = corrected_speed(N_des, Tt3), corrected_flow(W, Tt3, Pt3)

    cbase, tbase = synthetic_base_maps()

    # Reference point chosen off-node so the design point lands inside smooth
    # interpolation cells rather than on a grid kink.
    cmap = scale_map(cbase; Nc_des=Nc_c, Wc_des=Wc_c, PR_des=PR_c, eta_des=η_c,
                     Nc_ref=0.93, Wc_ref=0.87)
    tmap = scale_map(tbase; Nc_des=Nc_t, Wc_des=Wc_t, PR_des=PR_t, eta_des=η_t,
                     Nc_ref=0.93, Wc_ref=0.87)

    function build_offdesign_net(; TtExit, N0)
        comp  = Compressor("Comp"; η_poly=η_c, map=cmap, mode=:off_design)
        heat  = HeatSource("Reactor"; TtExit=TtExit, dPqP=0.02)
        turb  = Turbine("Turb"; η_poly=η_t, map=tmap, mode=:off_design)
        shaft = Shaft("Main"; N=N0, mode=:off_design)

        net = FlowNetwork()
        add!(net, comp, heat, turb)
        connect!(net, comp => heat => turb)
        add_shaft!(net, shaft; drives=comp, driven_by=turb)
        set_state!(net, comp; Pt=P1, Tt=T1, W=W, fluid=fluid)
        (net, comp, turb, shaft)
    end

    @testset "design-point reproduction" begin
        net, comp, turb, shaft = build_offdesign_net(TtExit=TIT, N0=0.95 * N_des)
        sol = solve!(net)

        @test sol.status == :success
        @test shaft.N      ≈ N_des rtol=1e-3
        @test comp.PR      ≈ PR_c  rtol=1e-3
        @test comp.η_poly  ≈ η_c   rtol=1e-3
        @test comp.Wc_map  ≈ Wc_c  rtol=1e-3
        @test turb.Wc_map  ≈ Wc_t  rtol=1e-3
        @test turb.inlet[].Pt / turb.outlet[].Pt ≈ PR_t rtol=1e-3
        @test abs(power_balance(shaft)) < 1e-3 * w_comp * W

        # Regression for the Wc_map = Wc_act bug: the flow-match residual must
        # be a real constraint, not structurally zero.
        Wc_saved = comp.Wc_map
        comp.Wc_map = 1.05 * Wc_saved
        @test abs(residuals(comp)[1]) > 1e-3
        comp.Wc_map = Wc_saved
    end

    @testset "reduced turbine inlet temperature" begin
        net, comp, turb, shaft = build_offdesign_net(TtExit=1000.0, N0=N_des)
        sol = solve!(net)

        @test sol.status == :success
        @test abs(residuals(comp)[1])  < 1e-6
        @test abs(residuals(turb)[1])  < 1e-6
        @test abs(power_balance(shaft)) < 1e-3 * w_comp * W
        # Less turbine specific work available → shaft settles below design speed
        @test shaft.N < N_des
        @test shaft.N > 0.5 * N_des
        @test comp.PR < PR_c
    end
end

@testset "Off-design closed loop (recuperated BRU-like)" begin
    fluid = IdealGasFluid(M_molar = 83.8)

    # BRU-like parameters (NASA TN D-5815, ideal-gas stand-in for HeXe84)
    T1, P1, W = 300.0, 163.4e3, 0.6
    PR_c, η_c = 1.9, 0.80
    η_t       = 0.87
    ε_recup   = 0.95
    dPqP_cold, dPqP_hot, dPqP_heat = 0.011, 0.022, 0.027
    TIT       = 1144.0
    N_des     = 36_000.0

    # Turbine exhausts so the loop pressure closes back to P1 through the
    # recuperator hot side.
    P_exit = P1 / (1 - dPqP_hot)

    function build_loop(; comp_kw=(;), turb_kw=(;), shaft_kw=(;), TtExit=TIT)
        comp  = Compressor("Comp"; η_poly=η_c, comp_kw...)
        recup = HeatExchanger("Recup"; ε=ε_recup,
                              dPqP_hot=dPqP_hot, dPqP_cold=dPqP_cold)
        heat  = HeatSource("Heater"; TtExit=TtExit, dPqP=dPqP_heat)
        turb  = Turbine("Turb"; η_poly=η_t, turb_kw...)
        shaft = Shaft("Shaft"; shaft_kw...)

        net = FlowNetwork()
        add!(net, comp, recup, heat, turb)
        connect!(net, comp => recup => heat => turb => comp)
        add_shaft!(net, shaft; drives=comp, driven_by=turb)
        add_hx_pair!(net, recup; hot=turb)
        set_state!(net, comp; Pt=P1, Tt=T1, W=W, fluid=fluid)
        (net, comp, recup, heat, turb, shaft)
    end

    # ── Design solve (back-edge Newton; fixed PR, pressure closure) ──────────
    net_d, comp_d, recup_d, heat_d, turb_d, shaft_d =
        build_loop(comp_kw=(PR=PR_c,),
                   turb_kw=(mode=:pressure_closure, P_exit=P_exit),
                   shaft_kw=(N=N_des,))
    sol_d = solve!(net_d)
    @test sol_d.status == :success

    P_net_des = net_power(sol_d)
    PR_t_des  = pressure_ratio(turb_d)
    @test P_net_des > 0

    s_t = turb_d.inlet[]
    Nc_c, Wc_c = corrected_speed(N_des, T1),     corrected_flow(W, T1, P1)
    Nc_t, Wc_t = corrected_speed(N_des, s_t.Tt), corrected_flow(W, s_t.Tt, s_t.Pt)

    cbase, tbase = synthetic_base_maps()
    cmap = scale_map(cbase; Nc_des=Nc_c, Wc_des=Wc_c, PR_des=PR_c, eta_des=η_c,
                     Nc_ref=0.93, Wc_ref=0.87)
    tmap = scale_map(tbase; Nc_des=Nc_t, Wc_des=Wc_t, PR_des=PR_t_des, eta_des=η_t,
                     Nc_ref=0.93, Wc_ref=0.87)

    @testset "design-point reproduction with generator load" begin
        net, comp, recup, heat, turb, shaft =
            build_loop(comp_kw=(map=cmap, mode=:off_design),
                       turb_kw=(map=tmap, mode=:off_design),
                       shaft_kw=(N=0.97 * N_des, mode=:off_design, P_load=P_net_des))
        sol = solve!(net)

        @test sol.status == :success
        @test shaft.N ≈ N_des rtol=1e-3
        @test pressure_ratio(comp) ≈ PR_c     rtol=1e-3
        @test pressure_ratio(turb) ≈ PR_t_des rtol=1e-3
        @test comp.Wc_map ≈ Wc_c rtol=1e-3
        @test turb.Wc_map ≈ Wc_t rtol=1e-3
        @test power_balance(shaft) ≈ P_net_des rtol=1e-3
    end

    @testset "TIT sweep at constant speed (plan phase 7)" begin
        # Alternator-locked shaft: N fixed at design (no shaft residual);
        # the map operating points and the loop back-edge are the unknowns.
        net, comp, recup, heat, turb, shaft =
            build_loop(comp_kw=(map=cmap, mode=:off_design),
                       turb_kw=(map=tmap, mode=:off_design),
                       shaft_kw=(N=N_des,))

        fracs  = collect(1.0:-0.05:0.6)
        powers = Float64[]
        for frac in fracs
            heat.TtExit = frac * TIT
            sol = solve!(net)
            @test sol.status == :success
            @test abs(residuals(comp)[1]) < 1e-6
            @test abs(residuals(turb)[1]) < 1e-6
            push!(powers, net_power(sol))
        end

        # 100% TIT must reproduce the design point; power falls monotonically
        # as TIT is throttled back.
        @test powers[1] ≈ P_net_des rtol=1e-3
        @test all(diff(powers) .< 0)
    end
end
