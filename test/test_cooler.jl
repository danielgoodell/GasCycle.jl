using Test
using GasCycle
import GasCycle: cp

"""
Cold-end physics (roadmap item 6): ConstantPropertyLiquid coolants,
UA-mode (ε-NTU) heat exchanger, Radiator element, and loop closure through
the cold end so the compressor inlet state responds to the operating point.
"""

const _FPT_DIR = joinpath(@__DIR__, "..", "data")

@testset "ConstantPropertyLiquid backend" begin
    water = ConstantPropertyLiquid(joinpath(_FPT_DIR, "H2O.fpt"))
    oil   = ConstantPropertyLiquid(joinpath(_FPT_DIR, "Oil.fpt"))

    # H2O.fpt: Cp = 1 BTU/(lbm·R), rho = 62.37 lbm/ft³
    @test water.cp_val ≈ btulbmR_to_JkgK(1.0)
    @test water.rho    ≈ lbmft3_to_kgm3(62.37)
    # Oil.fpt: Cp delegates to Cpt(T) = 0.8
    @test oil.cp_val ≈ btulbmR_to_JkgK(0.8)
    @test oil.rho    ≈ lbmft3_to_kgm3(62.424)

    direct = ConstantPropertyLiquid(cp = 4186.8, rho = 999.0, name = "water")
    @test cp(direct, 350.0, 1e5) == 4186.8
    @test density(direct, 350.0, 1e5) == 999.0
    @test gamma(direct, 350.0, 1e5) == 1.0
    @test enthalpy(direct, 350.0, 1e5) ≈ 4186.8 * 350.0
    # Closed-form inversions round-trip exactly
    @test T_from_h(direct, enthalpy(direct, 312.0, 1e5), 1e5) ≈ 312.0
    @test T_from_s(direct, entropy(direct, 312.0, 1e5), 1e5) ≈ 312.0

    @test_throws Exception ConstantPropertyLiquid(cp = -1.0, rho = 999.0)
end

@testset "HeatExchanger UA mode (ε-NTU)" begin
    gas   = IdealGasFluid(M_molar = 83.8)
    water = ConstantPropertyLiquid(joinpath(_FPT_DIR, "H2O.fpt"))

    Th, Tc = 440.0, 290.0
    Wg, Ww = 0.6, 1.0
    C_gas = Wg * cp(gas, Th, 165e3)
    C_w   = Ww * cp(water, Tc, 600e3)
    UA    = 400.0

    net  = FlowNetwork()
    duct = Duct("Inlet"; dPqP = 0.0)
    cool = HeatExchanger("Cooler"; UA = UA, dPqP_hot = 0.01, dPqP_cold = 0.005)
    @test cool.mode == :UA
    add!(net, duct, cool)
    connect_port!(net, duct, :outlet, cool, :hot_inlet)
    set_state!(net, duct; Pt = 165e3, Tt = Th, W = Wg, fluid = gas)
    set_boundary!(net, cool, :cold_inlet; Pt = 600e3, Tt = Tc, W = Ww, fluid = water)
    sol = solve!(net)
    @test sol.status == :success

    # Counter-flow ε-NTU prediction (gas is C_min)
    Cr  = C_gas / C_w
    NTU = UA / C_gas
    e   = exp(-NTU * (1 - Cr))
    ε   = (1 - e) / (1 - Cr * e)
    @test cool.ε ≈ ε rtol = 1e-12
    @test cool.hot_outlet[].Tt ≈ Th - ε * (Th - Tc) rtol = 1e-10

    # Energy balance across both streams
    Q_gas = C_gas * (Th - cool.hot_outlet[].Tt)
    Q_w   = C_w   * (cool.cold_outlet[].Tt - Tc)
    @test Q_gas ≈ Q_w rtol = 1e-10
    @test Q_transferred(cool) ≈ Q_gas rtol = 1e-10

    # More coolant flow → higher effectiveness → colder gas outlet
    Tt_out_1 = cool.hot_outlet[].Tt
    set_boundary!(net, cool, :cold_inlet; Pt = 600e3, Tt = Tc, W = 2Ww, fluid = water)
    solve!(net)
    @test cool.hot_outlet[].Tt < Tt_out_1

    # Balanced-stream limit ε = NTU/(1+NTU)
    @test GasCycle._effectiveness_NTU_counterflow(2.0, 1.0) ≈ 2 / 3
end

@testset "Radiator element" begin
    gas  = HeXeIdealGas(0.72)
    s_in = FluidState(165e3, 437.0, 0.59, gas)

    # Size for a target exit temperature, then reproduce it at fixed area
    rad = Radiator("Rad"; mode = :fixed_TtExit, TtExit = 320.0, T_sink = 230.0,
                   emissivity = 0.85, dPqP = 0.01)
    out = compute!(rad, Port(s_in))
    @test out[].Tt == 320.0
    @test out[].Pt ≈ s_in.Pt * 0.99
    @test rad.A > 0
    @test Q_rejected(rad) ≈ s_in.W * cp(gas, 400.0, 1e5) * (437.0 - 320.0) rtol = 1e-10

    rad_od = Radiator("RadOD"; A = rad.A, T_sink = 230.0,
                      emissivity = 0.85, dPqP = 0.01)
    @test compute!(rad_od, Port(s_in))[].Tt ≈ 320.0 atol = 0.05

    # Hotter inlet at fixed area rejects more heat but exits hotter
    Q_des = Q_rejected(rad_od)
    out_hot = compute!(rad_od, Port(update(s_in; Tt = 500.0)))
    @test out_hot[].Tt > 320.0
    @test Q_rejected(rad_od) > Q_des

    # Segment-count convergence (Heun marching is second order)
    rad_fine = Radiator("RadFine"; A = rad.A, T_sink = 230.0,
                        emissivity = 0.85, dPqP = 0.01, N_seg = 400)
    @test compute!(rad_fine, Port(s_in))[].Tt ≈ compute!(rad_od, Port(s_in))[].Tt atol = 0.01

    # Entering below the sink temperature warms the stream toward the sink
    out_cold = compute!(rad_od, Port(update(s_in; Tt = 210.0)))
    @test 210.0 < out_cold[].Tt < 230.0

    # Sizing-mode guards
    bad = Radiator("Bad"; mode = :fixed_TtExit, TtExit = 220.0, T_sink = 230.0)
    @test_throws Exception compute!(bad, Port(s_in))
    @test_throws Exception Radiator("BadMode"; mode = :nonsense)
    @test_throws Exception Radiator("NoArea")   # :fixed_area needs A > 0
end

# ── Closed loop with floating compressor inlet ────────────────────────────────
#
# BRU-like recuperated loop where the cold end is modeled instead of pinning
# the compressor inlet: recuperator hot-out feeds a cooler (or radiator)
# whose outlet closes the loop through a back-edge.  set_state! supplies
# only mass flow, fluid, and the initial guess.

@testset "Closed loop through the cold end" begin
    fluid = IdealGasFluid(M_molar = 83.8)
    water = ConstantPropertyLiquid(joinpath(_FPT_DIR, "H2O.fpt"))

    T1g, P1, W = 300.0, 163.4e3, 0.6      # initial guess for the comp inlet
    PR_c, η_c  = 1.9, 0.80
    η_t        = 0.87
    ε_recup    = 0.95
    dPqP_cold, dPqP_hot, dPqP_heat, dPqP_cool = 0.011, 0.022, 0.027, 0.01
    TIT        = 1144.0
    N_des      = 36_000.0

    # Turbine exhaust closes the pressure loop: recup hot side then cooler
    P_exit = P1 / ((1 - dPqP_hot) * (1 - dPqP_cool))

    function build_loop(cold_end; comp_kw = (;), turb_kw = (;), shaft_kw = (;))
        comp  = Compressor("Comp"; η_poly = η_c, comp_kw...)
        recup = HeatExchanger("Recup"; ε = ε_recup,
                              dPqP_hot = dPqP_hot, dPqP_cold = dPqP_cold)
        heat  = HeatSource("Heater"; TtExit = TIT, dPqP = dPqP_heat)
        turb  = Turbine("Turb"; η_poly = η_t, turb_kw...)
        shaft = Shaft("Shaft"; shaft_kw...)

        net = FlowNetwork()
        add!(net, comp, recup, heat, turb, cold_end)
        connect!(net, comp => recup => heat => turb)
        add_shaft!(net, shaft; drives = comp, driven_by = turb)
        add_hx_pair!(net, recup; hot = turb)

        if cold_end isa HeatExchanger
            connect_port!(net, recup, :hot_outlet, cold_end, :hot_inlet)
            connect_port!(net, cold_end, :hot_outlet, comp, :inlet; back_edge = true)
            set_boundary!(net, cold_end, :cold_inlet;
                          Pt = 600e3, Tt = 290.0, W = 1.0, fluid = water)
        else  # Radiator
            connect_port!(net, recup, :hot_outlet, cold_end, :inlet)
            connect_port!(net, cold_end, :outlet, comp, :inlet; back_edge = true)
        end
        set_state!(net, comp; Pt = P1, Tt = T1g, W = W, fluid = fluid)
        (net, comp, recup, heat, turb, shaft)
    end

    @testset "design solve with water cooler (UA mode)" begin
        cool = HeatExchanger("Cooler"; UA = 400.0,
                             dPqP_hot = dPqP_cool, dPqP_cold = 0.005)
        net, comp, recup, heat, turb, shaft =
            build_loop(cool; comp_kw = (PR = PR_c,),
                       turb_kw = (mode = :pressure_closure, P_exit = P_exit),
                       shaft_kw = (N = N_des,))
        sol = solve!(net)
        @test sol.status == :success

        # Loop is closed: compressor inlet equals the cooler gas outlet
        s1 = comp.inlet[]
        @test s1.Tt ≈ cool.hot_outlet[].Tt rtol = 1e-8
        @test s1.Pt ≈ cool.hot_outlet[].Pt rtol = 1e-8
        @test s1.Pt ≈ P1 rtol = 1e-8           # pressure closure by construction
        @test s1.Tt != T1g                     # not the seed guess — solved

        # Cooler energy balance gas ↔ water
        C_w   = 1.0 * cp(water, 290.0, 600e3)
        Q_gas = Q_transferred(cool)
        @test C_w * (cool.cold_outlet[].Tt - 290.0) ≈ Q_gas rtol = 1e-8
        @test 290.0 < s1.Tt < 440.0

        # Cycle energy closes: reactor heat = net work + cooler rejection
        Q_in  = heat.Q
        P_net = net_power(sol)
        @test Q_in ≈ P_net + Q_gas rtol = 1e-6

        # Station/component report renders with the cooler in the loop
        rep = sprint(summary, sol)
        @test occursin("Cooler", rep)
        @test occursin("UA=", rep)

        # ── Off-design: maps scaled at this design point ─────────────────────
        s_t  = turb.inlet[]
        PR_t = pressure_ratio(turb)
        T1, P1d = s1.Tt, s1.Pt
        P_des = net_power(sol)

        Nc_ax = collect(0.5:0.05:1.3)
        Wc_ax = collect(0.4:0.05:1.4)
        cbase = PerformanceMap(Nc_ax, Wc_ax,
            [1.0 + 1.5 * n^2 * (1.3 - 0.5 * w) for n in Nc_ax, w in Wc_ax],
            [0.83 - 0.3 * (w - n)^2            for n in Nc_ax, w in Wc_ax])
        tbase = PerformanceMap(Nc_ax, Wc_ax,
            [1.0 + 2.0 * w * sqrt(n)           for n in Nc_ax, w in Wc_ax],
            [0.88 - 0.2 * (w - n)^2            for n in Nc_ax, w in Wc_ax])
        cmap = scale_map(cbase; Nc_des = corrected_speed(N_des, T1),
                         Wc_des = corrected_flow(W, T1, P1d),
                         PR_des = PR_c, eta_des = η_c, Nc_ref = 0.93, Wc_ref = 0.87)
        tmap = scale_map(tbase; Nc_des = corrected_speed(N_des, s_t.Tt),
                         Wc_des = corrected_flow(W, s_t.Tt, s_t.Pt),
                         PR_des = PR_t, eta_des = η_t, Nc_ref = 0.93, Wc_ref = 0.87)

        cool_od = HeatExchanger("Cooler"; UA = 400.0,
                                dPqP_hot = dPqP_cool, dPqP_cold = 0.005)
        net_od, comp_od, recup_od, heat_od, turb_od, shaft_od =
            build_loop(cool_od; comp_kw = (map = cmap, mode = :off_design),
                       turb_kw = (map = tmap, mode = :off_design),
                       shaft_kw = (N = N_des,))   # alternator-locked speed

        # Design-point reproduction with the cold end in the loop
        sol_od = solve!(net_od)
        @test sol_od.status == :success
        @test comp_od.inlet[].Tt ≈ T1  rtol = 1e-3
        @test comp_od.inlet[].Pt ≈ P1d rtol = 1e-3
        @test pressure_ratio(comp_od) ≈ PR_c rtol = 1e-3
        @test net_power(sol_od) ≈ P_des rtol = 1e-3

        # Throttle TIT: less heat into the loop → cooler gas reaches the
        # compressor — the cold end responds instead of staying pinned.
        # The shift is small (the high-ε water cooler anchors the gas outlet
        # near the water inlet temperature) but must be clearly resolved.
        heat_od.TtExit = 0.9 * TIT
        sol_th = solve!(net_od)
        @test sol_th.status == :success
        @test T1 - comp_od.inlet[].Tt > 0.2
        @test comp_od.inlet[].Tt ≈ cool_od.hot_outlet[].Tt rtol = 1e-6
        @test net_power(sol_th) < P_des
    end

    @testset "design sizing + off-design response with radiator" begin
        # Design: size the radiator for a 300 K compressor inlet, 200 K sink
        rad = Radiator("Radiator"; mode = :fixed_TtExit, TtExit = 300.0,
                       T_sink = 200.0, emissivity = 0.85, dPqP = dPqP_cool)
        net, comp, recup, heat, turb, shaft =
            build_loop(rad; comp_kw = (PR = PR_c,),
                       turb_kw = (mode = :pressure_closure, P_exit = P_exit),
                       shaft_kw = (N = N_des,))
        sol = solve!(net)
        @test sol.status == :success
        @test comp.inlet[].Tt ≈ 300.0 rtol = 1e-8
        A_des = rad.A
        @test A_des > 0
        @test Q_rejected(rad) > 0
        @test occursin("Radiator", sprint(summary, sol))

        # Off-design: freeze the area; the same loop must reproduce the
        # design point, then respond when the sink warms up.
        rad.mode = :fixed_area
        sol2 = solve!(net)
        @test sol2.status == :success
        @test comp.inlet[].Tt ≈ 300.0 atol = 0.1

        rad.T_sink = 250.0
        sol3 = solve!(net)
        @test sol3.status == :success
        @test comp.inlet[].Tt > 301.0   # warmer sink → warmer compressor inlet
    end
end
