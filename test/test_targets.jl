using Test
using GasCycle
using ForwardDiff

isdefined(@__MODULE__, :synthetic_comp_map) ||
    include(joinpath(@__DIR__, "synthetic_maps.jl"))

# Recuperated design-point cycle (same as the ForwardDiff test in test_elements).
function recup_net(; η_c = 0.88, W = 10.0, TIT = 1100.0, PR = 2.5)
    fluid  = HeXeIdealGas(0.47)
    T0, P0 = 400.0, 500e3
    net    = FlowNetwork()
    comp   = Compressor("C"; PR = PR, η_poly = η_c)
    recup  = HeatExchanger("R"; ε = 0.90, dPqP_hot = 0.01, dPqP_cold = 0.01)
    heater = HeatSource("H"; TtExit = TIT, dPqP = 0.02)
    turb   = Turbine("T"; mode = :pressure_closure, P_exit = P0 / (0.99 * 0.99), η_poly = 0.90)
    add!(net, comp, recup, heater, turb)
    connect!(net, comp => recup => heater => turb => comp)
    add_hx_pair!(net, recup; hot = turb)
    set_state!(net, comp; Pt = P0, Tt = T0, W = W, fluid = fluid)
    (net, comp, heater, turb)
end

@testset "Targets and free variables — design" begin
    net, = recup_net()
    sol0 = solve!(net)
    P0, W0 = net_power(sol0), 10.0

    @testset "seed mass flow to hit net power" begin
        net, comp, = recup_net()
        add_free!(net, :W)
        add_target!(net, net_power, 1.5 * P0)
        sol = solve!(net)
        @test sol.status == :success
        @test net_power(sol) ≈ 1.5 * P0 rtol = 1e-6
        # Ideal gas: power is linear in flow at fixed temperatures/pressures.
        @test comp.inlet[].W ≈ 1.5 * W0 rtol = 1e-6
    end

    @testset "two pairs: (W, PR) → (power, compressor exit T)" begin
        net, comp, = recup_net()
        add_free!(net, :W)
        add_free!(net, comp, :PR)
        add_target!(net, net_power, 2.0 * P0; name = "power")
        add_target!(net, n -> comp.outlet[].Tt, 520.0; name = "T2")
        sol = solve!(net)
        @test sol.status == :success
        @test net_power(sol) ≈ 2.0 * P0 rtol = 1e-6
        @test comp.outlet[].Tt ≈ 520.0 rtol = 1e-6
        @test comp.PR != 2.5
    end

    @testset "accessor-closure free variable" begin
        net, comp, heater = recup_net()
        add_free!(net, "TIT"; get = () -> heater.TtExit, set! = x -> (heater.TtExit = x))
        add_target!(net, net_power, 1.2 * P0)
        sol = solve!(net)
        @test sol.status == :success
        @test net_power(sol) ≈ 1.2 * P0 rtol = 1e-6
        @test heater.TtExit > 1100.0
    end

    @testset "bounds: unreachable target reports failure at the bound" begin
        net, comp, = recup_net()
        add_free!(net, :W; upper = 12.0)
        add_target!(net, net_power, 1.5 * P0)
        sol = @test_logs (:warn,) match_mode = :any solve!(net; maxiter = 50)
        @test sol.status == :failed
        @test comp.inlet[].W ≈ 12.0
    end

    @testset "errors" begin
        net, comp, = recup_net()
        add_target!(net, net_power, P0)
        @test_throws ErrorException solve!(net)                     # not square
        clear_targets!(net)
        @test solve!(net).status == :success
        @test_throws ErrorException add_free!(net, comp, :nonexistent)
        @test_throws ErrorException add_free!(net, :fluid)
        @test_throws ErrorException add_free!(net, :W; lower = 2.0, upper = 1.0)
        other = Compressor("X")
        @test_throws ErrorException add_free!(net, other, :PR)       # not in network
    end
end

@testset "Targets and free variables — AD through the solve" begin
    # W(η) such that net power = P_req; compare dW/dη with central differences.
    function W_for_power(η; P_req = 3.0e6)
        net, comp, = recup_net(η_c = η)
        add_free!(net, :W)
        add_target!(net, net_power, P_req)
        solve!(net; tol = 1e-10)
        comp.inlet[].W
    end
    η0, h = 0.88, 1e-5
    d_ad = ForwardDiff.derivative(W_for_power, η0)
    d_fd = (W_for_power(η0 + h) - W_for_power(η0 - h)) / 2h
    @test d_ad < 0                      # better compressor → less flow for same power
    @test d_ad ≈ d_fd rtol = 1e-6

    # Second derivative (nested Duals through the solve).
    d2_ad = ForwardDiff.derivative(η -> ForwardDiff.derivative(W_for_power, η), η0)
    d2_fd = (ForwardDiff.derivative(W_for_power, η0 + 1e-4) -
             ForwardDiff.derivative(W_for_power, η0 - 1e-4)) / 2e-4
    @test d2_ad ≈ d2_fd rtol = 1e-5

    # Derivative w.r.t. the target value itself: W ∝ P_req for an ideal gas.
    dW_dP = ForwardDiff.derivative(P -> W_for_power(0.88; P_req = P), 3.0e6)
    @test dW_dP ≈ W_for_power(0.88) / 3.0e6 rtol = 1e-8

    # A Float64-typed design field cannot hold the outer gradient's Duals.
    function bad(η)
        net, comp, heater = recup_net(η_c = η)
        add_free!(net, heater, :TtExit)
        add_target!(net, net_power, 3.0e6)
        solve!(net)
        net_power(net)
    end
    err = try ForwardDiff.derivative(bad, 0.88); nothing catch e; e end
    @test err isa ErrorException
    @test occursin("TtExit", err.msg) && occursin("one(eltype(p))", err.msg)
end

@testset "Targets and free variables — off-design" begin
    fluid = IdealGasFluid(M_molar = 83.8)
    P1, T1, W = 2.0e5, 400.0, 1.0
    PR_c, η_c, η_t, TIT, N_des = 2.0, 0.85, 0.90, 1100.0, 40_000.0

    # Design point (fixed-PR chain, bisect turbine PR for shaft balance).
    comp_d = Compressor("Comp"; PR = PR_c, η_poly = η_c)
    heat_d = HeatSource("Reactor"; TtExit = TIT, dPqP = 0.02)
    p2 = compute!(comp_d, Port(FluidState(P1, T1, W, fluid)))
    p3 = compute!(heat_d, p2)
    w_comp = specific_work(comp_d)
    turb_work(PR) = (t = Turbine("T"; PR = PR, η_poly = η_t); compute!(t, Port(p3[])); specific_work(t))
    lo, hi = 1.2, 4.0
    for _ in 1:60
        mid = (lo + hi) / 2
        turb_work(mid) < w_comp ? (lo = mid) : (hi = mid)
    end
    PR_t = (lo + hi) / 2
    cmap = synthetic_comp_map(Nc_des = corrected_speed(N_des, T1), Wc_des = corrected_flow(W, T1, P1),
                              PR_des = PR_c, eta_des = η_c)
    tmap = synthetic_turb_map(Np_des = corrected_speed(N_des, p3[].Tt),
                              Wp_des = corrected_flow(W, p3[].Tt, p3[].Pt),
                              PR_des = PR_t, eta_des = η_t)

    function build(; TtExit = TIT, N0 = 0.95 * N_des)
        comp  = Compressor("Comp"; η_poly = η_c, map = cmap, mode = :off_design)
        heat  = HeatSource("Reactor"; TtExit = TtExit, dPqP = 0.02)
        turb  = Turbine("Turb"; PR = PR_t, η_poly = η_t, map = tmap, mode = :off_design)
        shaft = Shaft("Main"; N = N0, mode = :off_design)
        net = FlowNetwork()
        add!(net, comp, heat, turb)
        connect!(net, comp => heat => turb)
        add_shaft!(net, shaft; drives = comp, driven_by = turb)
        set_state!(net, comp; Pt = P1, Tt = T1, W = W, fluid = fluid)
        (net, comp, heat, turb, shaft)
    end

    @testset "TIT to hold design shaft speed" begin
        net, comp, heat, turb, shaft = build(TtExit = 1000.0)
        add_free!(net, heat, :TtExit)
        add_target!(net, n -> shaft.N, N_des; name = "speed")
        sol = solve!(net)
        @test sol.status == :success
        @test shaft.N ≈ N_des rtol = 1e-6
        @test heat.TtExit ≈ TIT rtol = 1e-3       # maps pass through the design point
        @test comp.PR ≈ PR_c rtol = 1e-3
    end

    @testset "off-design solve is differentiable" begin
        # Shaft speed as a function of TIT (no targets): AD vs central differences.
        function speed(T)
            net, comp, heat, turb, shaft = build(TtExit = T)
            solve!(net; tol = 1e-10)
            shaft.N
        end
        d_ad = ForwardDiff.derivative(speed, 1050.0)
        d_fd = (speed(1050.0 + 0.01) - speed(1050.0 - 0.01)) / 0.02
        @test d_ad > 0
        @test d_ad ≈ d_fd rtol = 1e-5

        # With a target: TIT needed for a given speed, d(TIT)/d(N_req).
        function tit_for_speed(Nreq)
            net, comp, heat, turb, shaft = build(TtExit = 1000.0 * one(Nreq))
            add_free!(net, heat, :TtExit)
            add_target!(net, n -> shaft.N, Nreq)
            solve!(net; tol = 1e-10)
            heat.TtExit
        end
        g_ad = ForwardDiff.derivative(tit_for_speed, 0.98 * N_des)
        g_fd = (tit_for_speed(0.98 * N_des + 1.0) - tit_for_speed(0.98 * N_des - 1.0)) / 2
        @test g_ad > 0
        @test g_ad ≈ g_fd rtol = 1e-5
    end
end
