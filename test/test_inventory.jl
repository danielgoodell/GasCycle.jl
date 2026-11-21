using Test
using GasCycle
using ForwardDiff

isdefined(@__MODULE__, :synthetic_comp_map) ||
    include(joinpath(@__DIR__, "synthetic_maps.jl"))

const _R_HEXE84 = 8.314462618 / 0.0838

@testset "Volume mass integration" begin
    fluid = IdealGasFluid(M_molar = 83.8)
    # Heating 400 → 800 K at constant pressure: linear T, exact m = P·V/(R·T_lm).
    function heated(n)
        h = HeatSource("H"; TtExit = 800.0, dPqP = 0.0)
        net = FlowNetwork(); add!(net, h)
        set_state!(net, h; Pt = 500e3, Tt = 400.0, W = 1.0, fluid = fluid)
        add_volume!(net, h; V = 0.01, n = n)
        solve!(net)
        volume_mass(net.volumes[1])
    end
    T_lm = (800 - 400) / log(2)
    exact = 500e3 * 0.01 / (_R_HEXE84 * T_lm)
    @test heated(5) ≈ exact rtol = 1e-6
    @test heated(6) ≈ exact rtol = 1e-8
    errs = [abs(heated(n) - exact) for n in 1:6]
    @test issorted(errs; rev = true)                   # converges with points
    @test errs[1] / exact > 1e-2                       # midpoint alone is poor at 2:1

    # Uniform state: m = ρ·V exactly; HX sides are separate volumes.
    d = Duct("D"; dPqP = 0.0)
    net = FlowNetwork(); add!(net, d)
    set_state!(net, d; Pt = 300e3, Tt = 350.0, W = 1.0, fluid = fluid)
    add_volume!(net, d; V = 0.2)
    solve!(net)
    @test total_mass(net) ≈ 0.2 * density(fluid, 350.0, 300e3) rtol = 1e-12
    @test_throws ErrorException add_volume!(net, d; V = 1.0)              # twice
    @test_throws ErrorException add_volume!(net, Duct("X"); V = 1.0)      # not in network
    hx = HeatExchanger("HX"; ε = 0.9)
    add!(net, hx)
    @test_throws ErrorException add_volume!(net, hx; V = 1.0)             # needs side
end

# Closed loop with fixed temperatures: compressor and turbine on maps at fixed
# speed, cooler modeled as a negative-duty heat source to a fixed exit T.
function inventory_loop(; offdesign = nothing, W = 1.0, P1 = 500e3)
    fluid = IdealGasFluid(M_molar = 83.8)
    dh, dc, dq = 0.02, 0.01, 0.01
    if isnothing(offdesign)
        comp = Compressor("Comp"; PR = 2.0, η_poly = 0.85)
        turb = Turbine("Turb"; mode = :pressure_closure, P_exit = P1 / ((1 - dh) * (1 - dq)), η_poly = 0.90)
    else
        comp = Compressor("Comp"; η_poly = 0.85, map = offdesign.cmap, mode = :off_design)
        turb = Turbine("Turb"; PR = offdesign.PR_t, η_poly = 0.90, map = offdesign.tmap, mode = :off_design)
    end
    recup = HeatExchanger("Recup"; ε = 0.92, dPqP_hot = dh, dPqP_cold = dc)
    heat  = HeatSource("Heater"; TtExit = 1100.0, dPqP = 0.02)
    cool  = HeatSource("Cooler"; TtExit = 400.0, dPqP = dq)          # rejects heat
    shaft = Shaft("Shaft"; N = 36000.0)                               # alternator-locked speed
    net = FlowNetwork()
    add!(net, comp, recup, heat, turb, cool)
    connect!(net, comp => recup => heat => turb)
    add_hx_pair!(net, recup; hot = turb)
    connect_port!(net, recup, :hot_outlet, cool, :inlet)
    connect_port!(net, cool, :outlet, comp, :inlet; back_edge = true)
    add_shaft!(net, shaft; drives = comp, driven_by = turb)
    set_state!(net, comp; Pt = P1, Tt = 400.0, W = W, fluid = fluid)
    add_volume!(net, comp; V = 0.002)
    add_volume!(net, recup; V = 0.02, side = :cold)
    add_volume!(net, recup; V = 0.02, side = :hot)
    add_volume!(net, heat; V = 0.01)
    add_volume!(net, turb; V = 0.003)
    add_volume!(net, cool; V = 0.01)
    (; net, comp, turb, recup, cool)
end

@testset "Inventory control" begin
    d = inventory_loop()
    solve!(d.net; tol = 1e-12)
    s1, s3 = d.comp.inlet[], d.turb.inlet[]
    od = (cmap = synthetic_comp_map(Nc_des = corrected_speed(36000.0, s1.Tt),
                                    Wc_des = corrected_flow(s1.W, s1.Tt, s1.Pt),
                                    PR_des = 2.0, eta_des = 0.85),
          tmap = synthetic_turb_map(Np_des = corrected_speed(36000.0, s3.Tt),
                                    Wp_des = corrected_flow(s3.W, s3.Tt, s3.Pt),
                                    PR_des = pressure_ratio(d.turb), eta_des = 0.90),
          PR_t = pressure_ratio(d.turb))

    # Off-design at the design flow reproduces the design point; record its charge.
    m0 = inventory_loop(; offdesign = od)
    sol0 = solve!(m0.net; tol = 1e-12)
    @test sol0.status == :success
    @test m0.comp.inlet[].Pt ≈ 500e3 rtol = 1e-6
    M0 = loop_mass(m0.net, m0.comp)
    P0 = net_power(sol0)
    @test M0 ≈ total_mass(m0.net)
    @test occursin("total", sprint(inventory_report, sol0))

    # Same charge, wrong flow guess: inventory mode recovers the design flow.
    m1 = inventory_loop(; offdesign = od, W = 0.7)
    set_inventory!(m1.net, m1.comp; mass = M0)
    @test solve!(m1.net; tol = 1e-12).status == :success
    @test m1.comp.inlet[].W ≈ 1.0 rtol = 1e-7
    @test loop_mass(m1.net, m1.comp) ≈ M0 rtol = 1e-9

    # 80 % charge: corrected conditions unchanged, so flow, pressures and power
    # scale by 0.8 and temperatures stay put — inventory control.
    m2 = inventory_loop(; offdesign = od, W = 1.0)
    set_inventory!(m2.net, m2.comp; mass = 0.8 * M0)
    sol2 = solve!(m2.net; tol = 1e-12)
    @test sol2.status == :success
    @test m2.comp.inlet[].W ≈ 0.8 rtol = 1e-6
    @test m2.comp.inlet[].Pt ≈ 0.8 * 500e3 rtol = 1e-6
    @test m2.turb.inlet[].Tt ≈ m0.turb.inlet[].Tt rtol = 1e-8
    @test net_power(sol2) ≈ 0.8 * P0 rtol = 1e-6
    @test abs(loss_ledger(sol2).closure) < 1e-8 * loss_ledger(sol2).Q_in

    # Differentiable: dW/dM = W/M for this scale-invariant loop.
    function W_of(M)
        m = inventory_loop(; offdesign = od)
        set_inventory!(m.net, m.comp; mass = M)
        solve!(m.net; tol = 1e-12)
        m.comp.inlet[].W
    end
    @test ForwardDiff.derivative(W_of, 0.9 * M0) ≈ W_of(0.9 * M0) / (0.9 * M0) rtol = 1e-6

    # A loop seeded only through a boundary guess is found too.
    fluid = IdealGasFluid(M_molar = 83.8)
    a, b = Duct("A"; dPqP = 0.0), Duct("B"; dPqP = 0.0)
    n2 = FlowNetwork(); add!(n2, a, b)
    set_state!(n2, a; Pt = 1e5, Tt = 300.0, W = 1.0, fluid = fluid)
    set_boundary!(n2, b, :inlet; Pt = 2e5, Tt = 300.0, W = 0.5, fluid = fluid)
    get, set! = GasCycle._loop_seed(n2, b, :W)
    @test get() == 0.5
    set!(0.25)
    @test n2.boundaries[1][3].W == 0.25 && n2.seed_state.W == 1.0
end

@testset "Valve" begin
    fluid = IdealGasFluid(M_molar = 83.8)
    function valve_net(opening)
        v = Valve("V"; A = 1e-3, Cd = 0.8, opening = opening)
        net = FlowNetwork(); add!(net, v)
        set_state!(net, v; Pt = 300e3, Tt = 400.0, W = 0.5, fluid = fluid)
        (net, v)
    end
    net, v = valve_net(0.5)
    solve!(net)
    s = v.inlet[]
    ρ = density(s)
    @test s.Pt - v.outlet[].Pt ≈ 0.5^2 / (2ρ * (0.8 * 1e-3 * 0.5)^2) rtol = 1e-12
    @test v.outlet[].Tt == s.Tt
    # Opening as a control: throttle to a 230 kPa outlet (fully open gives ~274 kPa).
    add_free!(net, v, :opening; lower = 0.01, upper = 1.0)
    add_target!(net, n -> v.outlet[].Pt, 230e3)
    @test solve!(net).status == :success
    @test v.outlet[].Pt ≈ 230e3 rtol = 1e-8
    @test 0.5 < v.opening < 1.0
    # Overdemand → infeasible (NaN), not a silent negative pressure.
    net3, v3 = valve_net(0.01)
    one_pass!(net3)
    @test isnan(v3.outlet[].Pt)
end
