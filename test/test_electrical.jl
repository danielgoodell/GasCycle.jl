using Test
using GasCycle
using ForwardDiff

isdefined(@__MODULE__, :synthetic_comp_map) ||
    include(joinpath(@__DIR__, "synthetic_maps.jl"))

const _ef! = GasCycle._e_forward!
const _eb! = GasCycle._e_backward!

@testset "Electrical elements" begin
    @testset "DC cable" begin
        c = Cable("C"; R = 0.1)
        _ef!(c, 100.0, 0.0)
        P_in = _eb!(c, 500.0)
        I = c.e.I
        @test 100I - 2 * 0.1 * I^2 ≈ 500.0 rtol = 1e-12
        @test P_in - 500.0 ≈ 2 * 0.1 * I^2 rtol = 1e-12
        _ef!(c, 100.0, 0.0)                       # voltage sweep with the new current
        @test c.e.V_out ≈ 100.0 - 2 * 0.1 * I
        @test c.e.V_out * I ≈ 500.0 rtol = 1e-12  # delivered power = V_out·I
        @test isnan(_eb!(c, 1e6))                 # more than the line can carry
    end

    @testset "AC cable and cable from geometry" begin
        c = Cable("C"; kind = :ac, length = 20.0, area = 10e-6)
        @test c.R ≈ 1.72e-8 * 20 / 10e-6
        hot = Cable("C"; kind = :ac, length = 20.0, area = 10e-6, T = 393.15)
        @test hot.R ≈ c.R * (1 + 0.00393 * 100)
        _ef!(c, 400.0, 1200.0)
        _eb!(c, 10e3)
        _ef!(c, 400.0, 1200.0)
        @test sqrt(3) * c.e.V_out * c.e.I ≈ 10e3 rtol = 1e-12
        @test c.e.P_loss ≈ 3 * c.e.I^2 * c.R rtol = 1e-12
        @test c.e.f == 1200.0
    end

    @testset "converters" begin
        r = Rectifier("R"; η = 0.97)
        @test _ef!(r, 400.0, 1200.0) == (1.35 * 400.0, 0.0)  # passive bridge, DC out
        @test _eb!(r, 9700.0) ≈ 9700.0 / 0.97
        d = DCDCConverter("D"; V_out = 120.0, P_noload = 50.0, k_lin = 0.02, k_quad = 1e-6)
        @test _ef!(d, 540.0, 0.0)[1] == 120.0
        @test _eb!(d, 1000.0) ≈ 1000 + 50 + 20 + 1
        t = Transformer("T"; ratio = 3.0, η = 0.99)
        @test _ef!(t, 400.0, 1200.0) == (1200.0, 1200.0)       # AC keeps frequency
        box = PowerConverter("Box"; V_out = 28.0, η = 0.9)
        @test (GasCycle.e_kind_in(box), GasCycle.e_kind_out(box)) == (:ac, :dc)
        @test_throws ErrorException Converter("X"; input = :ac, output = :dc)   # no V_out/ratio
        @test_throws ErrorException Rectifier("X"; η = 1.5)
    end

    @testset "wiring checks" begin
        net = FlowNetwork()
        alt = Alternator("A"; V_rated = 400.0, N_rated = 36000.0, η = 0.95)
        @test_throws ErrorException connect_electrical!(net, alt => Cable("dc"; R = 0.1))   # :ac → :dc
        load = ElectricalLoad("L"; P = 1.0, kind = :ac)
        @test_throws ErrorException connect_electrical!(net, load => Cable("x"; kind = :ac, R = 0.1))
        bus = ElectricalBus("B"; kind = :ac)
        connect_electrical!(net, alt => bus)
        @test_throws ErrorException connect_electrical!(net, Alternator("A2"; V_rated = 1.0, N_rated = 1.0, η = 0.9) => bus)
    end
end

# Recuperated cycle, alternator cooled by the compressor-inlet gas (coldest
# point of the loop), and the
# electrical chain alt → feeder → rectifier → DC link → {PLR, DC-DC → bus → user}.
function elec_cycle(; offdesign = nothing, W = 1.0, user_P = 15e3, plr_P = 1000.0,
                     η_rect = 0.97, N = 36000.0)
    fluid  = HeXeIdealGas(0.47)
    T0, P0 = 400.0, 500e3
    if isnothing(offdesign)
        comp = Compressor("C"; PR = 2.5, η_poly = 0.88)
        turb = Turbine("T"; mode = :pressure_closure, P_exit = P0 / (0.99 * 0.99), η_poly = 0.90)  # closes to the alternator inlet
        shaft = Shaft("Shaft"; N = N)
    else
        comp = Compressor("C"; η_poly = 0.88, map = offdesign.cmap, mode = :off_design)
        turb = Turbine("T"; PR = offdesign.PR_t, η_poly = 0.90, map = offdesign.tmap, mode = :off_design)
        shaft = Shaft("Shaft"; N = 0.97 * N, mode = :off_design)
    end
    alt    = Alternator("Alt"; V_rated = 400.0, N_rated = 36000.0, poles = 4, R_phase = 0.02,
                        P_core_rated = 300.0, P_windage_rated = 200.0, dPqP = 0.005)
    recup  = HeatExchanger("R"; ε = 0.90, dPqP_hot = 0.01, dPqP_cold = 0.01)
    heater = HeatSource("H"; TtExit = 1100.0, dPqP = 0.02)
    net = FlowNetwork()
    add!(net, comp, alt, recup, heater, turb)
    connect!(net, alt => comp => recup => heater => turb => alt)
    add_hx_pair!(net, recup; hot = turb)
    add_shaft!(net, shaft; drives = [comp, alt], driven_by = turb)
    set_state!(net, alt; Pt = P0, Tt = T0, W = W, fluid = fluid)

    feeder = Cable("Feeder"; kind = :ac, length = 20.0, area = 10e-6)
    rect   = Rectifier("Rect"; η = η_rect)
    link   = ElectricalBus("DCLink")
    plr    = ParasiticLoad("PLR"; P = plr_P)
    dcdc   = DCDCConverter("DCDC"; V_out = 120.0, η = 0.95, P_noload = 50.0)
    bus    = ElectricalBus("Bus120")
    user   = ElectricalLoad("User"; P = user_P)
    connect_electrical!(net, alt => feeder => rect => link)
    connect_electrical!(net, link => plr)
    connect_electrical!(net, link => dcdc => bus => user)
    (; net, comp, turb, alt, recup, heater, shaft, feeder, rect, link, plr, dcdc, user)
end

@testset "Electrical network in the cycle — design" begin
    m = elec_cycle()
    add_free!(m.net, :W)
    add_target!(m.net, n -> power_balance(m.shaft), 0.0; scale = 1e3)
    sol = solve!(m.net; tol = 1e-10)
    @test sol.status == :success
    @test abs(power_balance(m.shaft)) < 1e-3                  # W
    @test shaft_power(m.alt) ≈ net_power(sol) rtol = 1e-8     # turbine − compressor feeds the alternator

    # Electrical bookkeeping
    @test electrical_frequency(m.alt) ≈ 4 * 36000 / 120
    @test m.alt.e.V_out ≈ 400.0 - sqrt(3) * m.alt.e.I * 0.02 rtol = 1e-12
    @test m.link.e.V_out ≈ 1.35 * m.feeder.e.V_out rtol = 1e-12
    @test m.user.e.V_in == 120.0
    @test electrical_power(sol) == 15e3
    @test parasitic_power(sol) == 1000.0
    eloss = m.feeder.e.P_loss + m.rect.e.P_loss + m.dcdc.e.P_loss
    @test m.alt.e.P_out ≈ 15e3 + 1000 + eloss rtol = 1e-12
    @test m.alt.P_mech ≈ m.alt.e.P_out + m.alt.e.P_loss rtol = 1e-12

    # Alternator losses heat the gas
    a = m.alt
    @test a.inlet[].W * (enthalpy(a.outlet[]) - enthalpy(a.inlet[])) ≈ a.e.P_loss rtol = 1e-9
    @test a.outlet[].Pt ≈ a.inlet[].Pt * (1 - 0.005)

    # First and second law with the electrical chain
    eb = energy_balance(sol)
    @test eb.P_alternator ≈ a.P_mech
    @test eb.P_delivered == 15e3 && eb.P_parasitic == 1000.0
    @test eb.P_elec_losses ≈ eloss rtol = 1e-12
    @test abs(eb.shaft_surplus) < 1e-3
    @test abs(eb.residual) < 1e-9 * eb.Q_in
    L = loss_ledger(sol)
    @test L.electrical && L.P_delivered == 15e3
    @test abs(L.closure) < 1e-9 * L.Q_in
    kinds = first.(losses_by_kind(L))
    @test :electrical in kinds && :parasitic in kinds
    alt_destroyed = sum(l.exergy for l in L.losses if l.element == "Alt" && l.kind === :electrical)
    @test 0 < alt_destroyed < a.e.P_loss          # part of the loss heat is recovered by the cycle
    @test occursin("delivered electric", sprint(show, MIME"text/plain"(), L))

    rep = sprint(electrical_report, sol)
    @test occursin("DCDC", rep) && occursin("ParasiticLoad", rep)
    rs = sprint(summary, sol)
    @test occursin("P_elec=", rs) && occursin("Electric delivered", rs)
end

@testset "Electrical network in the cycle — PLR speed control (off-design)" begin
    # Design point, then maps through it.
    d = elec_cycle()
    add_free!(d.net, :W)
    add_target!(d.net, n -> power_balance(d.shaft), 0.0; scale = 1e3)
    solve!(d.net; tol = 1e-10)
    W_des = d.comp.inlet[].W
    s1, s3 = d.comp.inlet[], d.turb.inlet[]
    cmap = synthetic_comp_map(Nc_des = corrected_speed(36000.0, s1.Tt),
                              Wc_des = corrected_flow(s1.W, s1.Tt, s1.Pt),
                              PR_des = 2.5, eta_des = 0.88)
    PR_t = pressure_ratio(d.turb)
    tmap = synthetic_turb_map(Np_des = corrected_speed(36000.0, s3.Tt),
                              Wp_des = corrected_flow(s3.W, s3.Tt, s3.Pt),
                              PR_des = PR_t, eta_des = 0.90)
    od = (; cmap, tmap, PR_t)

    function run(user_P; plr_guess = 500.0)
        mm = elec_cycle(; offdesign = od, W = W_des, user_P = user_P, plr_P = plr_guess)
        add_speed_control!(mm.net, mm.plr, mm.shaft; N = 36000.0, ref = 1000.0)
        (solve!(mm.net; tol = 1e-10), mm)
    end

    sol, m = run(15e3)
    @test sol.status == :success
    @test m.shaft.N ≈ 36000.0 rtol = 1e-8
    @test m.plr.P ≈ 1000.0 rtol = 1e-4          # reproduces the design split
    @test abs(power_balance(m.shaft)) < 1e-3

    # Shed 1 kW of user load: same speed and gas states, so the PLR takes up
    # exactly what the DC-DC stage no longer draws.
    sol2, m2 = run(14e3)
    @test sol2.status == :success
    @test m2.shaft.N ≈ 36000.0 rtol = 1e-8
    @test m2.plr.P - m.plr.P ≈ 1000 / 0.95 rtol = 1e-6

    # Asking for more than the unit can make drives the PLR to its bound.
    sol3, m3 = @test_logs (:warn,) match_mode = :any run(20e3)
    @test sol3.status == :failed
    @test m3.plr.P ≈ 0.0 atol = 1e-9
end

@testset "Electrical network — differentiable" begin
    # Mass flow needed to balance the shaft, as a function of rectifier efficiency.
    function W_for(η)
        m = elec_cycle(; η_rect = η)
        add_free!(m.net, :W)
        add_target!(m.net, n -> power_balance(m.shaft), 0.0; scale = 1e3)
        solve!(m.net; tol = 1e-12)
        m.comp.inlet[].W
    end
    d_ad = ForwardDiff.derivative(W_for, 0.97)
    d_fd = (W_for(0.97 + 1e-6) - W_for(0.97 - 1e-6)) / 2e-6
    @test d_ad < 0                     # better rectifier → less flow needed
    @test d_ad ≈ d_fd rtol = 1e-5
end
