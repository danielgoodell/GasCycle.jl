using Test
using GasCycle
using ForwardDiff
import GasCycle: cp

const _ACC_DATA = joinpath(@__DIR__, "..", "data")

# Open-loop recuperated cycle: the compressor inlet is seeded and the
# recuperator hot outlet dangles (implicit cooler).
function acc_open_loop(fluid; η_c = 0.88)
    T0, P0 = 400.0, 500e3
    net    = FlowNetwork()
    comp   = Compressor("C"; PR = 2.5, η_poly = η_c)
    recup  = HeatExchanger("R"; ε = 0.90, dPqP_hot = 0.01, dPqP_cold = 0.01)
    heater = HeatSource("H"; TtExit = 1100.0, dPqP = 0.02)
    turb   = Turbine("T"; mode = :pressure_closure, P_exit = P0 / (0.99 * 0.99), η_poly = 0.90)
    add!(net, comp, recup, heater, turb)
    connect!(net, comp => recup => heater => turb => comp)
    add_hx_pair!(net, recup; hot = turb)
    set_state!(net, comp; Pt = P0, Tt = T0, W = 10.0, fluid = fluid)
    (net, comp, recup, heater, turb)
end

# Closed loop through a water cooler or a radiator (as in test_cooler.jl).
function acc_closed_loop(cold_end)
    fluid = IdealGasFluid(M_molar = 83.8)
    water = ConstantPropertyLiquid(joinpath(_ACC_DATA, "H2O.fpt"))
    P1, W = 163.4e3, 0.6
    dPqP_hot, dPqP_cool = 0.022, 0.01
    comp  = Compressor("Comp"; PR = 1.9, η_poly = 0.80)
    recup = HeatExchanger("Recup"; ε = 0.95, dPqP_hot = dPqP_hot, dPqP_cold = 0.011)
    heat  = HeatSource("Heater"; TtExit = 1144.0, dPqP = 0.027)
    turb  = Turbine("Turb"; η_poly = 0.87, mode = :pressure_closure,
                    P_exit = P1 / ((1 - dPqP_hot) * (1 - dPqP_cool)))
    net = FlowNetwork()
    add!(net, comp, recup, heat, turb, cold_end)
    connect!(net, comp => recup => heat => turb)
    add_hx_pair!(net, recup; hot = turb)
    if cold_end isa HeatExchanger
        connect_port!(net, recup, :hot_outlet, cold_end, :hot_inlet)
        connect_port!(net, cold_end, :hot_outlet, comp, :inlet; back_edge = true)
        set_boundary!(net, cold_end, :cold_inlet; Pt = 600e3, Tt = 290.0, W = 1.0, fluid = water)
    else
        connect_port!(net, recup, :hot_outlet, cold_end, :inlet)
        connect_port!(net, cold_end, :outlet, comp, :inlet; back_edge = true)
    end
    set_state!(net, comp; Pt = P1, Tt = 300.0, W = W, fluid = fluid)
    net
end

lossof(L, el, kind) = sum((l.exergy for l in L.losses if l.element == el && l.kind === kind); init = 0.0)

@testset "Energy balance and loss ledger — open loop" begin
    for fluid in (HeXeIdealGas(0.47), HeXe(40.0))   # ideal and real gas
        net, comp, recup, heater, turb = acc_open_loop(fluid)
        sol = solve!(net)

        eb = energy_balance(sol)
        @test eb.P_shaft ≈ net_power(sol) rtol = 1e-12
        @test eb.Q_in ≈ heater.Q rtol = 1e-10
        @test eb.Q_out == 0
        @test abs(eb.imbalance) < 1e-9 * eb.Q_in        # HX works in enthalpy: exact
        @test abs(eb.residual)  < 1e-9 * eb.Q_in

        L = loss_ledger(sol)
        @test L.T0 == 400.0                              # seed inlet is the coldest point
        @test L.P_shaft ≈ net_power(sol) rtol = 1e-12
        @test abs(L.closure) < 1e-9 * L.Q_in            # Gouy–Stodola identity
        Ex = sum(last, L.exergy_in)
        @test 0 < Ex < L.Q_in                            # exergy < heat
        @test Ex ≈ heater.Q * (1 - L.T0 / 800) rtol = 0.2  # Carnot-ish sanity
        for (el, kind) in (("C", :aero), ("T", :aero), ("R", :heat_transfer),
                           ("R", :pressure), ("H", :pressure))
            @test lossof(L, el, kind) > 0               # second law
        end
        # The open loop's implicit cooler is one netted rejection line.
        rej = [l for l in L.losses if l.kind === :heat_rejection]
        @test length(rej) == 1
        @test rej[1].element == "R.hot_outlet"
        @test rej[1].exergy > 0
        # Grouping and display
        bk = losses_by_kind(L)
        @test first.(bk) ⊇ [:aero, :pressure, :heat_transfer, :heat_rejection]
        @test sum(last, bk) ≈ Ex - L.P_shaft rtol = 1e-9
        txt = sprint(show, MIME"text/plain"(), L)
        @test occursin("net shaft power", txt) && occursin("C  (aero)", txt)

        # Dead state is user-selectable; a colder T0 values the rejected heat more.
        L2 = loss_ledger(sol; T0 = 300.0)
        @test abs(L2.closure) < 1e-9 * L2.Q_in
        @test sum(last, L2.exergy_in) > Ex
    end
end

@testset "Energy balance and loss ledger — closed loops" begin
    @testset "water cooler (boundary coolant stream)" begin
        cool = HeatExchanger("Cooler"; UA = 400.0, dPqP_hot = 0.01, dPqP_cold = 0.005)
        net = acc_closed_loop(cool)
        sol = solve!(net; tol = 1e-10)
        eb = energy_balance(sol)
        @test eb.Q_open ≈ Q_transferred(cool) rtol = 1e-8   # coolant carries the rejection
        @test abs(eb.residual) < 1e-7 * eb.Q_in
        L = loss_ledger(sol)
        @test L.T0 == 290.0                                  # coolant inlet
        @test abs(L.closure) < 1e-7 * L.Q_in
        @test lossof(L, "Cooler", :heat_transfer) > 0
        rej = [l for l in L.losses if l.kind === :heat_rejection]
        @test length(rej) == 1 && rej[1].element == "Cooler.cold_outlet"
    end

    @testset "radiator" begin
        rad = Radiator("Rad"; mode = :fixed_TtExit, TtExit = 320.0, T_sink = 230.0, dPqP = 0.01)
        net = acc_closed_loop(rad)
        sol = solve!(net; tol = 1e-10)
        eb = energy_balance(sol)
        @test eb.Q_open ≈ 0 atol = 1e-6 * eb.Q_in           # no open streams
        @test eb.Q_out ≈ eb.Q_in - eb.P_shaft rtol = 1e-7
        L = loss_ledger(sol)
        @test L.T0 == 230.0                                  # radiator sink
        @test abs(L.closure) < 1e-7 * L.Q_in
        @test lossof(L, "Rad", :heat_rejection) > 0
        @test lossof(L, "Rad", :pressure) > 0
    end
end

@testset "Loss ledger — mixing" begin
    # Split a stream, heat one branch, remix: the mixer generates entropy.
    fluid = IdealGasFluid(M_molar = 83.8)
    d0  = Duct("In"; dPqP = 0.0)
    spl = Splitter("S"; fracs = [0.6, 0.4])
    hot = HeatSource("Hot"; TtExit = 700.0, dPqP = 0.0)
    byp = Duct("Bypass"; dPqP = 0.0)
    mix = Mixer("M")
    net = FlowNetwork()
    add!(net, d0, spl, hot, byp, mix)
    connect!(net, d0 => spl)
    connect_port!(net, spl, :outlet, hot, :inlet)
    connect_port!(net, spl, :bleed_outlet, byp, :inlet)
    connect_port!(net, hot, :outlet, mix, :inlet)
    connect_port!(net, byp, :outlet, mix, :bleed_inlet)
    set_state!(net, d0; Pt = 300e3, Tt = 400.0, W = 1.0, fluid = fluid)
    sol = solve!(net)
    L = loss_ledger(sol)
    @test lossof(L, "M", :mixing) > 0
    @test abs(lossof(L, "S", :mixing)) < 1e-9
    @test abs(L.closure) < 1e-9 * L.Q_in
    @test L.P_shaft == 0
end

@testset "Loss ledger — HeatSource with negative duty is a rejection" begin
    fluid = IdealGasFluid(M_molar = 83.8)
    heat = HeatSource("Heat"; TtExit = 900.0, dPqP = 0.01)
    cool = HeatSource("Cool"; TtExit = 500.0, dPqP = 0.01)    # cools: Q < 0
    net = FlowNetwork()
    add!(net, heat, cool)
    connect!(net, heat => cool)
    set_state!(net, heat; Pt = 300e3, Tt = 400.0, W = 1.0, fluid = fluid)
    sol = solve!(net)
    eb = energy_balance(sol)
    @test eb.Q_in ≈ heat.Q rtol = 1e-12
    @test eb.Q_out ≈ -cool.Q rtol = 1e-12
    L = loss_ledger(sol)
    @test first.(L.exergy_in) == ["Heat"]
    @test lossof(L, "Cool", :heat_rejection) > 0
    @test abs(L.closure) < 1e-9 * L.Q_in
end

@testset "Loss ledger — differentiable" begin
    aero(η) = (s = acc_open_loop(HeXeIdealGas(0.47); η_c = η);
               loss_ledger(solve!(s[1]; tol = 1e-12)).losses |> ls -> sum(l.exergy for l in ls if l.element == "C"))
    d_ad = ForwardDiff.derivative(aero, 0.88)
    d_fd = (aero(0.88 + 1e-5) - aero(0.88 - 1e-5)) / 2e-5
    @test d_ad < 0                    # better compressor → less loss
    @test d_ad ≈ d_fd rtol = 1e-5
end
