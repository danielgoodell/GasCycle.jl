using Test
using GasCycle
using ForwardDiff
import GasCycle: cp

const _TH_DATA = joinpath(@__DIR__, "..", "data")

Tmean(d) = (d.inlet[].Tt + d.outlet[].Tt) / 2

# Two independent streams: hot gas (heated, then a duct) and a cold boundary
# stream through its own duct.
function two_streams(; UA = 5.0, T_hot = 800.0)
    fluid = IdealGasFluid(M_molar = 83.8)
    heat = HeatSource("H"; TtExit = T_hot, dPqP = 0.0)
    dA = Duct("DuctA"; dPqP = 0.01)
    dB = Duct("DuctB"; dPqP = 0.01)
    net = FlowNetwork()
    add!(net, heat, dA, dB)
    connect!(net, heat => dA)
    set_state!(net, heat; Pt = 300e3, Tt = 400.0, W = 0.5, fluid = fluid)
    set_boundary!(net, dB, :inlet; Pt = 200e3, Tt = 300.0, W = 0.3, fluid = fluid)
    (; net, heat, dA, dB, fluid)
end

@testset "Heat links — stream to stream" begin
    m = two_streams()
    add_conductance!(m.net, m.dA, m.dB; UA = 5.0, name = "crossover wall")
    sol = solve!(m.net; tol = 1e-12)
    @test sol.status == :success
    l = m.net.heat_links[1]
    @test l.Q ≈ 5.0 * (Tmean(m.dA) - Tmean(m.dB)) rtol = 1e-10
    @test l.Q > 0
    hA = m.dA.inlet[].W * (enthalpy(m.dA.outlet[]) - enthalpy(m.dA.inlet[]))
    hB = m.dB.inlet[].W * (enthalpy(m.dB.outlet[]) - enthalpy(m.dB.inlet[]))
    @test hA ≈ -l.Q rtol = 1e-10
    @test hB ≈ l.Q rtol = 1e-10
    @test m.dA.outlet[].Pt ≈ m.dA.inlet[].Pt * 0.99           # pressure loss still applied

    eb = energy_balance(sol)
    @test abs(eb.imbalance) < 1e-9 * eb.Q_in
    @test abs(eb.residual) < 1e-9 * eb.Q_in
    L = loss_ledger(sol)
    @test abs(L.closure) < 1e-9 * L.Q_in
    xfer = sum(x.exergy for x in L.losses if x.element == "crossover wall")
    @test xfer > 0
    # Lumped estimate T0·Q·(1/T_B − 1/T_A) of the destroyed exergy
    @test xfer ≈ L.T0 * l.Q * (1 / Tmean(m.dB) - 1 / Tmean(m.dA)) rtol = 0.02
    @test occursin("crossover wall", sprint(heat_link_report, sol))

    @test_throws ErrorException add_conductance!(m.net, m.heat, m.dB; UA = 1.0)  # not a Duct
end

@testset "Heat links — environment nodes" begin
    m = two_streams()
    space  = ThermalBoundary("Space"; T = 4.0)
    struc  = ThermalBoundary("Structure"; T = 250.0)
    add_radiation!(m.net, m.dA, space; area = 0.05, emissivity = 0.8)
    add_conductance!(m.net, m.dA, struc; UA = 0.5)
    sol = solve!(m.net; tol = 1e-12)
    T = Tmean(m.dA)
    q_rad = 5.670374419e-8 * 0.8 * 0.05 * (T^4 - 4.0^4)
    q_con = 0.5 * (T - 250.0)
    env = Dict(environment_heat(sol))
    @test env["Space"] ≈ q_rad rtol = 1e-10
    @test env["Structure"] ≈ q_con rtol = 1e-10
    eb = energy_balance(sol)
    @test eb.Q_env ≈ q_rad + q_con rtol = 1e-10
    @test abs(eb.residual) < 1e-9 * eb.Q_in
    L = loss_ledger(sol)
    @test abs(L.closure) < 1e-9 * L.Q_in
    @test :heat_leak in first.(losses_by_kind(L))

    # An external thermal model updates the sink temperatures by name.
    set_environment!(m.net, Dict("Structure" => 350.0))
    solve!(m.net; tol = 1e-12)
    @test Dict(environment_heat(m.net))["Structure"] < env["Structure"]
    @test_throws ErrorException set_environment!(m.net, Dict("Nowhere" => 1.0))

    # Derivative of the leak to structure w.r.t. its temperature (AD vs FD).
    function leak(Ts)
        mm = two_streams()
        node = ThermalBoundary("Structure"; T = Ts)
        add_conductance!(mm.net, mm.dA, node; UA = 0.5)
        solve!(mm.net; tol = 1e-12)
        mm.net.heat_links[1].Q
    end
    d_ad = ForwardDiff.derivative(leak, 250.0)
    d_fd = (leak(250.01) - leak(249.99)) / 0.02
    @test d_ad < 0
    @test d_ad ≈ d_fd rtol = 1e-6
end

@testset "Radiator rejecting to an environment node" begin
    fluid = IdealGasFluid(M_molar = 83.8)
    node = ThermalBoundary("RadSink"; T = 230.0)
    heat = HeatSource("H"; TtExit = 500.0, dPqP = 0.0)
    rad  = Radiator("Rad"; A = 2.0, emissivity = 0.85, sink = node, dPqP = 0.01)
    net = FlowNetwork()
    add!(net, heat, rad)
    connect!(net, heat => rad)
    set_state!(net, heat; Pt = 300e3, Tt = 400.0, W = 0.3, fluid = fluid)
    solve!(net)
    T_cold = rad.outlet[].Tt
    @test Dict(environment_heat(net))["RadSink"] ≈ Q_rejected(rad) rtol = 1e-12
    set_environment!(net, ["RadSink" => 280.0])
    solve!(net)
    @test rad.outlet[].Tt > T_cold                     # warmer sink, less rejection
    @test loss_ledger(net).T0 == 280.0                 # dead state follows the node
end

# Gas cycle with an alternator, plus a pumped water loop that cools the power
# electronics (rectifier and DC-DC losses into a cold-plate duct).
function cooled_electronics(; wired = true)
    gas = HeXeIdealGas(0.47)
    water = ConstantPropertyLiquid(joinpath(_TH_DATA, "H2O.fpt"))
    P0 = 500e3
    alt    = Alternator("Alt"; V_rated = 400.0, N_rated = 36000.0, poles = 4, η = 0.95)
    comp   = Compressor("C"; PR = 2.5, η_poly = 0.88)
    recup  = HeatExchanger("R"; ε = 0.90, dPqP_hot = 0.01, dPqP_cold = 0.01)
    heater = HeatSource("H"; TtExit = 1100.0, dPqP = 0.02)
    turb   = Turbine("T"; mode = :pressure_closure, P_exit = P0 / (0.99 * 0.99), η_poly = 0.90)
    shaft  = Shaft("Shaft"; N = 36000.0)
    pump   = Pump("Pump"; ΔP = 50e3, η = 0.4)
    plate  = Duct("ColdPlate"; dPqP = 0.0)
    net = FlowNetwork()
    add!(net, alt, comp, recup, heater, turb, pump, plate)
    connect!(net, alt => comp => recup => heater => turb => alt)
    add_hx_pair!(net, recup; hot = turb)
    add_shaft!(net, shaft; drives = [comp, alt], driven_by = turb)
    set_state!(net, alt; Pt = P0, Tt = 400.0, W = 1.0, fluid = gas)
    connect!(net, pump => plate)
    set_boundary!(net, pump, :inlet; Pt = 200e3, Tt = 290.0, W = 0.2, fluid = water)

    rect = Rectifier("Rect"; η = 0.97)
    link = ElectricalBus("DCLink")
    dcdc = DCDCConverter("DCDC"; V_out = 120.0, η = 0.95, P_noload = 50.0)
    bus  = ElectricalBus("Bus120")
    user = ElectricalLoad("User"; P = 15e3)
    connect_electrical!(net, alt => rect => link => dcdc => bus => user)
    wired && connect_electrical!(net, bus => pump)
    add_loss_heat!(net, rect, plate)
    add_loss_heat!(net, dcdc, plate; fraction = 0.8)
    add_free!(net, :W)
    add_target!(net, n -> power_balance(shaft), 0.0; scale = 1e3)
    (; net, alt, pump, plate, rect, dcdc, user, water)
end

@testset "Electronics cooling: loss heat into a pumped coolant loop" begin
    for wired in (true, false)
        m = cooled_electronics(; wired)
        sol = solve!(m.net; tol = 1e-10)
        @test sol.status == :success
        s = m.pump.inlet[]
        @test m.pump.P_elec ≈ s.W * 50e3 / (density(s) * 0.4) rtol = 1e-9
        @test m.pump.outlet[].Pt ≈ s.Pt + 50e3
        Q_in_plate = m.rect.e.P_loss + 0.8 * m.dcdc.e.P_loss
        @test m.plate.Q_ext ≈ Q_in_plate rtol = 1e-9
        rise = m.plate.outlet[].W * (enthalpy(m.plate.outlet[]) - enthalpy(m.pump.inlet[]))
        @test rise ≈ m.pump.P_elec + Q_in_plate rtol = 1e-8   # pump power + electronics heat
        if wired
            @test m.alt.e.P_out ≈ 15e3 + m.pump.P_elec + m.rect.e.P_loss + m.dcdc.e.P_loss rtol = 1e-9
        end

        eb = energy_balance(sol)
        @test eb.P_pumps ≈ m.pump.P_elec
        @test eb.P_pumps_external ≈ (wired ? 0.0 : m.pump.P_elec) atol = 1e-9
        @test eb.Q_loss_to_streams ≈ Q_in_plate rtol = 1e-9
        @test abs(eb.imbalance) < 1e-8 * eb.Q_in
        @test abs(eb.residual) < 1e-8 * eb.Q_in

        L = loss_ledger(sol)
        @test abs(L.closure) < 1e-8 * L.Q_in
        @test any(x -> x[1] == "Pump (external power)", L.exergy_in) == !wired
        rect_net = sum(x.exergy for x in L.losses if x.element == "Rect")
        @test 0 < rect_net < m.rect.e.P_loss             # part recovered as coolant exergy
        @test :pumping in first.(losses_by_kind(L))
    end
    @test_throws ErrorException add_loss_heat!(cooled_electronics().net,
                                               ElectricalLoad("X"; P = 1.0), Duct("D"))
end
