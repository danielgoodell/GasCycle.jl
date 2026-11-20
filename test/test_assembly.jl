using Test
using GasCycle

const _AS_DATA = joinpath(@__DIR__, "..", "data")

# Turbine–alternator–compressor on one shaft; gas enters at the alternator.
function tac_assembly()
    tac   = Assembly("TAC")
    alt   = Alternator("Alt"; V_rated = 400.0, N_rated = 36000.0, poles = 4, η = 0.95)
    comp  = Compressor("Comp"; PR = 1.9, η_poly = 0.80)
    turb  = Turbine("Turb"; η_poly = 0.87, mode = :pressure_closure,
                    P_exit = 163.4e3 / ((1 - 0.022) * (1 - 0.01)))
    shaft = Shaft("Shaft"; N = 36000.0)
    add!(tac, alt, comp, turb)
    connect!(tac, alt => comp)
    add_shaft!(tac, shaft; drives = [comp, alt], driven_by = turb)
    expose!(tac, :inlet, alt)
    expose!(tac, :comp_out, comp, :outlet)
    expose!(tac, :turb_in, turb, :inlet)
    expose!(tac, :outlet, turb)
    expose!(tac, :elec, alt)
    tac
end

# Power conversion unit: TAC + recuperator + heater + water cooler, closed
# loop, its own rectifier and user load sized to balance the shaft.
function pcu_assembly(name)
    gas   = IdealGasFluid(M_molar = 83.8)
    water = ConstantPropertyLiquid(joinpath(_AS_DATA, "H2O.fpt"))
    pcu   = Assembly(name)
    tac   = tac_assembly()
    recup = HeatExchanger("Recup"; ε = 0.95, dPqP_hot = 0.022, dPqP_cold = 0.011)
    heat  = HeatSource("Heater"; TtExit = 1144.0, dPqP = 0.027)
    cool  = HeatExchanger("Cooler"; UA = 400.0, dPqP_hot = 0.01, dPqP_cold = 0.005)
    add!(pcu, tac, recup, heat, cool)
    connect_port!(pcu, tac, :comp_out, recup, :inlet)
    connect!(pcu, recup => heat)
    connect_port!(pcu, heat, :outlet, tac, :turb_in)
    add_hx_pair!(pcu, recup; hot = tac)
    connect_port!(pcu, recup, :hot_outlet, cool, :hot_inlet)
    connect_port!(pcu, cool, :hot_outlet, tac, :inlet; back_edge = true)
    set_boundary!(pcu, cool, :cold_inlet; Pt = 600e3, Tt = 290.0, W = 1.0, fluid = water)
    # Seed guess for this loop (overridden by the loop closure once solved).
    set_boundary!(pcu, tac["Alt"], :inlet; Pt = 163.4e3, Tt = 300.0, W = 0.6, fluid = gas)

    rect = Rectifier("Rect"; η = 0.97)
    user = ElectricalLoad("User"; P = 5e3)
    connect_electrical!(pcu, tac => rect => user)
    shaft = tac["Shaft"]
    add_free!(pcu, user, :P)
    add_target!(pcu, n -> power_balance(shaft), 0.0; scale = 1e3, name = "shaft balance")
    expose!(pcu, :elec, user)
    pcu
end

@testset "Assemblies" begin
    pcu1, pcu2 = pcu_assembly("PCU1"), pcu_assembly("PCU2")
    net = FlowNetwork()
    add!(net, pcu1)
    add!(net, pcu2)
    # The network needs one seed; the second loop is seeded by its boundary guess.
    set_state!(net, pcu1["TAC"]; Pt = 163.4e3, Tt = 300.0, W = 0.6,
               fluid = IdealGasFluid(M_molar = 83.8))
    sol = solve!(net; tol = 1e-10)
    @test sol.status == :success

    @testset "names and lookup" begin
        c1 = sol["PCU1.TAC.Comp"]
        @test c1 === pcu1["TAC"]["Comp"] === pcu1["TAC.Comp"]
        @test c1 !== sol["PCU2.TAC.Comp"]
        @test pcu1["TAC.Shaft"].name == "PCU1.TAC.Shaft"
        @test pcu1["User"].name == "PCU1.User"
        @test assemblies(net) == ["PCU1", "PCU1.TAC", "PCU2", "PCU2.TAC"]
        m = members(net, "PCU1.TAC")
        @test length(m) == 4 && pcu1["TAC.Shaft"] in m       # alt, comp, turb, shaft
        @test any(fv -> fv.name == "PCU1.User.P", net.free_vars)
        @test any(t -> t.name == "PCU1.shaft balance", net.targets)
        @test_throws ErrorException pcu1["Nope"]
        @test occursin("PCU2.TAC.Alt", sprint(summary, sol))
    end

    @testset "identical units solve identically" begin
        u1, u2 = pcu1["User"], pcu2["User"]
        @test u1.P > 0
        @test u1.P ≈ u2.P rtol = 1e-8
        @test abs(power_balance(pcu1["TAC.Shaft"])) < 1e-3
        @test abs(power_balance(pcu2["TAC.Shaft"])) < 1e-3
        @test pcu1["TAC.Comp"].inlet[].Tt ≈ pcu2["TAC.Comp"].inlet[].Tt rtol = 1e-8
        @test electrical_power(sol) ≈ u1.P + u2.P rtol = 1e-12
    end

    @testset "bookkeeping by assembly" begin
        eb = energy_balance(sol)
        @test abs(eb.residual) < 1e-7 * eb.Q_in
        L = loss_ledger(sol)
        @test L.T0 == 290.0                                   # coolant inlets, not the seed guesses
        @test abs(L.closure) < 1e-7 * L.Q_in
        by1 = Dict(losses_by_assembly(L, net))
        @test by1["PCU1"] ≈ by1["PCU2"] rtol = 1e-7
        @test !haskey(by1, "(top level)")
        @test sum(values(by1)) ≈ sum(l.exergy for l in L.losses)
        by2 = Dict(losses_by_assembly(L, net; depth = 2))
        @test haskey(by2, "PCU1.TAC") && haskey(by2, "PCU1")
        @test by2["PCU1.TAC"] + by2["PCU1"] ≈ by1["PCU1"] rtol = 1e-10
    end

    @testset "errors" begin
        @test_throws ErrorException add!(net, pcu1)                       # merged already
        a = Assembly("A")
        @test_throws ErrorException connect!(FlowNetwork(), a => Duct("D"))  # no :outlet exposed
        b = Assembly("B"); d = Duct("D"); add!(b, d); expose!(b, :inlet, d)
        set_state!(b.net, d; Pt = 1e5, Tt = 300.0, W = 1.0, fluid = IdealGasFluid(M_molar = 83.8))
        n2 = FlowNetwork(); add!(n2, b)                                   # seed adopted
        c = Assembly("C"); d2 = Duct("D2"); add!(c, d2)
        set_state!(c.net, d2; Pt = 1e5, Tt = 300.0, W = 1.0, fluid = IdealGasFluid(M_molar = 83.8))
        @test_throws ErrorException add!(n2, c)                           # second seed
    end
end
