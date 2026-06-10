using Test
using GasCycle

@testset "summary(sol) and stations" begin
    fluid = IdealGasFluid(M_molar = 83.8)

    net    = FlowNetwork()
    comp   = Compressor("Comp"; PR = 1.9, η_poly = 0.80)
    bsplit = Splitter("BldSplit"; fracs = [0.98, 0.02])
    recup  = HeatExchanger("Recup"; ε = 0.95, dPqP_hot = 0.022, dPqP_cold = 0.011)
    bcool  = HeatSource("BldCool"; TtExit = 310.0, mode = :fixed_TtExit)
    heater = HeatSource("Heater"; TtExit = 1144.0, dPqP = 0.027)
    bmix   = Mixer("BldMix")
    turb   = Turbine("Turb"; mode = :pressure_closure, P_exit = 170.0e3, η_poly = 0.87)
    shaft  = Shaft("Shaft"; N = 36000.0)

    add!(net, comp, bsplit, recup, bcool, heater, bmix, turb)
    connect!(net, comp => bsplit => recup => heater => bmix => turb => comp)
    add_shaft!(net, shaft; drives = comp, driven_by = turb)
    add_hx_pair!(net, recup; hot = turb)
    connect_port!(net, bsplit, :bleed_outlet, bcool, :inlet)
    connect_port!(net, bcool, :outlet, bmix, :bleed_inlet)
    set_state!(net, comp; Pt = 163.4e3, Tt = 300.0, W = 0.6, fluid = fluid)

    sol = solve!(net)
    @test sol.status == :success

    # Station list comes out in physical flow order, main loop then branches.
    sts    = stations(sol)
    labels = first.(sts)
    @test labels == ["Comp.in", "Comp.out", "BldSplit.out", "Recup.cold_out",
                     "Heater.out", "BldMix.out", "Turb.out", "Recup.hot_out",
                     "BldSplit.bleed_out", "BldCool.out"]
    @test allunique(labels)

    # Station states are the solved ones, not stale copies.
    d = Dict(sts)
    @test d["Comp.in"].Tt ≈ 300.0
    @test d["Heater.out"].Tt ≈ 1144.0
    @test d["Turb.out"].Pt ≈ 170.0e3
    @test d["BldSplit.bleed_out"].W ≈ 0.02 * 0.6
    @test d["BldMix.out"].W ≈ 0.6

    # summary(io, sol) renders without error and contains the key sections.
    txt = sprint(summary, sol)
    @test occursin("GasCycle solution — success", txt)
    @test occursin("Station", txt)
    @test occursin("Recup.hot_out", txt)
    @test occursin("HeatExchanger", txt)
    @test occursin("N=36000 rpm", txt)
    @test occursin("Net shaft power", txt)
    @test occursin("Cycle efficiency", txt)

    # Open chain (no loop, no branches) also works.
    net2  = FlowNetwork()
    duct  = Duct("Inlet"; dPqP = 0.01)
    comp2 = Compressor("Comp"; PR = 2.0, η_poly = 0.85)
    add!(net2, duct, comp2)
    connect!(net2, duct => comp2)
    set_state!(net2, duct; Pt = 100e3, Tt = 300.0, W = 1.0, fluid = fluid)
    sol2 = solve!(net2)
    @test first.(stations(sol2)) == ["Inlet.in", "Inlet.out", "Comp.out"]
    @test occursin("Net shaft power", sprint(summary, sol2))
end
