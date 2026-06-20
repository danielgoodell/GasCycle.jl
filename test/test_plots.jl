using Test
using GasCycle

# RecipesBase is a GasCycle dependency; reach it through the package so the
# test target doesn't need its own entry.
const RB = GasCycle.RecipesBase

@testset "Plot recipes" begin
    fluid = IdealGasFluid(M_molar = 83.8)

    net    = FlowNetwork()
    comp   = Compressor("Comp"; PR = 1.9, η_poly = 0.80)
    recup  = HeatExchanger("Recup"; ε = 0.95, dPqP_hot = 0.022, dPqP_cold = 0.011)
    heater = HeatSource("Heater"; TtExit = 1144.0, dPqP = 0.027)
    turb   = Turbine("Turb"; mode = :pressure_closure, P_exit = 170.0e3, η_poly = 0.87)
    shaft  = Shaft("Shaft"; N = 36000.0)

    add!(net, comp, recup, heater, turb)
    connect!(net, comp => recup => heater => turb => comp)
    add_shaft!(net, shaft; drives = comp, driven_by = turb)
    add_hx_pair!(net, recup; hot = turb)
    set_state!(net, comp; Pt = 163.4e3, Tt = 300.0, W = 0.6, fluid = fluid)
    sol = solve!(net)
    @test sol.status == :success

    # ── T-s diagram ──────────────────────────────────────────────────────────
    series = RB.apply_recipe(Dict{Symbol,Any}(), GasCycle.TsDiagram((sol,)))
    # main path + dashed closure (closed loop has a back-edge) + labels
    @test length(series) == 3

    s_main, T_main = series[1].args
    @test length(s_main) == length(stations(sol; branches = false))
    @test T_main[1] ≈ 300.0          # seed station
    @test maximum(T_main) ≈ 1144.0   # TIT
    @test issorted([T_main[1], T_main[2]])  # compression heats the gas

    closure = series[2]
    @test closure.plotattributes[:linestyle] == :dash
    @test closure.args[1][2] == s_main[1] && closure.args[2][2] == T_main[1]

    @test_throws Exception RB.apply_recipe(Dict{Symbol,Any}(),
                                           GasCycle.TsDiagram((42,)))

    # ── Performance map plot ─────────────────────────────────────────────────
    # Real R-line compressor map, scaled to a design point.
    cm0  = compressor_map(joinpath(@__DIR__, "..", "data", "compressor_argon.map"))
    pmap = scale_map(cm0; Nc_des = 36000.0, Wc_des = 0.581, PR_des = 1.9, eta_des = 0.795)
    nlines = length(pmap.flow.speeds[1])              # one speed line per NcMap node

    series = RB.apply_recipe(Dict{Symbol,Any}(), GasCycle.MapPlot((pmap,)))
    @test length(series) == nlines
    @test length(series[1].args[1]) == length(series[1].args[2])  # (Wc, PR) pair

    # Turbomachine with map and solved state adds an operating point
    comp_od = Compressor("CompOD"; map = pmap, mode = :off_design)
    comp_od.N_shaft = 36000.0
    compute!(comp_od, Port(FluidState(101325.0, 288.15, 0.5, fluid)))
    series = RB.apply_recipe(Dict{Symbol,Any}(), GasCycle.MapPlot((comp_od,)))
    @test length(series) == nlines + 1
    op = series[end]
    @test op.plotattributes[:seriestype] == :scatter
    @test length(op.args[1]) == 1                     # single operating point

    @test_throws Exception RB.apply_recipe(Dict{Symbol,Any}(),
                                           GasCycle.MapPlot((Duct("D"),)))
end
