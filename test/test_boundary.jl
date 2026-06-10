using Test
using GasCycle

@testset "Boundary-stream seeding (set_boundary!)" begin
    gas = IdealGasFluid(M_molar = 83.8)
    oil = IdealGasFluid(M_molar = 6.206)   # cp ≈ 3349 J/(kg·K) surrogate liquid

    # Open gas path through a heat-rejection HX whose cold side is a fixed
    # boundary stream: Duct → Sink(hot side); oil enters Sink.cold_inlet.
    net  = FlowNetwork()
    duct = Duct("Inlet"; dPqP = 0.0)
    sink = HeatExchanger("Sink"; ε = 0.946, dPqP_hot = 0.0173, dPqP_cold = 0.005)
    add!(net, duct, sink)
    # Gas enters the HOT side (heat rejection) — explicit port connection,
    # since the serial-chain default aliases :inlet to the HX cold side.
    connect_port!(net, duct, :outlet, sink, :hot_inlet)
    set_state!(net, duct; Pt = 165e3, Tt = 437.0, W = 0.59, fluid = gas)

    # Without the boundary, the HX can never become ready
    @test_throws Exception one_pass!(net)

    set_boundary!(net, sink, :cold_inlet;
                  Pt = 622e3, Tt = 293.0, W = 0.0635, fluid = oil)
    sol = solve!(net)
    @test sol.status == :success

    # Gas is C_min here: gas outlet follows directly from the ε definition
    Th, Tc = 437.0, 293.0
    C_gas = 0.59 * GasCycle.cp(gas, Th, 165e3)
    C_oil = 0.0635 * GasCycle.cp(oil, Tc, 622e3)
    @test C_gas < C_oil
    @test sink.hot_outlet[].Tt ≈ Th - 0.946 * (Th - Tc) rtol = 1e-10

    # Energy balance across both streams
    Q_gas = C_gas * (Th - sink.hot_outlet[].Tt)
    Q_oil = C_oil * (sink.cold_outlet[].Tt - Tc)
    @test Q_gas ≈ Q_oil rtol = 1e-10

    # Stations include the boundary inlet and the dangling oil outlet
    labels = first.(stations(sol))
    @test "Sink.cold_in"  in labels
    @test "Sink.cold_out" in labels
    @test labels[1:3] == ["Inlet.in", "Inlet.out", "Sink.hot_out"]

    # Replacing a boundary state takes effect (no duplicates)
    set_boundary!(net, sink, :cold_inlet;
                  Pt = 622e3, Tt = 300.0, W = 0.0635, fluid = oil)
    @test length(net.boundaries) == 1
    solve!(net)
    @test sink.hot_outlet[].Tt ≈ Th - 0.946 * (Th - 300.0) rtol = 1e-10

    # summary() renders with a boundary stream present
    @test occursin("Sink.cold_in", sprint(summary, sol))
end
