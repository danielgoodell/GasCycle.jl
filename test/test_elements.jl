using Test
using GasCycle
import GasCycle: cp

# Shared ideal-gas fluid for all element tests
const FL = HeXeIdealGas(0.47)
const T0 = 400.0    # K
const P0 = 500e3    # Pa
const W0 = 10.0     # kg/s
const S0 = FluidState(P0, T0, W0, FL)

make_port() = Port(S0)

@testset "Compressor" begin
    PR   = 2.5
    η    = 0.88
    comp = Compressor("C"; PR=PR, η_poly=η)
    out  = compute!(comp, make_port())

    # Pressure must increase by PR
    @test out[].Pt ≈ P0 * PR  rtol=1e-6

    # Exit temperature must be > isentropic exit (efficiency < 1)
    γ    = gamma(FL, T0, P0)
    T_is = T0 * PR^((γ-1)/γ)
    @test out[].Tt > T_is
    @test out[].Tt > T0

    # Specific work = Cp * (Tt_out - Tt_in) for ideal gas
    Cp_val = cp(FL, T0, P0)
    W_spec = specific_work(comp)
    @test W_spec ≈ Cp_val * (out[].Tt - T0)  rtol=1e-4

    # Design mode: no residuals, PR is a fixed parameter (not a solver unknown)
    @test n_residuals(comp) == 0
    @test length(indep_vars(comp)) == 0
end

@testset "Turbine" begin
    PR   = 2.3
    η    = 0.90
    turb = Turbine("T"; PR=PR, η_poly=η)

    # Start with hot inlet (turbine inlet temperature)
    hot_state = FluidState(P0 * 2.5, 1100.0, W0, FL)
    out = compute!(turb, Port(hot_state))

    # Pressure must drop by PR
    @test out[].Pt ≈ hot_state.Pt / PR  rtol=1e-6

    # Temperature must drop (turbine extracts work)
    @test out[].Tt < hot_state.Tt

    # Isentropic exit is colder than actual (efficiency < 1)
    γ    = gamma(FL, hot_state.Tt, hot_state.Pt)
    T_is = hot_state.Tt / PR^((γ-1)/γ)
    @test out[].Tt > T_is

    @test n_residuals(turb) == 0
end

@testset "Duct" begin
    dP = 0.03
    duct = Duct("D"; dPqP=dP)
    out  = compute!(duct, make_port())

    @test out[].Pt ≈ P0 * (1 - dP)  rtol=1e-6
    @test out[].Tt ≈ T0              rtol=1e-6   # isenthalpic
    @test out[].W  ≈ W0              rtol=1e-6
    @test n_residuals(duct) == 0
end

@testset "HeatSource :fixed_TtExit" begin
    TtExit = 1100.0
    hs = HeatSource("R"; TtExit=TtExit, dPqP=0.02, mode=:fixed_TtExit)
    out = compute!(hs, make_port())

    @test out[].Tt ≈ TtExit           rtol=1e-4
    @test out[].Pt ≈ P0 * 0.98       rtol=1e-6
    @test hs.Q > 0.0                  # heat was added
end

@testset "HeatSource :fixed_Q" begin
    Q_in = 500e3   # 500 kW
    hs = HeatSource("R"; Q=Q_in, mode=:fixed_Q)
    out = compute!(hs, make_port())

    Cp_val = cp(FL, T0, P0)
    ΔT_expected = Q_in / (W0 * Cp_val)
    @test out[].Tt ≈ T0 + ΔT_expected  rtol=1e-3
end

@testset "HeatExchanger effectiveness" begin
    ε = 0.90
    hx = HeatExchanger("HX"; ε=ε, dPqP_hot=0.01, dPqP_cold=0.01)

    # Hot stream: high T, same fluid and mass flow
    hot_state  = FluidState(P0 * 2.5, 900.0, W0, FL)
    cold_state = FluidState(P0,       500.0, W0, FL)

    hx.hot_inlet  = Port(hot_state)
    hx.cold_inlet = Port(cold_state)

    hot_out, cold_out = compute_hx!(hx)

    # Cold outlet must be warmer, hot outlet must be cooler
    @test cold_out[].Tt > cold_state.Tt
    @test hot_out[].Tt  < hot_state.Tt

    # Energy balance: Q_transferred from hot = Q_received by cold
    Cp_val = cp(FL, 700.0, P0)
    Q_hot  = hot_state.W  * Cp_val * (hot_state.Tt  - hot_out[].Tt)
    Q_cold = cold_state.W * Cp_val * (cold_out[].Tt - cold_state.Tt)
    @test Q_hot ≈ Q_cold  rtol=1e-3

    # Maximum possible Q (Cmin side limits)
    Q_max = min(hot_state.W, cold_state.W) * Cp_val * (hot_state.Tt - cold_state.Tt)
    @test Q_transferred(hx) ≈ ε * Q_max  rtol=1e-3
end

@testset "Splitter" begin
    sp = Splitter("Sp"; fracs=[0.98, 0.02])
    out_main = compute!(sp, make_port())

    # Fractions must sum to 1
    @test_throws ErrorException Splitter("Bad"; fracs=[0.5, 0.6])

    # Main outlet carries the right mass flow
    @test out_main[].W ≈ W0 * 0.98   rtol=1e-6
    @test sp.outlets[2][].W ≈ W0 * 0.02  rtol=1e-6

    # Same Tt and Pt on all outlets
    @test out_main[].Tt ≈ T0  rtol=1e-6
    @test out_main[].Pt ≈ P0  rtol=1e-6
    @test sp.outlets[2][].Tt ≈ T0  rtol=1e-6

    @test n_residuals(sp) == 0
    @test length(indep_vars(sp)) == 0
end

@testset "Mixer" begin
    # Two streams: different T, same fluid and P
    s1 = FluidState(P0, 800.0, W0 * 0.8, FL)
    s2 = FluidState(P0, 400.0, W0 * 0.2, FL)
    mx = Mixer("Mx")
    mx.inlets[1] = Port(s1)
    mx.inlets[2] = Port(s2)
    compute!(mx)

    # Mass balance
    @test mx.outlet[].W ≈ W0  rtol=1e-6

    # Energy balance: h_out = (W1*h1 + W2*h2) / W_total
    h1   = GasCycle.enthalpy(FL, s1.Tt, s1.Pt)
    h2   = GasCycle.enthalpy(FL, s2.Tt, s2.Pt)
    h_ex = (s1.W * h1 + s2.W * h2) / (s1.W + s2.W)
    h_ac = GasCycle.enthalpy(FL, mx.outlet[].Tt, mx.outlet[].Pt)
    @test h_ac ≈ h_ex  rtol=1e-4

    # Mixed Tt must be between the two inlet temperatures
    @test mx.outlet[].Tt > s2.Tt
    @test mx.outlet[].Tt < s1.Tt

    @test n_residuals(mx) == 0
    @test length(indep_vars(mx)) == 0
end

@testset "Shaft power balance" begin
    comp = Compressor("C"; PR=2.5, η_poly=0.88)
    turb = Turbine("T";    PR=2.1, η_poly=0.90)
    sh   = Shaft("S";      N=15000.0)

    hot_state = FluidState(P0 * 2.5, 1100.0, W0, FL)
    compute!(comp, make_port())
    compute!(turb, Port(hot_state))

    link!(sh, [turb], [comp])

    # Design mode: shaft balance is an output (not a solver constraint)
    @test n_residuals(sh) == 0
    @test length(indep_vars(sh)) == 0

    # power_balance returns a scalar (not necessarily zero for arbitrary PRs)
    pb = power_balance(sh)
    @test pb isa Float64
end

@testset "ForwardDiff: net power gradient" begin
    using ForwardDiff

    fluid = HeXeIdealGas(0.47)

    function cycle_power(x::AbstractVector)
        PR_comp, ε_recup = x[1], x[2]
        T0, P0, W = 400.0, 500e3, 10.0

        net    = FlowNetwork()
        comp   = Compressor("C"; PR=PR_comp, η_poly=0.88)
        recup  = HeatExchanger("R"; ε=ε_recup, dPqP_hot=0.01, dPqP_cold=0.01)
        heater = HeatSource("H"; TtExit=1100.0, dPqP=0.02)
        P_exit = P0 / ((1-0.01)*(1-0.01))
        turb   = Turbine("T"; mode=:pressure_closure, P_exit=P_exit, η_poly=0.90)

        add!(net, comp, recup, heater, turb)
        connect!(net, comp => recup => heater => turb => comp)
        add_hx_pair!(net, recup; hot=turb)
        set_state!(net, comp; Pt=P0, Tt=T0, W=W, fluid=fluid)
        sol = solve!(net; maxiter=200)
        net_power(sol) / 1000
    end

    x0 = [2.5, 0.90]
    ∇W = ForwardDiff.gradient(cycle_power, x0)

    # Gradient must be non-NaN and finite
    @test all(isfinite.(∇W))

    # ∂W/∂PR_comp > 0: higher pressure ratio → more turbine work (physics check)
    @test ∇W[1] > 0

    # Verify against forward-difference (tolerance 0.01%)
    h = 1e-4
    fd1 = (cycle_power([x0[1]+h*x0[1], x0[2]]) - cycle_power([x0[1]-h*x0[1], x0[2]])) / (2h*x0[1])
    @test ∇W[1] ≈ fd1  rtol=1e-4
end
