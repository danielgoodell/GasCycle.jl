using Test
using GasCycle
using ForwardDiff
import GasCycle: cp

"""
NobleGasMixture validation (roadmap item 4).

Oracles from Tournier, El-Genk & Gallo, AIAA 2006-4154 (reference/):
spot values at 2 MPa / 400 K, ideal-gas limits, virial-coefficient checks,
plus exact thermodynamic-consistency identities the correlations must obey
once assembled into an EOS.
"""

@testset "NobleGasMixture — virial coefficients" begin
    # B_He(300 K) ≈ +11.8 cm³/mol (Sengers data, paper Eq. 10)
    @test GasCycle._B_pure(HELIUM, 300.0) * 1e6 ≈ 11.8 atol = 0.2
    # Fig. 4: B/V* ≈ −0.33 at θ = 1 (any heavy gas, corresponding states)
    @test GasCycle._B_pure(XENON, XENON.Tcr) / (8.31441 * XENON.Tcr / XENON.Pcr) ≈
          -0.334 atol = 0.01
    # Boyle point at θ ≈ 2.7 (textbook LJ value; sanity on Eq. 9 exponents)
    @test GasCycle._ΨB2(2.6)[1] < 0 < GasCycle._ΨB2(2.9)[1]
    # He-Xe B12 from the LJ series lands in the paper's Fig. 6 data band
    @test 12.0 < GasCycle._B12(HELIUM, XENON, 400.0) * 1e6 < 25.0
    # LJ series reproduces Hirschfelder table values
    @test GasCycle._Bstar_LJ2(1.0)[1] ≈ -2.538 atol = 0.001
    @test GasCycle._Bstar_LJ2(10.0)[1] ≈ 0.46088 atol = 0.0005
end

@testset "NobleGasMixture — El-Genk spot values (2 MPa, 400 K)" begin
    for (gas, cp_ref, γ_ref) in ((ARGON, 21.257, 1.697),
                                 (KRYPTON, 21.72, 1.724),
                                 (XENON, 23.00, 1.786))
        fp = NobleGasFluid(gas)
        @test cp(fp, 400.0, 2e6) * fp.M ≈ cp_ref rtol = 5e-4
        @test gamma(fp, 400.0, 2e6) ≈ γ_ref rtol = 5e-4
    end
    # compressibility: Z(Xe) = 0.957, Z(Kr) = 0.987 at 2 MPa / 400 K
    Zof(g) = 2e6 * g.M / (density(NobleGasFluid(g), 400.0, 2e6) * 8.31441 * 400.0)
    @test Zof(XENON) ≈ 0.957 atol = 0.001
    @test Zof(KRYPTON) ≈ 0.987 atol = 0.001

    # ideal-gas limit at 0.1 MPa
    fp = NobleGasFluid(ARGON)
    @test cp(fp, 400.0, 1e5) * fp.M ≈ 2.5 * 8.31441 rtol = 2e-3
    @test gamma(fp, 400.0, 1e5) ≈ 5 / 3 rtol = 1.5e-3
end

@testset "NobleGasMixture — thermodynamic consistency" begin
    # All three cp paths must agree exactly: Eq. 13 analytic, ∂h/∂T (AD),
    # and T·∂s/∂T (AD).  Covers all three composition branches.
    for fp in (HeXe(83.8), NobleGasFluid(XENON), NobleGasMixture(ARGON, KRYPTON, 0.5))
        for (T, P) in ((300.0, 163.4e3), (400.0, 2e6), (1144.0, 297.9e3))
            cp_an = cp(fp, T, P)
            @test ForwardDiff.derivative(t -> enthalpy(fp, t, P), T) ≈ cp_an rtol = 1e-10
            @test T * ForwardDiff.derivative(t -> entropy(fp, t, P), T) ≈ cp_an rtol = 1e-10
        end
    end

    hexe = HeXe(83.8)
    # Maxwell relation: (∂h/∂P)_T = v − T(∂v/∂T)_P
    T, P = 400.0, 1e6
    dhdP = ForwardDiff.derivative(p -> enthalpy(hexe, T, p), P)
    v(t) = 1 / density(hexe, t, P)
    @test dhdP ≈ v(T) - T * ForwardDiff.derivative(v, T) rtol = 1e-8

    # dilute limit: ∂s/∂lnP → −R
    dsdlnP = ForwardDiff.derivative(p -> entropy(hexe, 400.0, p), 1e5) * 1e5
    @test dsdlnP ≈ -8.31441 / hexe.M rtol = 3e-3

    # inversion round-trips (Newton, exact derivatives)
    @test T_from_h(hexe, enthalpy(hexe, 700.0, P), P) ≈ 700.0 atol = 1e-8
    @test T_from_s(hexe, entropy(hexe, 700.0, P), P) ≈ 700.0 atol = 1e-6
    @test T_from_h(hexe, enthalpy(hexe, 700.0, P), P; T_guess = 400.0) ≈ 700.0 atol = 1e-8
    # T_from_s honors the caller's guess: same root from a near-exact guess,
    # a far-off guess (isentropic full-PR call pattern), and one bad enough
    # to trip the half-T step cap
    s700 = entropy(hexe, 700.0, P)
    @test T_from_s(hexe, s700, P; T_guess = 701.0) ≈ 700.0 atol = 1e-6
    @test T_from_s(hexe, s700, P; T_guess = 400.0) ≈ 700.0 atol = 1e-6
    @test T_from_s(hexe, s700, P; T_guess = 4000.0) ≈ 700.0 atol = 1e-6
end

@testset "NobleGasMixture — composition and AD through x₁" begin
    hexe = HeXe(83.8)
    @test hexe.M ≈ 0.0838 rtol = 1e-12
    @test hexe.x1 ≈ 0.3731 atol = 1e-3          # mole fraction He
    # paper: M = 40 g/mol ⇒ 72 mole% He
    @test HeXe(40.0).x1 ≈ 0.72 atol = 0.005

    # exact mixture-ratio sensitivity (ideal limit has a closed form)
    g = ForwardDiff.derivative(x -> cp(NobleGasMixture(HELIUM, XENON, x), 400.0, 1e5), 0.72)
    M72 = 0.72 * HELIUM.M + 0.28 * XENON.M
    @test g ≈ -2.5 * 8.31441 * (HELIUM.M - XENON.M) / M72^2 rtol = 1e-3

    @test_throws Exception NobleGasMixture(HELIUM, XENON, 1.2)
end

@testset "NobleGasMixture — vs HeXe84.fpt at cycle conditions" begin
    fpt_path = joinpath(@__DIR__, "..", "HeXe84.fpt")
    if !isfile(fpt_path)
        @warn "HeXe84.fpt not found; skipping comparison"
    else
        hexe = HeXe(83.8)
        fpt  = FPTFluid(fpt_path)
        for (T, P) in ((300.0, 163.4e3), (417.0, 310.5e3), (1144.4, 297.9e3))
            @test density(hexe, T, P) ≈ density(fpt, T, P) rtol = 0.01
            @test cp(hexe, T, P) ≈ cp(fpt, T, P) rtol = 0.02
            @test gamma(hexe, T, P) ≈ gamma(fpt, T, P) rtol = 0.01
        end
        # Known table artifact (see validation/RESULTS.md): near the Xe
        # corner the FPT Cp drops BELOW ideal, which is inconsistent with
        # its own (correct) density.  The virial model rises above ideal,
        # as thermodynamic consistency demands.  Densities still agree.
        @test density(hexe, 260.0, 1.5e6) ≈ density(fpt, 260.0, 1.5e6) rtol = 0.005
        cp_ideal = 2.5 * 8.31441 / 0.0838
        @test cp(hexe, 260.0, 1.5e6) > cp_ideal          # virial: above ideal
        @test cp(fpt, 260.0, 1.5e6) < cp_ideal           # table: artifact
    end
end

@testset "NobleGasMixture — drop-in cycle solve" begin
    # Recuperated closed-loop design solve with the new backend; the
    # ideal-gas backend brackets the answer (real-gas effects are small at
    # BRU conditions).
    function bru_net_power(fluid)
        comp  = Compressor("Comp"; PR = 1.9, η_poly = 0.80, η_type = :isentropic)
        recup = HeatExchanger("Recup"; ε = 0.95, dPqP_hot = 0.022, dPqP_cold = 0.011)
        heat  = HeatSource("Heater"; TtExit = 1144.0, dPqP = 0.027)
        turb  = Turbine("Turb"; mode = :pressure_closure,
                        P_exit = 163.4e3 / (1 - 0.022) / (1 - 0.0173),
                        η_poly = 0.87, η_type = :isentropic)
        shaft = Shaft("Shaft"; N = 36_000.0)
        net = FlowNetwork()
        add!(net, comp, recup, heat, turb)
        connect!(net, comp => recup => heat => turb)
        add_shaft!(net, shaft; drives = comp, driven_by = turb)
        add_hx_pair!(net, recup; hot = turb)
        set_state!(net, comp; Pt = 163.4e3, Tt = 300.0, W = 0.6, fluid = fluid)
        sol = solve!(net)
        @test sol.status == :success
        net_power(sol)
    end
    P_virial = bru_net_power(HeXe(83.8))
    P_ideal  = bru_net_power(IdealGasFluid(M_molar = 83.8))
    @test P_virial ≈ P_ideal rtol = 0.01
    @test P_virial != P_ideal   # real-gas terms are present, just small
end

@testset "NobleGasMixture — transport: pure-gas limits" begin
    # Eq. 7 makes Pr = cp·μ/λ = 2/3 exactly for a pure dilute monatomic gas
    for g in (HELIUM, NEON, ARGON, KRYPTON, XENON)
        @test prandtl(NobleGasFluid(g), 400.0, 1e4) ≈ 2 / 3 atol = 1e-3
    end
    # dilute viscosities vs handbook values
    @test viscosity(NobleGasFluid(HELIUM), 300.0, 1e4) ≈ 19.9e-6 rtol = 0.02
    @test viscosity(NobleGasFluid(XENON), 300.0, 1e4) ≈ 23.2e-6 rtol = 0.02

    # Eq. 23b/33b reproduce the Table 1 critical-point transport column:
    # μ* vs μcr = Δμcr/(1−1/2.3), λ* vs λcr·(1 + printed λ*cr deviation)
    Vst(g) = 8.31441 * g.Tcr / g.Pcr
    @test GasCycle._mustar(XENON.M, XENON.Tcr, Vst(XENON)) ≈ 52.25e-6 rtol = 0.015
    for (g, dev) in ((HELIUM, 0.065), (NEON, 0.0047), (ARGON, -0.0050),
                     (KRYPTON, 0.0099), (XENON, -0.0028))
        @test GasCycle._lamstar(g.M, g.Tcr, Vst(g)) ≈ g.λcr * (1 + dev) rtol = 1e-3
    end

    # El-Genk §V pressure effects at 300 K / 2 MPa (vs the dilute limit)
    xe = NobleGasFluid(XENON)
    @test viscosity(xe, 300.0, 2e6) / viscosity(xe, 300.0, 1e3) - 1 ≈ 0.043 atol = 0.005
    @test conductivity(xe, 300.0, 2e6) / conductivity(xe, 300.0, 1e3) - 1 ≈ 0.123 atol = 0.02
    m40 = HeXe(40.0)
    @test viscosity(m40, 300.0, 2e6) / viscosity(m40, 300.0, 1e3) - 1 ≈ 0.005 atol = 0.002
    @test conductivity(m40, 300.0, 2e6) / conductivity(m40, 300.0, 1e3) - 1 ≈ 0.008 atol = 0.002

    # composition endpoints collapse to the pure-gas path
    for (T, P) in ((400.0, 2e6),)
        @test viscosity(NobleGasMixture(HELIUM, XENON, 1.0), T, P) ==
              viscosity(NobleGasFluid(HELIUM), T, P)
        @test conductivity(NobleGasMixture(HELIUM, XENON, 0.0), T, P) ==
              conductivity(NobleGasFluid(XENON), T, P)
    end
end

@testset "NobleGasMixture — transport vs Johnson NASA/CR-2006-214394" begin
    # He-Xe oracles: Johnson Tables 4-6 (Hirschfelder LJ theory; μ first
    # order, k third order) and Table 7 (Taylor 1988 experiment).  El-Genk's
    # data-fitted correlations sit systematically above LJ theory with a
    # T-growing offset (μ to +5%, k to +9% at 1200 K) that CANCELS in Pr —
    # our Pr matches Taylor's measurements to ≤1.1% at all four molecular
    # weights, where Johnson's own method is high by 2.3-4.7%.
    P = 0.1e6   # dilute-limit comparison pressure
    johnson = (  # (M kg/kmol, T K, μ×10⁶ Pa·s, k(3rd order) W/m·K)
        (20.183, 400.0, 31.6340, 0.124434),
        (20.183, 800.0, 50.7860, 0.195756),
        (20.183, 1200.0, 66.4595, 0.254660),
        (39.94, 400.0, 33.1088, 0.080589),
        (39.94, 800.0, 54.5763, 0.127540),
        (39.94, 1200.0, 71.9178, 0.166432),
        (83.8, 400.0, 31.8624, 0.032369),
        (83.8, 800.0, 54.3585, 0.052095),
        (83.8, 1200.0, 72.2199, 0.068370),
    )
    for (M, T, μ_o, k_o) in johnson
        f = HeXe(M)
        @test viscosity(f, T, P) ≈ μ_o * 1e-6 rtol = 0.05
        @test conductivity(f, T, P) ≈ k_o rtol = 0.09
    end
    # Taylor 1988 experimental Prandtl numbers (Johnson Table 7)
    for (M, T, Pr_o) in ((14.5, 972.0, 0.301), (28.3, 982.0, 0.231),
                         (40.0, 941.0, 0.214), (83.8, 962.0, 0.251))
        @test prandtl(HeXe(M), T, P) ≈ Pr_o rtol = 0.015
    end
end

@testset "NobleGasMixture — transport AD" begin
    hexe = HeXe(83.8)
    T, P = 800.0, 1.5e6
    # T- and P-derivatives propagate and match central differences
    for f in (viscosity, conductivity, prandtl)
        dT = ForwardDiff.derivative(t -> f(hexe, t, P), T)
        @test dT ≈ (f(hexe, T + 0.5, P) - f(hexe, T - 0.5, P)) rtol = 1e-4
        dP = ForwardDiff.derivative(p -> f(hexe, T, p), P)
        @test dP ≈ (f(hexe, T, P + 500.0) - f(hexe, T, P - 500.0)) / 1e3 rtol = 1e-4
    end
    # derivative through the mixture ratio (Dual-valued x1 via HeXe(M))
    dμdM = ForwardDiff.derivative(m -> viscosity(HeXe(m), T, P), 83.8)
    @test dμdM ≈ (viscosity(HeXe(83.9), T, P) - viscosity(HeXe(83.7), T, P)) / 0.2 rtol = 1e-4
end

@testset "NobleGasMixture — h_from_s interface fallback" begin
    hexe = HeXe(83.8)
    s = entropy(hexe, 900.0, 2.0e5)
    h = h_from_s(hexe, s, 1.5e5)   # generic fallback: enthalpy ∘ T_from_s
    @test entropy(hexe, T_from_h(hexe, h, 1.5e5), 1.5e5) ≈ s rtol = 1e-8
end
