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
