using Test
using GasCycle
import GasCycle: cp   # explicit import to avoid ambiguity with Base.Filesystem.cp

@testset "IdealGasFluid — monatomic HeXe" begin
    x_He = 0.47
    fl = HeXeIdealGas(x_He)

    # Compute expected Cp using same constants as IdealGasFluid
    M_He, M_Xe = 4.002602, 131.293
    M_mix = x_He * M_He + (1.0 - x_He) * M_Xe
    R_universal = 8314.46261815324
    Cp_expected = 2.5 * R_universal / M_mix

    @test cp(fl, 500.0, 1e5)   ≈ Cp_expected  rtol=1e-6
    @test cp(fl, 1000.0, 5e5)  ≈ Cp_expected  rtol=1e-6   # constant for monatomic

    # h = Cp * T  (referenced to 0 K)
    @test enthalpy(fl, 500.0, 1e5) ≈ Cp_expected * 500.0  rtol=1e-6

    # Entropy increases with T, decreases with P
    s1 = entropy(fl, 500.0, 1e5)
    s2 = entropy(fl, 800.0, 1e5)
    s3 = entropy(fl, 500.0, 2e5)
    @test s2 > s1
    @test s3 < s1

    # Density follows ideal gas law ρ = P / (R * T)
    R_specific = R_universal / M_mix
    ρ = density(fl, 300.0, 1e5)
    @test ρ ≈ 1e5 / (R_specific * 300.0)  rtol=1e-6

    # Gamma = 5/3 for monatomic ideal gas
    @test gamma(fl, 500.0, 1e5) ≈ 5.0/3.0  rtol=1e-6

    # T_from_h is the exact inverse of enthalpy for ideal gas (analytic)
    T_test = 750.0
    h_test = enthalpy(fl, T_test, 2e5)
    @test T_from_h(fl, h_test, 2e5) ≈ T_test  rtol=1e-8

    # T_from_s is the exact inverse of entropy for ideal gas (analytic)
    s_test = entropy(fl, T_test, 2e5)
    @test T_from_s(fl, s_test, 2e5) ≈ T_test  rtol=1e-6
end

@testset "FPTFluid — HeXe84.fpt" begin
    fpt_path = joinpath(@__DIR__, "..", "data", "HeXe84.fpt")
    if !isfile(fpt_path)
        @warn "HeXe84.fpt not found; skipping FPT tests"
    else
        fl = FPTFluid(fpt_path)

        # Cp should be ~247 J/(kg·K) from the file (Cp≈0.05904 BTU/(lbm·R))
        Cp_expected = 0.05904 * 4186.8   # ≈ 247.2 J/(kg·K)

        # Test at a point well inside the grid
        T_test = 900.0   # K  (= 1620 R, well inside 450–2520 R range)
        P_test = 101325.0 * 2   # Pa (= ~2 atm ≈ 29.4 psia, inside 0.7–435 psia range)

        @test GasCycle.cp(fl, T_test, P_test) ≈ Cp_expected  rtol=0.05

        # h should increase with T
        h1 = enthalpy(fl, 700.0, P_test)
        h2 = enthalpy(fl, 900.0, P_test)
        h3 = enthalpy(fl, 1100.0, P_test)
        @test h1 < h2 < h3

        # s should increase with T at fixed P
        s1 = entropy(fl, 700.0, P_test)
        s2 = entropy(fl, 900.0, P_test)
        @test s2 > s1

        # s should decrease with P at fixed T
        s_lo = entropy(fl, 900.0, P_test)
        s_hi = entropy(fl, 900.0, P_test * 2)
        @test s_hi < s_lo

        # T_from_h round-trip
        h_val = enthalpy(fl, T_test, P_test)
        @test T_from_h(fl, h_val, P_test; T_guess=T_test) ≈ T_test  rtol=1e-3

        # T_from_s round-trip
        s_val = entropy(fl, T_test, P_test)
        @test T_from_s(fl, s_val, P_test; T_guess=T_test) ≈ T_test  rtol=1e-3

        # Out-of-table lookups fail fast by default with a useful bounds error.
        err = try
            enthalpy(fl, fl.Tt_max + 10.0, P_test)
        catch e
            e
        end
        @test err isa DomainError
        msg = sprint(showerror, err)
        @test occursin("outside table bounds", msg)
        @test occursin("enthalpy", msg)
        @test occursin("bounds=:warn", msg)
        @test occursin("bounds=:clamp", msg)

        # Compatibility modes are explicit: warn+clamp or silent clamp.
        fl_warn = FPTFluid(fpt_path; bounds=:warn)
        @test_logs (:warn, r"cp requested outside table bounds") begin
            GasCycle.cp(fl_warn, fl_warn.Tt_max + 10.0, P_test)
        end

        fl_clamp = FPTFluid(fpt_path; bounds=:clamp)
        @test GasCycle.cp(fl_clamp, fl_clamp.Tt_max + 10.0, P_test) ≈
              GasCycle.cp(fl_clamp, fl_clamp.Tt_max, P_test)

        @test_throws ErrorException FPTFluid(fpt_path; bounds=:extrapolate)

        # ── Entropy pressure-interpolation (validation/PLAN.md rung 0) ──────
        # s ≈ f(T) − R·ln P, so linear-in-P interpolation overshoots ∂s/∂lnP
        # mid-cell on the coarse Pt grid (~14 % at the BRU compressor inlet,
        # ⇒ +14 °R on the compressor outlet).  The default :log_pressure mode
        # detrends the log term and must recover ∂s/∂lnP ≈ −R everywhere.
        R_HeXe84 = 8314.46 / 83.8
        @test fl.s_interp == :log_pressure
        @test fl.R_s ≈ R_HeXe84 rtol = 0.005   # fitted from the table itself

        # Mid-cell pressure derivative (163 kPa sits between the 101.6 and
        # 198.2 kPa nodes — the worst case that produced the BRU offset)
        dsdlnP = (entropy(fl, 300.0, 180e3) - entropy(fl, 300.0, 150e3)) /
                 log(180 / 150)
        @test dsdlnP ≈ -R_HeXe84 rtol = 0.01

        fl_lin = FPTFluid(fpt_path; s_interp = :linear)
        dsdlnP_lin = (entropy(fl_lin, 300.0, 180e3) - entropy(fl_lin, 300.0, 150e3)) /
                     log(180 / 150)
        @test abs(dsdlnP_lin) > 1.1 * R_HeXe84   # legacy artifact, kept for NPSS-compat

        # Table-node values are mode-independent (detrend is exact at nodes)
        P_node = fl.itp_s.knots[1][3]
        @test entropy(fl, 900.0, P_node) ≈ entropy(fl_lin, 900.0, P_node) rtol = 1e-12

        # Isentropic step at the BRU compressor state now implies a physical
        # exponent consistent with the table's own γ (was γ_eff ≈ 1.74)
        s_bru = entropy(fl, 300.0, 163.4e3)
        T2s   = T_from_s(fl, s_bru, 310.5e3; T_guess = 380.0)
        γ_eff = 1 / (1 - log(T2s / 300.0) / log(310.5 / 163.4))
        @test γ_eff ≈ gamma(fl, 300.0, 163.4e3) rtol = 0.005

        @test_throws ErrorException FPTFluid(fpt_path; s_interp = :cubic)
    end
end
