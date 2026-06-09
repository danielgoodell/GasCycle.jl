using Test
using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
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
    fpt_path = joinpath(@__DIR__, "..", "HeXe84.fpt")
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
    end
end
