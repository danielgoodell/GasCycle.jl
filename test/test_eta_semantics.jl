using Test
using GasCycle

# Ideal monatomic gas closed forms (γ = 5/3, k = (γ−1)/γ = 0.4):
#   isentropic-η compressor:  T_out = T_in (1 + (PR^k − 1)/η)
#   isentropic-η turbine:     T_out = T_in (1 − η (1 − PR^−k))
#   polytropic-η compressor:  T_out = T_in PR^(k/η)
#   polytropic-η turbine:     T_out = T_in PR^(−kη)

@testset "Efficiency semantics (η_type)" begin
    fluid = IdealGasFluid(M_molar = 83.8)
    k     = 0.4
    T_in  = 300.0
    P_in  = 163.4e3
    inlet() = Port(FluidState(P_in, T_in, 1.0, fluid))

    @testset "Compressor" begin
        PR, η = 1.9, 0.80

        c_is = Compressor("C"; PR, η_poly = η, η_type = :isentropic)
        compute!(c_is, inlet())
        @test c_is.outlet[].Tt ≈ T_in * (1 + (PR^k - 1) / η) rtol = 1e-10

        c_po = Compressor("C"; PR, η_poly = η)   # default :polytropic
        compute!(c_po, inlet())
        @test c_po.outlet[].Tt ≈ T_in * PR^(k / η) rtol = 1e-3

        # Same numeric η: the polytropic model implies a lower adiabatic
        # efficiency for compression, hence a hotter outlet.
        @test c_po.outlet[].Tt > c_is.outlet[].Tt
        @test pressure_ratio(c_is) ≈ PR
    end

    @testset "Turbine" begin
        PR, η = 1.75, 0.87
        T_hot = 1144.0
        hot() = Port(FluidState(297.9e3, T_hot, 1.0, fluid))

        t_is = Turbine("T"; PR, η_poly = η, η_type = :isentropic)
        compute!(t_is, hot())
        @test t_is.outlet[].Tt ≈ T_hot * (1 - η * (1 - PR^(-k))) rtol = 1e-10

        t_po = Turbine("T"; PR, η_poly = η)
        compute!(t_po, hot())
        @test t_po.outlet[].Tt ≈ T_hot * PR^(-k * η) rtol = 1e-3

        # Same numeric η: the isentropic model extracts less work
        # (η_is > η_poly for the same physical expander).
        @test t_is.outlet[].Tt > t_po.outlet[].Tt

        # η_type also applies in :pressure_closure mode
        t_pc = Turbine("T"; mode = :pressure_closure, P_exit = 297.9e3 / PR,
                       η_poly = η, η_type = :isentropic)
        compute!(t_pc, hot())
        @test t_pc.outlet[].Tt ≈ t_is.outlet[].Tt rtol = 1e-12
    end

    @test_throws Exception Compressor("bad"; η_type = :adiabatic)
    @test_throws Exception Turbine("bad"; η_type = :secret)

    # FPT-backed fluid works through the isentropic path too (bisection
    # inversions).  Check the defining identity Δh_actual = Δh_is/η using the
    # table's own h and s, and sanity-check against the ideal-gas closed form
    # (HeXe84.fpt deviates from ideal monatomic by ~2% on the temperature
    # rise at these conditions, so only a loose agreement is expected).
    fpt = FPTFluid(joinpath(@__DIR__, "..", "HeXe84.fpt"))
    c = Compressor("C"; PR = 1.9, η_poly = 0.80, η_type = :isentropic)
    compute!(c, Port(FluidState(163.4e3, 300.0, 0.6, fpt)))
    s_in  = entropy(fpt, 300.0, 163.4e3)
    T_is  = T_from_s(fpt, s_in, 1.9 * 163.4e3; T_guess = 300.0)
    Δh_is = enthalpy(fpt, T_is, 1.9 * 163.4e3) - enthalpy(fpt, 300.0, 163.4e3)
    @test specific_work(c) ≈ Δh_is / 0.80 rtol = 1e-6
    @test c.outlet[].Tt ≈ 300.0 * (1 + (1.9^k - 1) / 0.80) atol = 10.0
end
