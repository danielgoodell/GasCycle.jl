using Test
using GasCycle
using ForwardDiff
import GasCycle: cp

"""
PolynomialLiquid: temperature-dependent coolant backend and the generated
data/Water.fpt, data/DC200.fpt, data/WaterEG50.fpt files.

Spot-value oracles are the fit anchor data documented in data/README.md
(CRC/IAPWS water; Clearco/Dow DC-200 2 cSt; Melinder/vendor 50/50 water-EG),
checked to tolerances just above the documented fit residuals.
"""

const _DATA = joinpath(@__DIR__, "..", "data")

@testset "PolynomialLiquid — file parsing and spot values" begin
    water = PolynomialLiquid(joinpath(_DATA, "Water.fpt"))
    dc    = PolynomialLiquid(joinpath(_DATA, "DC200.fpt"))
    eg    = PolynomialLiquid(joinpath(_DATA, "WaterEG50.fpt"))
    P = 3e5

    # Water vs CRC/IAPWS anchors
    @test cp(water, 293.15, P)      ≈ 4181.8  rtol = 2e-3
    @test cp(water, 353.15, P)      ≈ 4196.3  rtol = 2e-3
    @test density(water, 293.15, P) ≈ 998.21  rtol = 5e-4
    @test density(water, 363.15, P) ≈ 965.31  rtol = 5e-4
    @test conductivity(water, 313.15, P) ≈ 0.6305 rtol = 3e-3
    @test viscosity(water, 293.15, P) ≈ 1.002e-3 rtol = 4e-2
    @test prandtl(water, 293.15, P) ≈ 4181.8 * 1.002e-3 / 0.5984 rtol = 5e-2
    # cp(T) actually varies: the 0→35→100 °C dip-and-rise shape
    @test cp(water, 273.65, P) > cp(water, 308.15, P)
    @test cp(water, 372.65, P) > cp(water, 308.15, P)
    # Advisory range parsed from the header comment
    @test water.T_min ≈ 273.15 atol = 0.5
    @test water.T_max ≈ 373.15 atol = 0.5

    # DC-200 2 cSt vs Clearco TDS at 25 °C: cp 0.410 cal/(g·°C), ρ 873 kg/m³
    @test cp(dc, 298.15, P)      ≈ 0.410 * 4186.8 rtol = 1e-3
    @test density(dc, 298.15, P) ≈ 873.0 rtol = 1e-3
    @test conductivity(dc, 298.15, P) ≈ 0.109 rtol = 1e-3
    @test cp(dc, 358.15, P) > cp(dc, 298.15, P)        # cp rises with T
    @test density(dc, 358.15, P) < density(dc, 298.15, P)
    @test_throws Exception viscosity(dc, 298.15, P)     # no mu data in file

    # 50/50 water-EG vs Melinder-consistent cp and vendor SG table
    @test cp(eg, 293.15, P)      ≈ 3247.0 + 3.25 * 20 rtol = 1e-3
    @test density(eg, 299.85, P) ≈ 1064.0 rtol = 2e-3   # 26.7 °C vendor anchor
    @test density(eg, 255.35, P) ≈ 1088.0 rtol = 2e-3   # −17.8 °C
    @test_throws Exception viscosity(eg, 293.15, P)

    # Gas-side interface conventions
    @test gamma(water, 300.0, P) == 1.0
    @test enthalpy(water, 0.0, P) == 0.0                # 0 K reference
end

@testset "PolynomialLiquid — thermodynamic consistency and inversions" begin
    water = PolynomialLiquid(joinpath(_DATA, "Water.fpt"))
    P = 3e5

    # dh/dT = cp and ds/dT = cp/T (exact identities of the integration)
    for T in (280.0, 320.0, 360.0)
        dh = (enthalpy(water, T+1e-3, P) - enthalpy(water, T-1e-3, P)) / 2e-3
        @test dh ≈ cp(water, T, P) rtol = 1e-8
        ds = (entropy(water, T+1e-3, P) - entropy(water, T-1e-3, P)) / 2e-3
        @test ds ≈ cp(water, T, P) / T rtol = 1e-8
    end
    @test entropy(water, 298.15, P) == 0.0               # reference convention

    # Newton inversion round trips
    for T in (275.0, 310.0, 370.0)
        @test T_from_h(water, enthalpy(water, T, P), P; T_guess = 300.0) ≈ T rtol = 1e-9
        @test T_from_s(water, entropy(water, T, P), P; T_guess = 300.0) ≈ T rtol = 1e-9
    end

    # AD through the inversion: dT/dh = 1/cp via the implicit function theorem
    h0 = enthalpy(water, 330.0, P)
    dTdh = ForwardDiff.derivative(h -> T_from_h(water, h, P; T_guess = 300.0), h0)
    @test dTdh ≈ 1 / cp(water, 330.0, P) rtol = 1e-8

    # Constant-coefficient degenerate case matches ConstantPropertyLiquid
    cpl  = ConstantPropertyLiquid(cp = 4186.8, rho = 999.0)
    poly = PolynomialLiquid(cp = 4186.8, rho = 999.0, name = "const")
    @test enthalpy(poly, 350.0, P) == enthalpy(cpl, 350.0, P)
    @test entropy(poly, 350.0, P)  ≈ entropy(cpl, 350.0, P) rtol = 1e-12

    # Legacy constant files load through the polynomial parser too (Oil.fpt
    # delegates Cp → Cpt → 0.8)
    oil = PolynomialLiquid(joinpath(_DATA, "Oil.fpt"))
    @test cp(oil, 400.0, P) ≈ btulbmR_to_JkgK(0.8)
    @test density(oil, 400.0, P) ≈ lbmft3_to_kgm3(62.424)

    # Constructor guards
    @test_throws Exception PolynomialLiquid(cp = [-1.0], rho = 999.0)
    @test_throws Exception PolynomialLiquid(cp = [5000.0, -20.0], rho = 999.0,
                                            T_min = 233.15, T_max = 423.15)
end

@testset "PolynomialLiquid — coolant in a closed loop" begin
    gas = IdealGasFluid(M_molar = 83.8)
    eg  = PolynomialLiquid(joinpath(_DATA, "WaterEG50.fpt"))

    net  = FlowNetwork()
    duct = Duct("Inlet"; dPqP = 0.0)
    cool = HeatExchanger("Cooler"; UA = 400.0, dPqP_hot = 0.01, dPqP_cold = 0.005)
    add!(net, duct, cool)
    connect_port!(net, duct, :outlet, cool, :hot_inlet)
    set_state!(net, duct; Pt = 165e3, Tt = 440.0, W = 0.6, fluid = gas)
    set_boundary!(net, cool, :cold_inlet; Pt = 3e5, Tt = 280.0, W = 1.0, fluid = eg)
    sol = solve!(net)
    @test sol.status == :success

    # Energy balance: gas enthalpy drop equals coolant enthalpy rise,
    # with the coolant side evaluated through the T-dependent polynomials.
    Q_gas = Q_transferred(cool)
    T_c_out = cool.cold_outlet[].Tt
    Q_cool = 1.0 * (enthalpy(eg, T_c_out, 3e5) - enthalpy(eg, 280.0, 3e5))
    @test Q_gas ≈ Q_cool rtol = 1e-8
    @test 280.0 < T_c_out < 440.0
end
