using Test
using GasCycle

@testset "Unit conversions" begin
    # Round trips
    @test K_to_R(R_to_K(540.0)) ≈ 540.0
    @test Pa_to_psia(psia_to_Pa(23.7)) ≈ 23.7
    @test kgps_to_lbps(lbps_to_kgps(1.32)) ≈ 1.32
    @test kg_to_lbm(lbm_to_kg(1.0)) ≈ 1.0
    @test Jkg_to_btulbm(btulbm_to_Jkg(100.0)) ≈ 100.0
    @test JkgK_to_btulbmR(btulbmR_to_JkgK(0.06)) ≈ 0.06
    @test kgm3_to_lbmft3(lbmft3_to_kgm3(0.5)) ≈ 0.5
    @test radps_to_rpm(rpm_to_radps(36000.0)) ≈ 36000.0
    @test W_to_hp(hp_to_W(10.0)) ≈ 10.0

    # Known values (BRU design point and standard definitions)
    @test R_to_K(540.0) ≈ 300.0
    @test R_to_K(2060.0) ≈ 1144.44 atol = 0.01
    @test psia_to_Pa(14.695948775) ≈ 101325.0 atol = 0.1   # 1 standard atm
    @test psia_to_Pa(23.7) ≈ 163405.7 atol = 0.5
    @test lbps_to_kgps(1.32) ≈ 0.59874 atol = 1e-4
    @test lbm_to_kg(1.0) ≈ 0.45359237
    @test btulbm_to_Jkg(1.0) ≈ 2326.0
    @test btulbmR_to_JkgK(1.0) ≈ 4186.8
    @test lbmft3_to_kgm3(1.0) ≈ 16.01846 atol = 1e-4
    @test rpm_to_radps(60.0) ≈ 2π   # 1 rev/s = 2π rad/s
    @test hp_to_W(1.0) ≈ 745.69987 atol = 1e-4
end
