using Test
using GasCycle

const _DATA = joinpath(@__DIR__, "..", "data")

@testset "NPSS map reader" begin
    # ── Real R-line compressor map ────────────────────────────────────────────
    @testset "compressor map parse + eval" begin
        cm = compressor_map(joinpath(_DATA, "compressor_argon.map"))
        @test cm isa CompressorMap
        @test cm.NcMapDes    == 1.0
        @test cm.RlineMapDes == 2.212
        @test cm.RlineStall  == 1.0
        @test cm.re !== nothing            # S_Re subelement parsed

        # 6 speed lines (NcMap 0.5 … 1.0), 11 R-line nodes each
        @test length(cm.flow.speeds[1]) == 6
        @test cm.flow.speeds[1] == [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        @test cm.flow.coords[1][1] == collect(1.0:0.2:3.0)

        # Forward eval reproduces tabulated nodes exactly (unscaled map, s_*=1).
        # NcMap=1.0, RlineMap=2.20 → WcMap 0.9977, PRmap 1.9137, effMap 0.7995.
        o = eval_map(cm, 1.0, 2.20)
        @test o.Wc  ≈ 0.9977 atol = 1e-4
        @test o.PR  ≈ 1.9137 atol = 1e-4
        @test o.eff ≈ 0.7995 atol = 1e-4

        # Scaling lands the design point exactly at the map's design anchors.
        sc = scale_map(cm; Nc_des = 36000.0, Wc_des = 0.581, PR_des = 1.9, eta_des = 0.795)
        d  = eval_map(sc, 36000.0, sc.RlineMapDes)
        @test d.Wc  ≈ 0.581 atol = 1e-9
        @test d.PR  ≈ 1.9   atol = 1e-9
        @test d.eff ≈ 0.795 atol = 1e-9
    end

    # ── Real PR-parameterized turbine map ─────────────────────────────────────
    @testset "turbine map parse + eval" begin
        tm = turbine_map(joinpath(_DATA, "turbine_argon_tt.map"))
        @test tm isa TurbineMap
        @test tm.NpMapDes == 1.0
        @test tm.PRmapDes == 1.645
        @test tm.re !== nothing

        @test length(tm.flow.speeds[1]) == 6     # NpMap 0.3 … 1.1
        # NpMap=1.0, PRmap=1.661 → WpMap 1.0062, effMap 0.9154.
        o = eval_map(tm, 1.0, 1.661)
        @test o.Wp  ≈ 1.0062 atol = 1e-4
        @test o.eff ≈ 0.9154 atol = 1e-4

        sc = scale_map(tm; Np_des = 30000.0, Wp_des = 0.486, PR_des = 1.645, eta_des = 0.91)
        d  = eval_map(sc, 30000.0, 1.645)
        @test d.Wp  ≈ 0.486 atol = 1e-9
        @test d.eff ≈ 0.91  atol = 1e-9
    end

    # ── Reynolds correction (parsed S_Re tables) ──────────────────────────────
    @testset "Reynolds correction" begin
        cm = compressor_map(joinpath(_DATA, "compressor_argon.map"))
        # rc === nothing ⇒ no correction (default element behavior).
        @test eval_map(cm, 1.0, 2.20, nothing).eff ≈ eval_map(cm, 1.0, 2.20).eff
        # At the table reference (RNI = 1.0) the factors are unity; at RNI=1.5461
        # the digitized s_effRe = 1.00594 multiplies efficiency, Wc unchanged.
        base = eval_map(cm, 1.0, 2.20, 1.0)
        bump = eval_map(cm, 1.0, 2.20, 1.5461)
        @test base.eff ≈ eval_map(cm, 1.0, 2.20).eff
        @test bump.eff / base.eff ≈ 1.00594 atol = 1e-5
        @test bump.Wc ≈ base.Wc                     # s_WcRe ≡ 1
    end

    # ── Name-agnostic: NPSS cares about hierarchy, not names ──────────────────
    @testset "positional / name-agnostic parse" begin
        mktemp() do path, io
            # Scrambled axis (ALPHA/SPED/R) and leaf (FLOW/PRES/EF) names.
            write(io, """
            Subelement CompressorRlineMap S_map {
              NcMapDes = 1.0;  RlineMapDes = 2.0;
              Table TB_Wc(real ALPHA, real SPED, real R) {
                ALPHA = 0.0 {
                  SPED = 1.0 { R = { 1.0, 2.0, 3.0 }  FLOW = { 0.5, 1.0, 1.5 } }
                }
                SPED.interp = "linear";  R.interp = "linear";
              }
              Table TB_PR(real ALPHA, real SPED, real R) {
                ALPHA = 0.0 {
                  SPED = 1.0 { R = { 1.0, 2.0, 3.0 }  PRES = { 2.0, 1.9, 1.8 } }
                }
                R.interp = "linear";
              }
              Table TB_eff(real ALPHA, real SPED, real R) {
                ALPHA = 0.0 {
                  SPED = 1.0 { R = { 1.0, 2.0, 3.0 }  EF = { 0.80, 0.85, 0.82 } }
                }
                R.interp = "linear";
              }
            }
            """)
            close(io)
            cm = compressor_map(path)
            @test cm.flow.coords[1][1] == [1.0, 2.0, 3.0]
            o = eval_map(cm, 1.0, 2.0)                # mid R-line node
            @test (o.Wc, o.PR, o.eff) == (1.0, 1.9, 0.85)
        end
    end

    # ── Error paths ───────────────────────────────────────────────────────────
    @testset "errors" begin
        mktemp() do path, io
            write(io, "Subelement X S { Table TB_Wc(real NcMap, real RlineMap) { " *
                      "NcMap = 1.0 { RlineMap = {1.0} WcMap = {0.5} } } }")
            close(io)
            @test_throws Exception compressor_map(path)   # missing TB_PR/TB_eff
        end
    end
end
