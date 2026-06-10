using Test
using GasCycle

# Analytic speed lines the fixture file was generated from
_Wc(n, r)  = n * (1.1 - 0.12 * (r - 2))
_PR(n, r)  = 1 + 1.8 * n^2 * (1 + 0.18 * (r - 2))
_eff(n, r) = 0.85 - 0.05 * (r - 2)^2 - 0.1 * (n - 1)^2

@testset "NPSS .map reader" begin
    path = joinpath(@__DIR__, "data", "synthetic_compressor.map")
    tables = read_npss_map(path)

    # ── Parse fidelity ───────────────────────────────────────────────────────
    @test sort(collect(keys(tables))) == ["TB_PR", "TB_Wc", "TB_eff"]

    tWc = tables["TB_Wc"]
    @test tWc.argnames == ["ALPHA", "NcMap", "RlineMap"]
    @test tWc.quantity == "Wc"
    @test tWc.alphas == [0.0]
    @test tWc.speeds[1] == [0.6, 0.8, 1.0, 1.1]
    @test tWc.coords[1][1] == [1.0, 1.5, 2.0, 2.5, 3.0]   # wrapped across lines

    for (t, f) in [(tables["TB_Wc"], _Wc), (tables["TB_PR"], _PR),
                   (tables["TB_eff"], _eff)]
        for (i, n) in enumerate(t.speeds[1]), (j, r) in enumerate(t.coords[1][i])
            @test t.values[1][i][j] ≈ f(n, r) atol = 1e-6
        end
    end

    # ── Conversion to PerformanceMap ─────────────────────────────────────────
    pmap = to_performance_map(tables)
    @test pmap isa PerformanceMap
    @test pmap.Nc_axis == [0.6, 0.8, 1.0, 1.1]

    # At tabulated Rline nodes the conversion must be exact (the fixture
    # rounds to 6 digits, hence the tolerance and the clamp at the corner).
    inmap(w) = clamp(w, first(pmap.Wc_axis), last(pmap.Wc_axis))
    for n in pmap.Nc_axis, r in [1.0, 1.5, 2.0, 2.5, 3.0]
        PR_q, eff_q = query(pmap, n, inmap(_Wc(n, r)))
        @test PR_q ≈ _PR(n, r) atol = 1e-5
        @test eff_q ≈ _eff(n, r) atol = 1e-5
    end

    # Between Rline nodes the map file itself is the resolution limit: the
    # converted map must reproduce linear interpolation of the source data
    # (which is what NPSS computes from the same table).  PR and Wc are both
    # affine in Rline here, so PR matches the analytic formula too; η is
    # quadratic, so it matches the source chord, not the formula.
    lin(r, r0, r1, y0, y1) = y0 + (r - r0) / (r1 - r0) * (y1 - y0)
    for n in pmap.Nc_axis, (r, r0, r1) in [(1.2, 1.0, 1.5), (2.7, 2.5, 3.0)]
        PR_q, eff_q = query(pmap, n, inmap(_Wc(n, r)))
        @test PR_q ≈ _PR(n, r) atol = 1e-5
        @test eff_q ≈ lin(r, r0, r1, _eff(n, r0), _eff(n, r1)) atol = 1e-5
    end

    # Beyond a speed line's own flow range (but inside the map's), values
    # clamp to that line's end point.
    PR_hi, _ = query(pmap, 0.6, _Wc(0.6, 1.0) + 0.05)
    @test PR_hi ≈ _PR(0.6, 1.0) atol = 1e-9

    # The converted map plugs into the standard design-point scaling workflow
    scaled = scale_map(pmap; Nc_des = 36000.0, Wc_des = 0.5, PR_des = 1.9,
                       eta_des = 0.80, Nc_ref = 1.0, Wc_ref = _Wc(1.0, 2.0))
    PR_s, eta_s = query(scaled, 36000.0, 0.5)
    @test PR_s ≈ 1.9 atol = 1e-9
    @test eta_s ≈ 0.80 atol = 1e-9

    # ── Error paths ──────────────────────────────────────────────────────────
    @test_throws Exception to_performance_map(Dict("TB_Wc" => tables["TB_Wc"]))

    # ── 2-level tables (no ALPHA) normalize to a single alpha ────────────────
    mktemp() do tmppath, io
        write(io, """
        // two-argument table, single-line arrays
        Table TB_Wc(real NcMap, real RlineMap) {
            NcMap = 1.0 { RlineMap = { 1.0, 2.0 }  Wc = { 0.9, 0.8 } }
        }
        """)
        close(io)
        t2 = read_npss_map(tmppath)["TB_Wc"]
        @test t2.alphas == [0.0]
        @test t2.speeds[1] == [1.0]
        @test t2.values[1][1] == [0.9, 0.8]
    end
end
