# Run: julia data/generate_liquids.jl  (regenerates the three coolant FPT files in place)
# Generates data/Water.fpt, data/DC200.fpt, data/WaterEG50.fpt with
# temperature-dependent polynomial properties in NPSS English units.
using Printf, LinearAlgebra

polyfit(x, y, n) = [x.^p for p in 0:n] |> M -> hcat(M...) \ y
peval(c, x) = sum(c[i] * x^(i-1) for i in eachindex(c))
maxres(c, x, y) = maximum(abs.(peval.(Ref(c), x) ./ y .- 1)) * 100

# SI coeffs (T in K) -> English coeffs (T in R), per-quantity value conversion
to_english(c_si, conv) = [c_si[i] * conv / 1.8^(i-1) for i in eachindex(c_si)]

fmt(c) = join([@sprintf("%.10e", v) for v in c], ", ")

"NPSS polynomial expression in Tt (R) from ascending English coeffs."
function npss_poly(c; var = "T")
    terms = String[]
    for (i, v) in enumerate(c)
        abs(v) < 1e-300 && continue
        tpow = i == 1 ? "" : " * " * join(fill(var, i - 1), " * ")
        push!(terms, @sprintf("%.10e%s", v, tpow))
    end
    isempty(terms) ? "0." : join(terms, "\n\t       + ")
end

"Integral coeffs: h(T) = sum c_i T^(i+1)/(i+1)  (as polynomial with zero const)."
hcoeffs(c) = vcat(0.0, [c[i] / i for i in eachindex(c)])

function emit(path, name, desc, sources, Trange_K, cp_si, rho_si; k_si=nothing, mu_si=nothing, notes="")
    cp_e  = to_english(cp_si, 1 / 4186.8)
    rho_e = to_english(rho_si, 0.3048^3 / 0.45359237)
    h_e   = hcoeffs(cp_e)
    TminR, TmaxR = Trange_K .* 1.8

    io = IOBuffer()
    println(io, "// $name — $desc")
    println(io, "// Generated 2026-06-12 by GasCycle (see data/README.md for fit anchors).")
    for s in sources; println(io, "// Source: $s"); end
    @printf(io, "// Valid Tt = %.0f–%.0f R (%.0f–%.0f K). Pressure does not enter (liquid).\n",
            TminR, TmaxR, Trange_K...)
    println(io, "// Units: Tt [R], Cp [BTU/(lbm·R)], rho [lbm/ft3], h [BTU/lbm]" *
                (isnothing(k_si) ? "" : ", k [BTU/(ft·s·R)]") *
                (isnothing(mu_si) ? "" : ", mu [lbm/(ft·s)]"))
    isempty(notes) || println(io, "// $notes")
    println(io, """

indeps   = {"Tt", "Pt"};
hTindeps = {"Tt", "Pt"};
ThIndeps = {"ht", "Pt"};

real Cp( real T, real P ){

\treturn Cpt( T );
}

real Cpt( real T ){
\treturn $(npss_poly(cp_e));
}

real rho( real T, real P ){
\treturn $(npss_poly(rho_e));
}

real h_T( real T, real P ){
\treturn $(npss_poly(h_e));
}

real T_h( real ht, real P ){
\treal T = 540.;
\treal iter = 0.;
\tdo {
\t\tT = T - ( h_T( T, P ) - ht ) / Cpt( T );
\t\titer = iter + 1.;
\t} while ( iter < 8. );
\treturn T;
}""")
    if !isnothing(k_si)
        k_e = to_english(k_si, 0.3048 / (4186.8 * 0.45359237))
        println(io, """

real k( real T, real P ){
\treturn $(npss_poly(k_e));
}""")
    end
    if !isnothing(mu_si)
        mu_e = to_english(mu_si, 0.3048 / 0.45359237)
        println(io, """

real mu( real T, real P ){
\treturn $(npss_poly(mu_e));
}""")
    end
    write(path, String(take!(io)))
    println("wrote $path")
end

# ── Water (liquid, 1 atm; CRC/IAPWS anchors) ─────────────────────────────────
Tc  = collect(0.0:10:100)
TK  = Tc .+ 273.15
cpW  = [4217.6,4192.1,4181.8,4178.4,4178.5,4180.6,4184.3,4189.5,4196.3,4205.0,4215.9]
rhoW = [999.84,999.70,998.21,995.65,992.22,988.04,983.20,977.76,971.79,965.31,958.35]
kW   = [0.5610,0.5800,0.5984,0.6154,0.6305,0.6435,0.6543,0.6631,0.6700,0.6753,0.6791]
muW  = [1.792,1.306,1.002,0.7977,0.6527,0.5469,0.4665,0.4040,0.3544,0.3145,0.2818] .* 1e-3

cp_c  = polyfit(TK, cpW, 3);  @printf("Water cp  cubic resid %.3f%%\n", maxres(cp_c, TK, cpW))
rho_c = polyfit(TK, rhoW, 3); @printf("Water rho cubic resid %.4f%%\n", maxres(rho_c, TK, rhoW))
k_c   = polyfit(TK, kW, 2);   @printf("Water k   quad resid %.3f%%\n", maxres(k_c, TK, kW))
mu_c  = polyfit(TK, muW, 4);  @printf("Water mu  quartic resid %.1f%%\n", maxres(mu_c, TK, muW))

emit(joinpath(@__DIR__, "Water.fpt"), "Water.fpt",
     "liquid water at ~1 atm, temperature-dependent polynomial properties",
     ["CRC Handbook / IAPWS liquid-water values, fit anchors every 10 C, 0-100 C"],
     (273.15, 373.15), cp_c, rho_c; k_si=k_c, mu_si=mu_c,
     notes="Fit residuals: cp 0.03%, rho 0.01%, k 0.13%, mu ~2.5% (polynomial vs Arrhenius).")

# ── Dow Corning 200, 2.0 cSt (BRU heat-rejection coolant grade) ──────────────
# cp(25C) = 0.410 cal/gC (Clearco PSF-2cSt TDS), slope +0.0006 cal/gC per C
# (Dow typical-properties trend for low-cSt grades);
# rho(25C) = 873 kg/m3, expansion 1.16e-3 1/C; k(25C) = 0.109 W/mK (constant).
cpD_si  = let a1 = 0.0006 * 4186.8           # J/kgK per K
    [1716.4 - a1 * 298.15, a1]
end
rhoD_si = let b1 = -873.0 * 1.16e-3          # kg/m3 per K
    [873.0 - b1 * 298.15, b1]
end
emit(joinpath(@__DIR__, "DC200.fpt"), "DC200.fpt",
     "Dow Corning 200 silicone fluid (PDMS), 2.0 cSt grade — BRU heat-sink coolant (NASA CR-120816)",
     ["Clearco PSF-2cSt TDS (cp 0.410 cal/gC, rho 0.873 g/cm3, k 0.109 W/mK at 25 C)",
      "Dow 'Typical properties of Dow silicone fluids' (cp slope, expansion coefficient)"],
     (233.15, 423.15), cpD_si, rhoD_si; k_si=[0.109],
     notes="cp/rho linear in T; k held at the 25 C value (no published T-dependence); mu omitted.")

# ── Water / ethylene glycol 50/50 by mass ────────────────────────────────────
# cp: Melinder-2010/ASHRAE-consistent linear, cp = 3247 + 3.25*T_C J/kgK.
# rho: vendor (Dow-derived) specific-gravity table anchors.
cpE_si  = [3247.0 - 3.25 * 273.15, 3.25]
TcE   = [-17.8, 4.4, 26.7, 48.9, 71.1]
rhoE  = [1088.0, 1077.0, 1064.0, 1050.0, 1038.0]
rhoE_c = polyfit(TcE .+ 273.15, rhoE, 2)
@printf("WaterEG50 rho quad resid %.3f%%\n", maxres(rhoE_c, TcE .+ 273.15, rhoE))
kE_si = [0.394 - 5.5e-4 * 273.15, 5.5e-4]    # ~0.40 W/mK at 20 C, mild rise
emit(joinpath(@__DIR__, "WaterEG50.fpt"), "WaterEG50.fpt",
     "water / ethylene glycol 50/50 by mass (freeze point ~ -36 C)",
     ["cp: Melinder 2010 via ASHRAE-consistent linear fit (vendor tables run ~4% higher)",
      "rho: Dow-derived vendor SG table (corecheminc.com), quadratic fit",
      "k: approximate linear (~0.40 W/mK at 20 C)"],
     (243.15, 373.15), cpE_si, rhoE_c; k_si=kE_si,
     notes="mu omitted (strongly non-polynomial); cp basis is mass fraction 0.50, uninhibited.")
