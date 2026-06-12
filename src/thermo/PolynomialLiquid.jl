"""
Polynomial-property incompressible-liquid backend — temperature-dependent
generalization of `ConstantPropertyLiquid` for coolant loops.

  cp(T) = Σ aᵢ Tⁱ                       [J/(kg·K)], T in K
  h(T)  = Σ aᵢ Tⁱ⁺¹/(i+1)               (referenced to 0 K)
  s(T)  = a₀·ln(T/298.15) + Σᵢ≥₁ aᵢ (Tⁱ − 298.15ⁱ)/i
  ρ(T), μ(T), k(T) = polynomials        (μ/k optional)
  γ     = 1

Pressure does not enter any property. All evaluations and the Newton
inversions are plain arithmetic, so ForwardDiff Dual numbers propagate
exactly (the inversions converge in Dual arithmetic, which carries the
implicit-function derivative).

Construct directly in SI units (coefficients ascending in T [K], scalars
mean constant):

    eg = PolynomialLiquid(cp = [2359.2, 3.25], rho = [1208.7, -0.63, ...],
                          name = "WaterEG50")

or read an NPSS function-style FPT file whose `Cp`/`Cpt`, `rho` (and
optionally `mu`, `k`) functions return polynomials in Tt [°R] with English
units (BTU/(lbm·°R), lbm/ft³, lbm/(ft·s), BTU/(ft·s·°R)):

    eg = PolynomialLiquid("data/WaterEG50.fpt")

The generated files `data/Water.fpt`, `data/DC200.fpt`, `data/WaterEG50.fpt`
follow this format (provenance in `data/README.md`).
"""
struct PolynomialLiquid <: FluidProperties
    cp_c::Vector{Float64}   # cp coefficients, ascending powers of T [K] → J/(kg·K)
    rho_c::Vector{Float64}  # density [kg/m³]
    mu_c::Vector{Float64}   # viscosity [Pa·s]; empty → not available
    k_c::Vector{Float64}    # conductivity [W/(m·K)]; empty → not available
    name::String
    T_min::Float64          # advisory validity range [K] (not enforced)
    T_max::Float64
end

_aspoly(x::Real) = [Float64(x)]
_aspoly(x::AbstractVector) = Vector{Float64}(x)

function PolynomialLiquid(; cp, rho, mu = Float64[], k = Float64[],
                          name::String = "liquid",
                          T_min::Real = 233.15, T_max::Real = 423.15)
    cp_c, rho_c = _aspoly(cp), _aspoly(rho)
    isempty(cp_c) && error("PolynomialLiquid: cp coefficients required")
    isempty(rho_c) && error("PolynomialLiquid: rho coefficients required")
    for Tq in (T_min, T_max)
        evalpoly(Tq, cp_c)  > 0 || error("PolynomialLiquid \"$name\": cp ≤ 0 at $Tq K")
        evalpoly(Tq, rho_c) > 0 || error("PolynomialLiquid \"$name\": rho ≤ 0 at $Tq K")
    end
    PolynomialLiquid(cp_c, rho_c, _aspoly(mu), _aspoly(k), name,
                     Float64(T_min), Float64(T_max))
end

cp(fp::PolynomialLiquid, T, P)      = evalpoly(T, fp.cp_c)
density(fp::PolynomialLiquid, T, P) = evalpoly(T, fp.rho_c)
gamma(fp::PolynomialLiquid, T, P)   = 1.0

"""h(T) = Σ aᵢ Tⁱ⁺¹/(i+1), referenced to 0 K like the other backends."""
enthalpy(fp::PolynomialLiquid, T, P) =
    T * evalpoly(T, [c / i for (i, c) in enumerate(fp.cp_c)])

const _PL_S_TREF = 298.15  # entropy reference [K], matches ConstantPropertyLiquid

function entropy(fp::PolynomialLiquid, T, P)
    s = fp.cp_c[1] * log(T / _PL_S_TREF)
    for i in 2:length(fp.cp_c)
        s += fp.cp_c[i] * (T^(i-1) - _PL_S_TREF^(i-1)) / (i - 1)
    end
    s
end

function viscosity(fp::PolynomialLiquid, T, P)
    isempty(fp.mu_c) && error("viscosity not available for PolynomialLiquid \"$(fp.name)\"")
    evalpoly(T, fp.mu_c)
end

function conductivity(fp::PolynomialLiquid, T, P)
    isempty(fp.k_c) && error("conductivity not available for PolynomialLiquid \"$(fp.name)\"")
    evalpoly(T, fp.k_c)
end

prandtl(fp::PolynomialLiquid, T, P) =
    cp(fp, T, P) * viscosity(fp, T, P) / conductivity(fp, T, P)

# ── Newton inversions (cp > 0 ⟹ h and s strictly increasing in T) ────────────

_pl_val(x::Real) = x
_pl_val(x::ForwardDiff.Dual) = _pl_val(ForwardDiff.value(x))

function T_from_h(fp::PolynomialLiquid, h_target, P; T_guess = 350.0)
    T = T_guess + zero(h_target)
    for _ in 1:50
        dT = (enthalpy(fp, T, P) - h_target) / cp(fp, T, P)
        T -= dT
        _pl_val(T) < 1.0 && (T = one(T))
        abs(_pl_val(dT)) < 1e-10 * max(abs(_pl_val(T)), 1.0) && break
    end
    T
end

function T_from_s(fp::PolynomialLiquid, s_target, P; T_guess = 350.0)
    T = T_guess + zero(s_target)
    for _ in 1:50
        dT = (entropy(fp, T, P) - s_target) * T / cp(fp, T, P)
        T -= dT
        _pl_val(T) < 1.0 && (T = one(T))
        abs(_pl_val(dT)) < 1e-10 * max(abs(_pl_val(T)), 1.0) && break
    end
    T
end

# ── Function-style FPT reader (polynomial returns) ───────────────────────────
#
# Reuses _fpt_function_bodies from ConstantPropertyLiquid.jl.  A property
# function must return a polynomial in its temperature argument, e.g.
#   return 2.31e-01 + 3.33e-04 * T;     (any identifier counts as T)
# or delegate with a single call (`Cp` → `Cpt(T)`).  English units as in
# the file header; converted to SI with T in K on load.

"""
    _fpt_poly_return(bodies, fname) -> Vector{Float64} (ascending in T) or nothing
"""
function _fpt_poly_return(bodies::Dict{String,String}, fname::String; _depth::Int = 0)
    _depth > 4 && return nothing
    haskey(bodies, fname) || return nothing
    m = match(r"return\s+([^;]+);", bodies[fname])
    isnothing(m) && return nothing
    expr = replace(m.captures[1], r"\s+" => "")

    occursin(r"^NaN$"i, expr) && return nothing
    md = match(r"^(\w+)\(\w+\)$", expr)   # delegation: Cpt(T)
    isnothing(md) || return _fpt_poly_return(bodies, String(md.captures[1]);
                                             _depth = _depth + 1)
    _parse_fpt_poly(expr)
end

"""
Parse `c0+c1*T+c2*T*T` (whitespace already stripped; `+-c` allowed) into
ascending coefficients.  Any non-numeric factor is treated as the
temperature variable.  Returns nothing if a term is unparseable.
"""
function _parse_fpt_poly(expr::AbstractString)
    # Split on top-level '+' (not part of an exponent like 1e+05)
    terms = String[]
    buf = IOBuffer()
    prev = ' '
    for c in expr
        if c == '+' && !(prev in ('e', 'E'))
            push!(terms, String(take!(buf)))
        else
            write(buf, c)
        end
        prev = c
    end
    push!(terms, String(take!(buf)))

    coeffs = Float64[]
    for t in terms
        isempty(t) && continue
        coeff, degree = 1.0, 0
        for (j, f) in enumerate(split(t, '*'))
            v = tryparse(Float64, f)
            if !isnothing(v)
                coeff *= v
            elseif occursin(r"^\w+$", f)
                degree += 1
            elseif j == 1 && f == "-"     # bare leading minus
                coeff *= -1.0
            else
                return nothing
            end
        end
        length(coeffs) < degree + 1 && append!(coeffs, zeros(degree + 1 - length(coeffs)))
        coeffs[degree + 1] += coeff
    end
    isempty(coeffs) ? nothing : coeffs
end

"""English coeffs (T in °R) → SI coeffs (T in K): cᵢ_SI = conv · cᵢ_E · 1.8ⁱ."""
_pl_to_si(c::Vector{Float64}, conv::Float64) =
    [conv * c[i] * 1.8^(i-1) for i in eachindex(c)]

"""
    PolynomialLiquid(path) -> PolynomialLiquid

Read polynomial `Cp` [BTU/(lbm·°R)] and `rho` [lbm/ft³] — plus `mu`
[lbm/(ft·s)] and `k` [BTU/(ft·s·°R)] when present — from an NPSS
function-style FPT file and convert to SI. The advisory T range is taken
from a `// Valid Tt = <lo>–<hi> R` header comment when present.
"""
function PolynomialLiquid(path::String)
    content = read(path, String)
    bodies  = _fpt_function_bodies(content)

    cp_e = _fpt_poly_return(bodies, "Cp")
    isnothing(cp_e) && error("PolynomialLiquid: no polynomial Cp function in $path")
    rho_e = _fpt_poly_return(bodies, "rho")
    isnothing(rho_e) && error("PolynomialLiquid: no polynomial rho function in $path")
    mu_e = _fpt_poly_return(bodies, "mu")
    k_e  = _fpt_poly_return(bodies, "k")

    mr = match(r"Valid Tt = (\d+)[^\d](\d+) R", content)
    T_min, T_max = isnothing(mr) ? (233.15, 423.15) :
                   (parse(Float64, mr.captures[1]) / 1.8,
                    parse(Float64, mr.captures[2]) / 1.8)

    PolynomialLiquid(cp  = _pl_to_si(cp_e, btulbmR_to_JkgK(1.0)),
                     rho = _pl_to_si(rho_e, lbmft3_to_kgm3(1.0)),
                     mu  = isnothing(mu_e) ? Float64[] :
                           _pl_to_si(mu_e, 0.45359237 / 0.3048),
                     k   = isnothing(k_e) ? Float64[] :
                           _pl_to_si(k_e, btulbmR_to_JkgK(1.0) * 0.45359237 / 0.3048),
                     name = splitext(basename(path))[1],
                     T_min = T_min, T_max = T_max)
end
