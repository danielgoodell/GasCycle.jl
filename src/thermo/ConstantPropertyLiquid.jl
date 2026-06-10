"""
Constant-property incompressible-liquid backend for coolant loops
(water, oil) on the cold side of a heat-rejection exchanger.

  h(T) = cp · T          (referenced to 0 K, like IdealGasFluid)
  s(T) = cp · ln(T / 298.15 K)
  ρ    = constant
  γ    = 1               (incompressible: Cp ≈ Cv)

Pressure does not enter any property — liquid coolant loops in this code
only carry sensible heat.  All methods are closed-form, so ForwardDiff
Dual numbers propagate exactly.

Construct directly in SI units:

    water = ConstantPropertyLiquid(cp = 4186.8, rho = 999.0, name = "H2O")

or read the constants from an NPSS function-style FPT file (the
collaborator's `H2O.fpt` / `Oil.fpt`, which define `real Cp(...)` /
`real rho(...)` returning constants in BTU/lbm·R and lbm/ft³):

    water = ConstantPropertyLiquid("H2O.fpt")
"""
struct ConstantPropertyLiquid <: FluidProperties
    cp_val::Float64  # specific heat [J/(kg·K)]
    rho::Float64     # density [kg/m³]
    name::String
end

function ConstantPropertyLiquid(; cp::Real, rho::Real, name::String = "liquid")
    cp > 0  || error("ConstantPropertyLiquid: cp must be positive, got $cp")
    rho > 0 || error("ConstantPropertyLiquid: rho must be positive, got $rho")
    ConstantPropertyLiquid(Float64(cp), Float64(rho), name)
end

cp(fp::ConstantPropertyLiquid, T, P)       = fp.cp_val
enthalpy(fp::ConstantPropertyLiquid, T, P) = fp.cp_val * T
entropy(fp::ConstantPropertyLiquid, T, P)  = fp.cp_val * log(T / 298.15)
density(fp::ConstantPropertyLiquid, T, P)  = fp.rho
gamma(fp::ConstantPropertyLiquid, T, P)    = 1.0

# Exact closed-form inversions — fully AD-compatible
T_from_h(fp::ConstantPropertyLiquid, h_target, P; T_guess=500.0) =
    h_target / fp.cp_val

T_from_s(fp::ConstantPropertyLiquid, s_target, P; T_guess=500.0) =
    298.15 * exp(s_target / fp.cp_val)

# ── Function-style FPT reader ─────────────────────────────────────────────────
#
# Unlike the Table-format files handled by FPTFluid, liquid coolant FPT
# files define properties as NPSS functions.  We only support the
# constant-property case: a function whose body returns a numeric literal,
# or delegates with a single call to another function that does (Oil.fpt's
# `Cp` returns `Cpt(T)`, and `Cpt` returns 0.8).

"""
    _fpt_function_bodies(content) -> Dict{String,String}

Extract `real name(args) { body }` blocks, tracking brace depth so nested
braces inside a body (do/while, if) don't truncate it.
"""
function _fpt_function_bodies(content::AbstractString)
    bodies = Dict{String,String}()
    for m in eachmatch(r"real\s+(\w+)\s*\([^)]*\)\s*\{", content)
        start = m.offset + ncodeunits(m.match)
        depth = 1
        i = start
        while i <= lastindex(content) && depth > 0
            c = content[i]
            c == '{' && (depth += 1)
            c == '}' && (depth -= 1)
            i = nextind(content, i)
        end
        bodies[m.captures[1]] = content[start:prevind(content, i, 2)]
    end
    bodies
end

"""
    _fpt_constant_return(bodies, fname; _depth) -> Float64 or nothing

Resolve a function to a constant: its first `return` statement is either a
numeric literal or a call `g(...)` to another function that resolves to a
constant.  Returns `nothing` if neither applies.
"""
function _fpt_constant_return(bodies::Dict{String,String}, fname::String; _depth::Int = 0)
    _depth > 4 && return nothing
    haskey(bodies, fname) || return nothing
    m = match(r"return\s+([^;]+);", bodies[fname])
    isnothing(m) && return nothing
    expr = strip(m.captures[1])

    val = tryparse(Float64, expr)
    isnothing(val) || return val

    mc = match(r"^(\w+)\s*\(", expr)
    isnothing(mc) && return nothing
    _fpt_constant_return(bodies, String(mc.captures[1]); _depth = _depth + 1)
end

"""
    ConstantPropertyLiquid(path) -> ConstantPropertyLiquid

Read constant `Cp` [BTU/(lbm·R)] and `rho` [lbm/ft³] from an NPSS
function-style FPT file and convert to SI.  Errors if either function is
missing or not reducible to a constant (a temperature-dependent liquid
table needs a different backend).
"""
function ConstantPropertyLiquid(path::String)
    content = read(path, String)
    bodies  = _fpt_function_bodies(content)

    cp_btu = _fpt_constant_return(bodies, "Cp")
    isnothing(cp_btu) && error(
        "ConstantPropertyLiquid: no constant Cp function found in $path")
    rho_lb = _fpt_constant_return(bodies, "rho")
    isnothing(rho_lb) && error(
        "ConstantPropertyLiquid: no constant rho function found in $path")

    ConstantPropertyLiquid(cp   = btulbmR_to_JkgK(cp_btu),
                           rho  = lbmft3_to_kgm3(rho_lb),
                           name = splitext(basename(path))[1])
end
