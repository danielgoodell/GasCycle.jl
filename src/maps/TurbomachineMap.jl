"""
Turbomachinery performance-map interface — the intermediate layer between a
map's stored data and the Compressor/Turbine elements (the NPSS `Subelement` /
pyCycle component role).

An element never reads table data directly; it calls `eval_map` and receives
the physical (scaled, Reynolds-corrected) outputs it needs.  This decouples the
element from *how* the map is stored, so a map can be a tabulation loaded from
any file format, or a plain Julia function ("map as a script").

Maps are evaluated **forward in their native coordinates** — no inversion:

  Compressor (R-line native):  eval_map(m, Nc, Rline, rc) -> (; Wc, PR, eff)
  Turbine    (PR native):      eval_map(m, Np, PR,    rc) -> (; Wp, eff)

where `Nc`/`Np` is the physical corrected speed, the third positional argument
is the native *line coordinate* (R-line for compressors, expansion ratio PR for
turbines), and `rc` is the value looked up in the map's Reynolds-correction
tables (see `ReynoldsModel`).  `rc === nothing` disables the correction.

Everything is parametric in the numeric type so ForwardDiff Duals flow through
`eval_map` (AD is a hard constraint — interpolation is hand-rolled, no `Float64`
annotations on values derived from the query coordinates).
"""

# ── Corrected-variable conventions (NPSS) ─────────────────────────────────────
# Nc = N / sqrt(Tt / T_std);  Wc = W * sqrt(Tt / T_std) / Pt * P_std
const T_STD = 288.15   # K
const P_STD = 101325.0 # Pa

"Physical shaft speed N [rpm], inlet Tt [K] → corrected speed."
corrected_speed(N, Tt) = N / sqrt(Tt / T_STD)

"Mass flow W [kg/s], inlet Tt [K], Pt [Pa] → corrected mass flow."
corrected_flow(W, Tt, Pt) = W * sqrt(Tt / T_STD) / Pt * P_STD

# ── 1-D interpolation (lagrange3 / linear, linear or clamped extrapolation) ────

"Per-axis interpolation/extrapolation method, mirroring the NPSS `.interp`/`.extrap` declarations."
struct MapAxis
    interp::Symbol   # :lagrange3 | :linear
    extrap::Symbol   # :linear | :clamp
end
MapAxis() = MapAxis(:lagrange3, :linear)

_lin(x0, x1, y0, y1, xq) = y0 + (xq - x0) / (x1 - x0) * (y1 - y0)

# Cubic Lagrange through the four points xs[i0 .. i0+3].
function _lagrange4(xs, ys, i0::Int, xq)
    s = zero(promote_type(typeof(xq), eltype(ys)))
    @inbounds for i in i0:i0+3
        Li = one(s)
        for j in i0:i0+3
            i == j && continue
            Li *= (xq - xs[j]) / (xs[i] - xs[j])
        end
        s += ys[i] * Li
    end
    s
end

"""
    _interp1(xs, ys, xq, axis) -> value

1-D interpolation of `ys` over the ascending grid `xs` at query `xq`.
`xs`/`ys` are concrete numeric vectors; `xq` may be a Dual.  The interval search
branches on the *value* of `xq` only (the selected interval is non-differentiable,
the value within it is differentiable — standard for table interpolation).
"""
function _interp1(xs::AbstractVector{<:Real}, ys::AbstractVector, xq, axis::MapAxis)
    n = length(xs)
    n == 1 && return ys[1]
    if xq <= xs[1]
        return axis.extrap == :clamp ? ys[1] : _lin(xs[1], xs[2], ys[1], ys[2], xq)
    elseif xq >= xs[n]
        return axis.extrap == :clamp ? ys[n] : _lin(xs[n-1], xs[n], ys[n-1], ys[n], xq)
    end
    k = searchsortedlast(xs, xq)            # xs[k] <= xq < xs[k+1]
    if axis.interp == :linear || n < 4
        return _lin(xs[k], xs[k+1], ys[k], ys[k+1], xq)
    end
    i0 = clamp(k - 1, 1, n - 3)             # 4-point stencil around the interval
    _lagrange4(xs, ys, i0, xq)
end

# ── Tabulated quantity over a nested (alpha, speed, line) grid ─────────────────

"""
One tabulated map quantity (e.g. Wc, PR, eff) stored as the NPSS nested table

    alphas[a]            outer-axis values (single 0.0 for 2-arg tables)
    speeds[a][s]         speed values within alpha block a
    coords[a][s][:]      line-coordinate grid for that speed line
    values[a][s][:]      the quantity along that line

`axes = (alpha, speed, line)` carries the per-axis interpolation methods.  The
grids are name-agnostic: the loader fills them positionally from the table's
argument hierarchy (last arg = line coordinate), never from fixed names.
"""
struct MapTable
    alphas::Vector{Float64}
    speeds::Vector{Vector{Float64}}
    coords::Vector{Vector{Vector{Float64}}}
    values::Vector{Vector{Vector{Float64}}}
    axes::NTuple{3,MapAxis}
end

"Evaluate the nested table at (alpha α, speed n, line ℓ) by inside-out 1-D interps."
function eval_table(t::MapTable, α, n, ℓ)
    a_axis, n_axis, l_axis = t.axes
    per_alpha = map(eachindex(t.alphas)) do ia
        per_speed = map(eachindex(t.speeds[ia])) do is
            _interp1(t.coords[ia][is], t.values[ia][is], ℓ, l_axis)
        end
        _interp1(t.speeds[ia], per_speed, n, n_axis)
    end
    _interp1(t.alphas, per_alpha, α, a_axis)
end

# ── Reynolds correction ───────────────────────────────────────────────────────

"""
Reynolds-correction tables (the NPSS `S_Re` subelement).  `coord` is the
*abstract* table index — NPSS has no built-in meaning for it; the user chose
what to feed it (these files use Re/Re_des, but raw Re works equally well).
`s_eff` and `s_flow` multiply efficiency and corrected flow respectively.
"""
struct ReynoldsTables
    coord::Vector{Float64}
    s_eff::Vector{Float64}
    s_flow::Vector{Float64}
    axis::MapAxis
end

# (s_eff, s_flow) factors at table coordinate `rc`; (1,1) when absent or rc===nothing.
_re_factors(::Nothing, rc) = (1.0, 1.0)
_re_factors(::ReynoldsTables, ::Nothing) = (1.0, 1.0)
function _re_factors(re::ReynoldsTables, rc)
    (_interp1(re.coord, re.s_eff, rc, re.axis),
     _interp1(re.coord, re.s_flow, rc, re.axis))
end

"""
Element-side Reynolds model: turns an inlet `FluidState` into the scalar `rc`
fed to the map's Reynolds tables.  This is where the (user-defined, machine-type
dependent) convention lives — NPSS leaves it to the user, so we keep it pluggable
and default it to *off* (`reynolds = nothing` on the element ⇒ no correction).

  ReDesIndex(Re_des; convention)  rc = Re(state)/Re_des   (RNI-style index)
  RawRe(; convention)             rc = Re(state)          (table indexed on raw Re)
  FunctionReynolds(f)             rc = f(state)           (full user control)

`Re(state)` is a geometry-free Reynolds index `sqrt(γ·Pt·ρ)/μ` (a Reynolds number
per unit length; the absolute length cancels in the `ReDesIndex` ratio).  The
`convention` symbol (`:compressor`/`:turbine`) is the standardization hook for
machine-specific definitions; both currently use this property group.
"""
abstract type ReynoldsModel end

struct ReDesIndex <: ReynoldsModel
    Re_des::Float64
    convention::Symbol
end
ReDesIndex(Re_des; convention::Symbol = :compressor) = ReDesIndex(Re_des, convention)

struct RawRe <: ReynoldsModel
    convention::Symbol
end
RawRe(; convention::Symbol = :compressor) = RawRe(convention)

struct FunctionReynolds{F} <: ReynoldsModel
    f::F
end

# Geometry-free Reynolds index per unit length: sqrt(γ Pt ρ)/μ.
function reynolds_index(s)
    ρ = density(s)
    μ = viscosity(s.fluid, s.Tt, s.Pt)
    sqrt(gamma(s) * s.Pt * ρ) / μ
end

re_coord(m::ReDesIndex, s) = reynolds_index(s) / m.Re_des
re_coord(m::RawRe, s)      = reynolds_index(s)
re_coord(m::FunctionReynolds, s) = m.f(s)
re_coord(::Nothing, s)     = nothing

# ── Map types ─────────────────────────────────────────────────────────────────

abstract type TurbomachineMap end

"""
R-line parameterized compressor map: native `(alpha, Nc, Rline)`, outputs
`(Wc, PR, eff)`.  Scale factors `s_*` are 1 on a freshly loaded map and set by
`scale_map` to pass through a design point.  `alpha` is held at `alphaMapDes`
(variable-vane support would pass it through).
"""
struct CompressorMap{T<:Real} <: TurbomachineMap
    flow::MapTable          # Wc(alpha, Nc, Rline)
    pr::MapTable            # PR
    eff::MapTable           # eff
    s_Nc::T
    s_Wc::T
    s_PR::T
    s_eff::T
    NcMapDes::Float64
    RlineMapDes::Float64
    alphaMapDes::Float64
    RlineStall::Float64
    re::Union{ReynoldsTables,Nothing}
end

function eval_map(m::CompressorMap, Nc, Rline, rc = nothing)
    α     = m.alphaMapDes
    NcMap = Nc / m.s_Nc
    Wc_b  = eval_table(m.flow, α, NcMap, Rline)
    PR_b  = eval_table(m.pr,   α, NcMap, Rline)
    eff_b = eval_table(m.eff,  α, NcMap, Rline)
    sf_eff, sf_flow = _re_factors(m.re, rc)
    (; Wc  = m.s_Wc * Wc_b * sf_flow,
       PR  = 1 + m.s_PR * (PR_b - 1),
       eff = m.s_eff * eff_b * sf_eff)
end

"""
PR parameterized turbine map: native `(Np, PR)`, outputs `(Wp, eff)` — PR is an
*input* (set by the cycle), flow is an output.
"""
struct TurbineMap{T<:Real} <: TurbomachineMap
    flow::MapTable          # Wp(Np, PR)
    eff::MapTable           # eff
    s_Np::T
    s_Wp::T
    s_PR::T
    s_eff::T
    NpMapDes::Float64
    PRmapDes::Float64
    re::Union{ReynoldsTables,Nothing}
end

function eval_map(m::TurbineMap, Np, PR, rc = nothing)
    NpMap = Np / m.s_Np
    PRmap = 1 + (PR - 1) / m.s_PR     # physical PR → map PR coordinate
    Wp_b  = eval_table(m.flow, 0.0, NpMap, PRmap)
    eff_b = eval_table(m.eff,  0.0, NpMap, PRmap)
    sf_eff, sf_flow = _re_factors(m.re, rc)
    (; Wp  = m.s_Wp * Wp_b * sf_flow,
       eff = m.s_eff * eff_b * sf_eff)
end

"""
"Map as a script" backend: wraps a closure `f(speed, line, rc) -> NamedTuple`.
The closure returns whatever fields the consuming element reads
(`(; Wc, PR, eff)` for a compressor, `(; Wp, eff)` for a turbine).  `line_des`
is the line-coordinate seed a compressor uses on its first off-design pass.
"""
struct FunctionMap{F} <: TurbomachineMap
    f::F
    line_des::Float64
end
FunctionMap(f; line_des::Real = 2.0) = FunctionMap{typeof(f)}(f, Float64(line_des))
eval_map(m::FunctionMap, speed, line, rc = nothing) = m.f(speed, line, rc)

"Line-coordinate seed for a compressor's first off-design pass (R-line for an R-line map)."
design_line(m::CompressorMap) = m.RlineMapDes
design_line(m::FunctionMap)   = m.line_des

# ── Native-coordinate scaling ─────────────────────────────────────────────────

"""
    scale_map(base::CompressorMap; Nc_des, Wc_des, PR_des, eta_des) -> CompressorMap
    scale_map(base::TurbineMap;    Np_des, Wp_des, PR_des, eta_des) -> TurbineMap

Solve the four map scale factors so the scaled map passes exactly through the
design point at the map's own design anchors (`NcMapDes`/`RlineMapDes` for a
compressor, `NpMapDes`/`PRmapDes` for a turbine):

  s_Nc = Nc_des / NcMapDes,   s_Wc = Wc_des / Wc(anchor)
  s_PR = (PR_des - 1)/(PR(anchor) - 1),   s_eff = eta_des / eff(anchor)
"""
function scale_map(base::CompressorMap; Nc_des, Wc_des, PR_des, eta_des)
    α, n, ℓ = base.alphaMapDes, base.NcMapDes, base.RlineMapDes
    Wc_d  = eval_table(base.flow, α, n, ℓ)
    PR_d  = eval_table(base.pr,   α, n, ℓ)
    eff_d = eval_table(base.eff,  α, n, ℓ)
    s_Nc  = Nc_des / base.NcMapDes
    s_Wc  = Wc_des / Wc_d
    s_PR  = (PR_des - 1) / (PR_d - 1)
    s_eff = eta_des / eff_d
    T = promote_type(typeof(s_Nc), typeof(s_Wc), typeof(s_PR), typeof(s_eff))
    CompressorMap{T}(base.flow, base.pr, base.eff,
                     T(s_Nc), T(s_Wc), T(s_PR), T(s_eff),
                     base.NcMapDes, base.RlineMapDes, base.alphaMapDes,
                     base.RlineStall, base.re)
end

# ── Plotting support: physical speed-line sweeps in native coordinates ─────────

"""
    speed_lines(m; npts) -> Vector{(speed, x, PR)}

Sweep each tabulated speed line over its native line coordinate and return the
physical operating curve: corrected flow `x` (Wc/Wp) vs pressure ratio `PR`.
Used by `mapplot`.
"""
function speed_lines(m::CompressorMap; npts::Int = 40)
    map(eachindex(m.flow.speeds[1])) do is
        Nc = m.flow.speeds[1][is] * m.s_Nc
        rmin, rmax = extrema(m.flow.coords[1][is])
        Wc = Float64[]; PR = Float64[]
        for r in range(rmin, rmax, length = npts)
            o = eval_map(m, Nc, r)
            push!(Wc, o.Wc); push!(PR, o.PR)
        end
        (Nc, Wc, PR)
    end
end

function speed_lines(m::TurbineMap; npts::Int = 40)
    map(eachindex(m.flow.speeds[1])) do is
        Np = m.flow.speeds[1][is] * m.s_Np
        pmin, pmax = extrema(m.flow.coords[1][is])
        Wp = Float64[]; PR = Float64[]
        for pm in range(pmin, pmax, length = npts)
            p = 1 + m.s_PR * (pm - 1)        # map PR coord → physical PR
            o = eval_map(m, Np, p)
            push!(Wp, o.Wp); push!(PR, p)
        end
        (Np, Wp, PR)
    end
end

function scale_map(base::TurbineMap; Np_des, Wp_des, PR_des, eta_des)
    n, p = base.NpMapDes, base.PRmapDes
    Wp_d  = eval_table(base.flow, 0.0, n, p)
    eff_d = eval_table(base.eff,  0.0, n, p)
    s_Np  = Np_des / base.NpMapDes
    s_Wp  = Wp_des / Wp_d
    s_PR  = (PR_des - 1) / (base.PRmapDes - 1)
    s_eff = eta_des / eff_d
    T = promote_type(typeof(s_Np), typeof(s_Wp), typeof(s_PR), typeof(s_eff))
    TurbineMap{T}(base.flow, base.eff,
                  T(s_Np), T(s_Wp), T(s_PR), T(s_eff),
                  base.NpMapDes, base.PRmapDes, base.re)
end
