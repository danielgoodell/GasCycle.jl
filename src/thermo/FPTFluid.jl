"""
Fluid property backend that reads NPSS-format FPT (Fluid Property Table) files.

## File format

Standard NPSS thermodynamic property table files produced by NPSS's
`createFluidPropTable()` function (or equivalent). The format is:

    indeps   = { "Pt", "Tt", ""}
    ...
    Table h_T( real Pt, real Tt ) {
        Pt = <value_psia> {
            Tt  = { R1, R2, ... }   // Rankine
            ht  = { v1, v2, ... }   // BTU/lbm
        }
        ...
    }

## Unit conventions (NPSS ENGLISH)
  Pt  — psia
  Tt  — Rankine
  ht  — BTU/lbm
  s   — BTU/(lbm·R)
  Cp  — BTU/(lbm·R)
  rho — lbm/ft³
  gam — dimensionless

All values are converted to SI on load.

## Tables used
  h_T  → enthalpy(Pt, Tt)        (primary forward lookup)
  s_T  → entropy(Pt, Tt)
  Cp   → cp(Pt, Tt)
  gam  → gamma(Pt, Tt)
  rho  → density(Pt, Tt)
  T_h  → T_from_h(Pt, ht)        (pre-computed ragged inverse)
  T_s  → T_from_s(Pt, s)         (pre-computed ragged inverse)
  h_s  → h_from_s(Pt, s)         (isentropic: h at same entropy, different P)
"""

using Interpolations
using ForwardDiff: ForwardDiff

# ── Unit conversion factors (ENGLISH → SI) ───────────────────────────────────
const _PSIA_TO_PA  = 6894.757
const _R_TO_K      = 5.0 / 9.0
const _BTU_LBM_TO_J_KG   = 2326.0            # BTU/lbm → J/kg
const _BTU_LBM_R_TO_SI   = 4186.8            # BTU/(lbm·R) → J/(kg·K)
const _LBM_FT3_TO_KG_M3  = 16.01846          # lbm/ft³ → kg/m³

# ── Main struct ───────────────────────────────────────────────────────────────
struct FPTFluid{H,S,CP,GAM,RHO,TH,TS,HS} <: FluidProperties
    name::String
    bounds::Symbol  # :error, :warn, or :clamp for out-of-table T/P lookups
    s_interp::Symbol  # :log_pressure (detrended, default) or :linear (legacy/NPSS-compat)
    R_s::Float64      # gas constant fitted from the s table [J/(kg·K)] (:log_pressure)
    P_ref::Float64    # detrending reference pressure [Pa]
    itp_h::H        # h(Pt, Tt)  [J/kg]
    itp_s::S        # s(Pt, Tt) [J/(kg·K)], or detrended σ = s + R_s·ln(Pt/P_ref)
    itp_cp::CP      # cp(Pt, Tt) [J/(kg·K)]
    itp_gam::GAM    # γ(Pt, Tt)  [-]
    itp_rho::RHO    # ρ(Pt, Tt)  [kg/m³]
    # Inverse / isentropic tables (optional; may be nothing)
    itp_Th::TH      # Tt(Pt, ht) [K] from pre-computed T_h table
    itp_Ts::TS      # Tt(Pt, s)  [K] from pre-computed T_s table
    itp_hs::HS      # ht(Pt, s)  [J/kg] from h_s table
    Pt_min::Float64; Pt_max::Float64
    Tt_min::Float64; Tt_max::Float64
end

struct _RaggedFPTTable
    Pt_axis::Vector{Float64}
    x_rows::Vector{Vector{Float64}}
    y_rows::Vector{Vector{Float64}}
end

_fpt_value(x) = x
_fpt_value(x::ForwardDiff.Dual) = _fpt_value(ForwardDiff.value(x))

# ── Fast line-by-line FPT parser ─────────────────────────────────────────────

"""
Parse a `{ v1, v2, ..., vN }` value list from a single line fragment.
Returns a pre-allocated Float64 vector, avoiding intermediate string splits.
"""
function _parse_inline_values(s::AbstractString, buf::Vector{Float64})
    empty!(buf)
    i = firstindex(s)
    n = lastindex(s)
    while i <= n
        c = s[i]
        (c == ' ' || c == '\t' || c == ',' || c == '{' || c == '}') && (i = nextind(s, i); continue)
        # start of a number token
        j = i
        while j <= n
            nc = s[j]
            (nc == ',' || nc == ' ' || nc == '\t' || nc == '}') && break
            j = nextind(s, j)
        end
        push!(buf, parse(Float64, @view s[i:prevind(s,j)]))
        i = j
    end
    copy(buf)
end

"""
Parse all Table blocks from FPT file content using a line-by-line state machine.
Returns Dict{String, Tuple{Vector{Float64}, Vector{Vector{Float64}}, Vector{Vector{Float64}}}}.
"""
function _parse_fpt_tables(content::AbstractString)
    tables  = Dict{String, Any}()
    scratch = Float64[]         # reusable parse buffer

    current_table = ""
    pt_vals   = Float64[]
    inner_arr = Vector{Float64}[]
    out_arr   = Vector{Float64}[]
    arrays_this_pt = 0
    inner_buf = Float64[]       # set by first key=value in each Pt block
    brace_depth = 0             # 0 = top level, 1 = inside Table, 2 = inside Pt block

    for raw_line in eachline(IOBuffer(content))
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line, "//") && continue

        # ── Table header: Table name( real Pt, real X ) { ──────────────────
        if startswith(line, "Table ")
            # flush any open table
            if !isempty(current_table) && !isempty(pt_vals)
                tables[current_table] = (copy(pt_vals), copy(inner_arr), copy(out_arr))
            end
            current_table = ""
            empty!(pt_vals); empty!(inner_arr); empty!(out_arr)
            brace_depth = 0

            m = match(r"^Table\s+(\w+)\s*\(", line)
            isnothing(m) && continue
            current_table = m.captures[1]
            endswith(rstrip(line), '{') && (brace_depth = 1)
            continue
        end

        isempty(current_table) && continue

        # Track opening/closing braces to know depth
        opens  = count(==( '{'), line)
        closes = count(==('}'), line)

        # ── Pt = value { ────────────────────────────────────────────────────
        if brace_depth == 1
            m = match(r"^Pt\s*=\s*([\d.eE+\-]+)", line)
            if !isnothing(m)
                push!(pt_vals, parse(Float64, m.captures[1]))
                arrays_this_pt = 0
                brace_depth = 2
                continue
            end
            # bare '{' opening the table block
            if line == "{"
                brace_depth = 1
                continue
            end
            # bare '}' closes the table
            if line == "}"
                brace_depth = 0
                continue
            end
        end

        # ── Inside a Pt block: key = { values } ─────────────────────────────
        if brace_depth == 2
            # key = { ... } on one line
            m = match(r"^\w+\s*=\s*(\{[^}]*\})", line)
            if !isnothing(m)
                vals = _parse_inline_values(m.captures[1], scratch)
                if arrays_this_pt == 0
                    inner_buf = vals
                    arrays_this_pt = 1
                else
                    push!(inner_arr, inner_buf)
                    push!(out_arr, vals)
                    arrays_this_pt = 0
                end
            end
            # closing brace for the Pt block
            if closes > opens
                brace_depth = 1
            end
            continue
        end

        # Top-level closing brace ends the Table block
        if brace_depth == 1 && line == "}"
            brace_depth = 0
        end
    end

    # flush last table
    if !isempty(current_table) && !isempty(pt_vals)
        tables[current_table] = (copy(pt_vals), copy(inner_arr), copy(out_arr))
    end

    tables
end

"""
Build a 2D Gridded interpolant from a parsed FPT table.
Assumes the inner variable grid is the same for every Pt level.
"""
function _build_itp(Pt_vals, inner_arr, out_arr,
                    Pt_fac::Float64, in_fac::Float64, out_fac::Float64)
    nPt = length(Pt_vals)
    nIn = length(inner_arr[1])
    data = Matrix{Float64}(undef, nPt, nIn)
    for i in 1:nPt
        length(out_arr[i]) == nIn ||
            error("FPTFluid: inconsistent row length at Pt index $i")
        @. data[i, :] = out_arr[i] * out_fac
    end
    Pt_axis = Pt_vals .* Pt_fac
    in_axis = inner_arr[1] .* in_fac
    interpolate((Pt_axis, in_axis), data, Gridded(Linear()))
end

function _build_ragged_itp(Pt_vals, inner_arr, out_arr,
                           Pt_fac::Float64, in_fac::Float64, out_fac::Float64;
                           in_shift = nothing)
    nPt = length(Pt_vals)
    Pt_axis = Pt_vals .* Pt_fac
    x_rows = Vector{Vector{Float64}}(undef, nPt)
    y_rows = Vector{Vector{Float64}}(undef, nPt)
    for i in 1:nPt
        length(inner_arr[i]) == length(out_arr[i]) ||
            error("FPTFluid: inconsistent row length at Pt index $i")
        x = inner_arr[i] .* in_fac
        if !isnothing(in_shift)
            @. x += in_shift(Pt_axis[i])
        end
        y = out_arr[i] .* out_fac
        if length(x) < 2
            error("FPTFluid: ragged inverse table needs at least two points per row")
        elseif x[1] > x[end]
            reverse!(x)
            reverse!(y)
        end
        issorted(x) || error("FPTFluid: ragged inverse table row $i is not monotone")
        x_rows[i] = x
        y_rows[i] = y
    end
    _RaggedFPTTable(Pt_axis, x_rows, y_rows)
end

function _bracket_index(axis::Vector{Float64}, x)
    xv = _fpt_value(x)
    xv <= axis[1] && return 1
    xv >= axis[end] && return length(axis) - 1
    clamp(searchsortedlast(axis, xv), 1, length(axis) - 1)
end

function _interp_row(x::Vector{Float64}, y::Vector{Float64}, target)
    tv = _fpt_value(target)
    tv <= x[1]   && return y[1]   + zero(target)
    tv >= x[end] && return y[end] + zero(target)
    j = clamp(searchsortedlast(x, tv), 1, length(x) - 1)
    w = (target - x[j]) / (x[j+1] - x[j])
    y[j] + w * (y[j+1] - y[j])
end

function _ragged_lookup(tbl::_RaggedFPTTable, P, x)
    i = _bracket_index(tbl.Pt_axis, P)
    y1 = _interp_row(tbl.x_rows[i],   tbl.y_rows[i],   x)
    y2 = _interp_row(tbl.x_rows[i+1], tbl.y_rows[i+1], x)
    wP = (P - tbl.Pt_axis[i]) / (tbl.Pt_axis[i+1] - tbl.Pt_axis[i])
    y1 + wP * (y2 - y1)
end

"""
    _fit_R_from_s(Pt_axis, s_data) -> R [J/(kg·K)]

Fit the specific gas constant from the s table itself: for an (ideal or
mildly real) gas, s(T, P) = f(T) − R·ln(P/P_ref), so every adjacent-Pt node
pair at fixed T gives R ≈ −Δs/Δln(Pt).  Take the median over all pairs —
robust to real-gas content in any corner of the table.
"""
function _fit_R_from_s(Pt_axis::Vector{Float64}, s_data::Matrix{Float64})
    slopes = Float64[]
    for i in 1:length(Pt_axis)-1, j in axes(s_data, 2)
        push!(slopes, -(s_data[i+1, j] - s_data[i, j]) / log(Pt_axis[i+1] / Pt_axis[i]))
    end
    sort!(slopes)
    slopes[(length(slopes) + 1) ÷ 2]
end

"""
    FPTFluid(path; bounds=:error, s_interp=:log_pressure) -> FPTFluid

Load a fluid property table from an NPSS FPT file.

`bounds` controls property lookups outside the table domain:
- `:error` throws an informative error immediately (default)
- `:warn` clamps to the table bounds and warns
- `:clamp` clamps silently, matching the legacy behavior

`s_interp` controls how entropy is interpolated between pressure nodes:
- `:log_pressure` (default) — interpolate the pressure-detrended entropy
  σ = s + R·ln(Pt/P_ref) and restore the −R·ln(Pt/P_ref) term analytically
  (R fitted from the table).  Exact in P for an ideal gas; removes the
  large mid-cell ∂s/∂P error of linear-in-P interpolation on coarse Pt
  grids (~14 % for HeXe84.fpt at the BRU compressor inlet, which showed up
  as +14 °R on the compressor outlet — see validation/PLAN.md, rung 0).
- `:linear` — legacy bilinear-in-P behavior, kept for apples-to-apples
  comparison with table interpolators that are linear in P (e.g. NPSS).
"""
function FPTFluid(path::String; bounds::Symbol = :error,
                  s_interp::Symbol = :log_pressure)
    bounds in (:error, :warn, :clamp) ||
        error("FPTFluid: bounds must be :error, :warn, or :clamp, got :$bounds")
    s_interp in (:log_pressure, :linear) ||
        error("FPTFluid: s_interp must be :log_pressure or :linear, got :$s_interp")

    content = read(path, String)
    name    = splitext(basename(path))[1]
    tables  = _parse_fpt_tables(content)

    isempty(tables) && error("FPTFluid: no Table blocks found in $path")

    function get_itp(tname, Pt_fac, in_fac, out_fac)
        haskey(tables, tname) || return nothing
        Pt_vals, inner_arr, out_arr = tables[tname]
        isempty(Pt_vals) && return nothing
        _build_itp(Pt_vals, inner_arr, out_arr, Pt_fac, in_fac, out_fac)
    end

    function get_ragged_itp(tname, Pt_fac, in_fac, out_fac; in_shift = nothing)
        haskey(tables, tname) || return nothing
        Pt_vals, inner_arr, out_arr = tables[tname]
        isempty(Pt_vals) && return nothing
        _build_ragged_itp(Pt_vals, inner_arr, out_arr, Pt_fac, in_fac, out_fac;
                          in_shift = in_shift)
    end

    itp_h   = get_itp("h_T", _PSIA_TO_PA, _R_TO_K, _BTU_LBM_TO_J_KG)
    itp_cp  = get_itp("Cp",  _PSIA_TO_PA, _R_TO_K, _BTU_LBM_R_TO_SI)
    itp_gam = get_itp("gam", _PSIA_TO_PA, _R_TO_K, 1.0)
    itp_rho = get_itp("rho", _PSIA_TO_PA, _R_TO_K, _LBM_FT3_TO_KG_M3)

    # Entropy: optionally interpolate the pressure-detrended σ = s + R·ln(Pt/P_ref)
    itp_s = nothing
    R_s   = 0.0
    P_ref = 1.0
    if haskey(tables, "s_T") && !isempty(tables["s_T"][1])
        Pt_vals, inner_arr, out_arr = tables["s_T"]
        nPt = length(Pt_vals)
        nIn = length(inner_arr[1])
        s_data = Matrix{Float64}(undef, nPt, nIn)
        for i in 1:nPt
            length(out_arr[i]) == nIn ||
                error("FPTFluid: inconsistent row length at Pt index $i")
            @. s_data[i, :] = out_arr[i] * _BTU_LBM_R_TO_SI
        end
        Pt_axis = Pt_vals .* _PSIA_TO_PA
        if s_interp == :log_pressure
            P_ref = Pt_axis[1]
            R_s   = _fit_R_from_s(Pt_axis, s_data)
            for i in 1:nPt
                @. s_data[i, :] += R_s * log(Pt_axis[i] / P_ref)
            end
        end
        itp_s = interpolate((Pt_axis, inner_arr[1] .* _R_TO_K), s_data, Gridded(Linear()))
    end

    isnothing(itp_h)  && error("FPTFluid: required table 'h_T' not found in $path")
    isnothing(itp_s)  && error("FPTFluid: required table 's_T' not found in $path")
    isnothing(itp_cp) && error("FPTFluid: table 'Cp' not found in $path")

    # Inverse tables are ragged: h/s nodes vary with pressure.  Use row-wise
    # interpolation in the inner coordinate, then interpolate the two row
    # results in pressure.
    s_shift = s_interp == :log_pressure ? P -> R_s * log(P / P_ref) : nothing
    itp_Th = get_ragged_itp("T_h", _PSIA_TO_PA, _BTU_LBM_TO_J_KG, _R_TO_K)
    itp_Ts = get_ragged_itp("T_s", _PSIA_TO_PA, _BTU_LBM_R_TO_SI, _R_TO_K;
                            in_shift = s_shift)
    itp_hs = get_ragged_itp("h_s", _PSIA_TO_PA, _BTU_LBM_R_TO_SI, _BTU_LBM_TO_J_KG;
                            in_shift = s_shift)

    Pt_vals_SI = tables["h_T"][1] .* _PSIA_TO_PA
    Tt_vals_SI = tables["h_T"][2][1] .* _R_TO_K

    FPTFluid(name, bounds, s_interp, R_s, P_ref,
             itp_h, itp_s, itp_cp, itp_gam, itp_rho,
             itp_Th, itp_Ts, itp_hs,
             minimum(Pt_vals_SI), maximum(Pt_vals_SI),
             minimum(Tt_vals_SI), maximum(Tt_vals_SI))
end

# ── FluidProperties interface ─────────────────────────────────────────────────

function _bounds_message(fp::FPTFluid, prop::Symbol, T, P)
    "FPTFluid $(fp.name): $(prop) requested outside table bounds. " *
    "Requested: T=$(T) K, P=$(P) Pa. " *
    "Valid: T=$(fp.Tt_min)..$(fp.Tt_max) K, P=$(fp.Pt_min)..$(fp.Pt_max) Pa. " *
    "Use FPTFluid(path; bounds=:warn) or bounds=:clamp to allow clamped lookup."
end

function _bounds_message_P(fp::FPTFluid, prop::Symbol, P)
    "FPTFluid $(fp.name): $(prop) requested outside table pressure bounds. " *
    "Requested: P=$(P) Pa. Valid: P=$(fp.Pt_min)..$(fp.Pt_max) Pa. " *
    "Use FPTFluid(path; bounds=:warn) or bounds=:clamp to allow clamped lookup."
end

function _clamp_PT(fp::FPTFluid, T, P, prop::Symbol)
    Pv = _fpt_value(P)
    Tv = _fpt_value(T)
    out = (Pv < fp.Pt_min || Pv > fp.Pt_max ||
           Tv < fp.Tt_min || Tv > fp.Tt_max)
    if out
        msg = _bounds_message(fp, prop, T, P)
        fp.bounds == :error && throw(DomainError((T=T, P=P), msg))
        fp.bounds == :warn && @warn msg
    end
    (clamp(P, fp.Pt_min, fp.Pt_max),
     clamp(T, fp.Tt_min, fp.Tt_max))
end

function _clamp_P(fp::FPTFluid, P, prop::Symbol)
    Pv = _fpt_value(P)
    out = Pv < fp.Pt_min || Pv > fp.Pt_max
    if out
        msg = _bounds_message_P(fp, prop, P)
        fp.bounds == :error && throw(DomainError(P, msg))
        fp.bounds == :warn && @warn msg
    end
    clamp(P, fp.Pt_min, fp.Pt_max)
end

function cp(fp::FPTFluid, T, P)
    Pc, Tc = _clamp_PT(fp, T, P, :cp)
    isnothing(fp.itp_cp) && return enthalpy(fp, T+0.5, P) - enthalpy(fp, T-0.5, P)
    fp.itp_cp(Pc, Tc)
end

function enthalpy(fp::FPTFluid, T, P)
    Pc, Tc = _clamp_PT(fp, T, P, :enthalpy)
    fp.itp_h(Pc, Tc)
end

function entropy(fp::FPTFluid, T, P)
    Pc, Tc = _clamp_PT(fp, T, P, :entropy)
    s = fp.itp_s(Pc, Tc)
    fp.s_interp == :log_pressure ? s - fp.R_s * log(Pc / fp.P_ref) : s
end

function density(fp::FPTFluid, T, P)
    isnothing(fp.itp_rho) && return P / (cp(fp,T,P) * (1.0 - 1.0/gamma(fp,T,P)) * T)
    Pc, Tc = _clamp_PT(fp, T, P, :density)
    fp.itp_rho(Pc, Tc)
end

function gamma(fp::FPTFluid, T, P)
    isnothing(fp.itp_gam) && return 5.0/3.0
    Pc, Tc = _clamp_PT(fp, T, P, :gamma)
    fp.itp_gam(Pc, Tc)
end

function T_from_h(fp::FPTFluid, h_target, P; T_guess=500.0)
    if !isnothing(fp.itp_Th)
        Pc = _clamp_P(fp, P, :T_from_h)
        return _ragged_lookup(fp.itp_Th, Pc, h_target)
    end
    # bisection fallback (returns primal; does not propagate AD derivatives)
    T_lo, T_hi = fp.Tt_min, fp.Tt_max
    for _ in 1:60
        T_mid = 0.5 * (T_lo + T_hi)
        enthalpy(fp, T_mid, P) < h_target ? (T_lo = T_mid) : (T_hi = T_mid)
        (T_hi - T_lo) < 1e-6 && break
    end
    0.5 * (T_lo + T_hi)
end

function T_from_s(fp::FPTFluid, s_target, P; T_guess=500.0)
    if !isnothing(fp.itp_Ts)
        Pc = _clamp_P(fp, P, :T_from_s)
        s_lookup = fp.s_interp == :log_pressure ?
                   s_target + fp.R_s * log(Pc / fp.P_ref) : s_target
        return _ragged_lookup(fp.itp_Ts, Pc, s_lookup)
    end
    # bisection fallback (returns primal; does not propagate AD derivatives)
    T_lo, T_hi = fp.Tt_min, fp.Tt_max
    for _ in 1:60
        T_mid = 0.5 * (T_lo + T_hi)
        entropy(fp, T_mid, P) < s_target ? (T_lo = T_mid) : (T_hi = T_mid)
        (T_hi - T_lo) < 1e-6 && break
    end
    0.5 * (T_lo + T_hi)
end

"""
    h_from_s(fp, s, Pt_out) -> h

Enthalpy at pressure Pt_out for an isentropic process from entropy s.
"""
function h_from_s(fp::FPTFluid, s, P)
    if !isnothing(fp.itp_hs)
        Pc = _clamp_P(fp, P, :h_from_s)
        s_lookup = fp.s_interp == :log_pressure ?
                   s + fp.R_s * log(Pc / fp.P_ref) : s
        return _ragged_lookup(fp.itp_hs, Pc, s_lookup)
    end
    T_exit = T_from_s(fp, s, P)
    enthalpy(fp, T_exit, P)
end
