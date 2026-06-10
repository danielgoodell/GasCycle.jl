"""
Reader for NPSS turbomachinery performance map files (NEO `Table` syntax).

A standard NPSS compressor map file contains three-level nested tables:

    Table TB_Wc(real ALPHA, real NcMap, real RlineMap) {
        ALPHA = 0.0 {
            NcMap = 0.60 {
                RlineMap = { 1.0, 1.5, 2.0, 2.5, 3.0 }
                Wc       = { 0.62, 0.58, 0.54, 0.49, 0.44 }
            }
            NcMap = 0.70 { ... }
        }
    }
    Table TB_PR(...)  { ... }
    Table TB_eff(...) { ... }

`read_npss_map(path)` parses every table generically (any names, 2- or
3-level nesting, value lists wrapped across lines, `//` and `#` comments).

`to_performance_map(tables)` converts the standard compressor-map form to the
rectangular `(Nc, Wc)` `PerformanceMap` used by the solver, by inverting each
speed line's Wc(Rline) onto a common corrected-flow grid.

!!! warning
    Near choke a speed line is nearly vertical (Wc constant while PR varies),
    so PR(Wc) is ill-conditioned there and the rectangular resampling clamps
    to the line ends.  An Rline-parameterized solver coordinate (the native
    NPSS formulation) is the eventual fix; this converter covers maps whose
    operating range stays off the vertical segment.
"""

# ── Tokenizer ─────────────────────────────────────────────────────────────────

const _MAP_TOKEN_RE =
    r"//[^\n]*|#[^\n]*|[A-Za-z_][A-Za-z_0-9.]*|[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|[={}(),;]"

function _map_tokens(content::AbstractString)
    toks = String[]
    for m in eachmatch(_MAP_TOKEN_RE, content)
        t = m.match
        (startswith(t, "//") || startswith(t, "#")) && continue
        t == ";" && continue
        push!(toks, t)
    end
    toks
end

_is_number(t::AbstractString) = !isnothing(tryparse(Float64, t))

# ── Recursive-descent block parser ────────────────────────────────────────────

"""One `var = value { ... }` level: numeric children or leaf arrays."""
struct _MapNode
    var::String
    val::Float64
    children::Vector{_MapNode}
    arrays::Vector{Pair{String,Vector{Float64}}}
end

# Parses the region between '{' and '}' starting at toks[i] == "{".
# Returns (children, arrays, next_index).
function _parse_block(toks::Vector{String}, i::Int)
    toks[i] == "{" || error("NPSS map parse: expected '{' at token $i, got '$(toks[i])'")
    i += 1
    children = _MapNode[]
    arrays   = Pair{String,Vector{Float64}}[]

    while toks[i] != "}"
        name = toks[i]
        toks[i+1] == "=" || error("NPSS map parse: expected '=' after '$name'")
        i += 2
        if toks[i] == "{"                       # leaf array: name = { v, v, … }
            i += 1
            vals = Float64[]
            while toks[i] != "}"
                t = toks[i]
                t == "," || push!(vals, parse(Float64, t))
                i += 1
            end
            i += 1
            push!(arrays, name => vals)
        elseif _is_number(toks[i])              # nested: name = number { … }
            v = parse(Float64, toks[i])
            i += 1
            sub_children, sub_arrays, i = _parse_block(toks, i)
            push!(children, _MapNode(name, v, sub_children, sub_arrays))
        else
            error("NPSS map parse: unexpected token '$(toks[i])' after '$name ='")
        end
    end
    (children, arrays, i + 1)
end

# ── Table structure ───────────────────────────────────────────────────────────

"""
Parsed NPSS map table, normalized to three levels.

  argnames                   table arguments, e.g. ["ALPHA","NcMap","RlineMap"]
  alphas                     outer-variable values (a single 0.0 for 2-level tables)
  speeds[a]                  speed values for alpha index a
  coords[a][n]               inner coordinate vector (Rline, PR, …)
  values[a][n]               table quantity at each coordinate point
"""
struct NPSSMapTable
    name::String
    argnames::Vector{String}
    quantity::String                          # leaf value-array name, e.g. "Wc"
    alphas::Vector{Float64}
    speeds::Vector{Vector{Float64}}
    coords::Vector{Vector{Vector{Float64}}}
    values::Vector{Vector{Vector{Float64}}}
end

function _leaf_to_arrays(node::_MapNode, coord_name::String, tname::String)
    isempty(node.arrays) &&
        error("NPSS map: $tname: no data arrays under $(node.var) = $(node.val)")
    ci = findfirst(p -> first(p) == coord_name, node.arrays)
    isnothing(ci) &&
        error("NPSS map: $tname: coordinate array '$coord_name' missing under " *
              "$(node.var) = $(node.val)")
    vi = findfirst(p -> first(p) != coord_name, node.arrays)
    isnothing(vi) &&
        error("NPSS map: $tname: value array missing under $(node.var) = $(node.val)")
    coord = last(node.arrays[ci])
    vals  = last(node.arrays[vi])
    length(coord) == length(vals) ||
        error("NPSS map: $tname: coordinate/value length mismatch under " *
              "$(node.var) = $(node.val) ($(length(coord)) vs $(length(vals)))")
    (first(node.arrays[vi]), coord, vals)
end

function _build_table(name::String, argnames::Vector{String},
                      children::Vector{_MapNode},
                      arrays::Vector{Pair{String,Vector{Float64}}})
    coord_name = argnames[end]
    quantity   = ""

    # Normalize 2-level tables (no ALPHA) to a single alpha = 0.0
    alpha_nodes = if length(argnames) >= 3
        children
    else
        [_MapNode(length(argnames) >= 1 ? argnames[1] : "ALPHA", 0.0,
                  children, arrays)]
    end

    alphas = Float64[]
    speeds = Vector{Float64}[]
    coords = Vector{Vector{Float64}}[]
    values = Vector{Vector{Float64}}[]

    for a in alpha_nodes
        push!(alphas, a.val)
        sp = Float64[]
        co = Vector{Float64}[]
        va = Vector{Float64}[]
        for n in a.children
            q, c, v = _leaf_to_arrays(n, coord_name, name)
            quantity = q
            push!(sp, n.val); push!(co, c); push!(va, v)
        end
        push!(speeds, sp); push!(coords, co); push!(values, va)
    end

    NPSSMapTable(name, argnames, quantity, alphas, speeds, coords, values)
end

# ── Public API: reader ────────────────────────────────────────────────────────

"""
    read_npss_map(path) -> Dict{String,NPSSMapTable}

Parse an NPSS performance-map file.  Returns one entry per `Table` block,
keyed by table name (`"TB_Wc"`, `"TB_PR"`, `"TB_eff"`, …).
"""
function read_npss_map(path::String)
    toks = _map_tokens(read(path, String))
    tables = Dict{String,NPSSMapTable}()

    i = 1
    while i <= length(toks)
        if toks[i] == "Table"
            tname = toks[i+1]
            toks[i+2] == "(" || error("NPSS map parse: expected '(' after Table $tname")
            i += 3
            argnames = String[]
            while toks[i] != ")"
                t = toks[i]
                (t == "," || t == "real" || t == "int") || push!(argnames, t)
                i += 1
            end
            i += 1   # past ')'
            children, arrays, i = _parse_block(toks, i)
            tables[tname] = _build_table(tname, argnames, children, arrays)
        else
            i += 1
        end
    end

    isempty(tables) && error("read_npss_map: no Table blocks found in $path")
    tables
end

# ── Conversion to PerformanceMap ──────────────────────────────────────────────

"""Linear interpolation with end clamping; x need not be sorted."""
function _interp1_clamped(x::Vector{Float64}, y::Vector{Float64}, xq::Float64)
    p = sortperm(x)
    xs, ys = x[p], y[p]
    xq <= xs[1]   && return ys[1]
    xq >= xs[end] && return ys[end]
    k = searchsortedlast(xs, xq)
    dx = xs[k+1] - xs[k]
    dx ≈ 0.0 && return ys[k]   # vertical (choked) segment
    ys[k] + (xq - xs[k]) / dx * (ys[k+1] - ys[k])
end

"""
    to_performance_map(tables; alpha=0.0, flow="TB_Wc", pr="TB_PR",
                       eff="TB_eff", nWc=25) -> PerformanceMap

Convert a parsed NPSS compressor map (Rline-parameterized speed lines) to the
rectangular `(Nc, Wc)` `PerformanceMap` the solver consumes.  For each speed
line, PR and η are reinterpolated from Rline onto a common corrected-flow
grid (`nWc` points spanning the map's full flow range); flow values outside a
given speed line's range clamp to that line's end point.

NPSS corrected flow is in lbm/s — pass the result through `scale_map` at the
design point (the standard workflow), which makes the absolute flow units
irrelevant.
"""
function to_performance_map(tables::Dict{String,NPSSMapTable};
                            alpha::Float64 = NaN,
                            flow::String = "TB_Wc",
                            pr::String   = "TB_PR",
                            eff::String  = "TB_eff",
                            nWc::Int     = 25)
    for t in (flow, pr, eff)
        haskey(tables, t) || error("to_performance_map: table '$t' not found; " *
                                   "available: $(join(sort(collect(keys(tables))), ", "))")
    end
    tWc, tPR, teff = tables[flow], tables[pr], tables[eff]

    ai = isnan(alpha) ? 1 : findfirst(≈(alpha), tWc.alphas)
    isnothing(ai) && error("to_performance_map: alpha=$alpha not in $(tWc.alphas)")

    Nc_axis = tWc.speeds[ai]
    (Nc_axis == tPR.speeds[ai] && Nc_axis == teff.speeds[ai]) ||
        error("to_performance_map: speed grids differ between tables")

    Wc_lines = tWc.values[ai]
    Wc_all   = reduce(vcat, Wc_lines)
    # Every source sample's flow value goes on the common axis, so each
    # tabulated speed line — piecewise linear in Wc — is represented exactly
    # (the clamp plateau beyond a line's range starts exactly at its end
    # point, and no grid cell straddles a source breakpoint).  The `nWc`
    # background grid only adds resolution for the cross-speed direction.
    Wc_axis = sort!(unique!(vcat(
        collect(range(minimum(Wc_all), maximum(Wc_all), length = nWc)), Wc_all)))

    nN = length(Nc_axis)
    PR_grid  = Matrix{Float64}(undef, nN, length(Wc_axis))
    eta_grid = Matrix{Float64}(undef, nN, length(Wc_axis))

    for i in 1:nN
        (tWc.coords[ai][i] == tPR.coords[ai][i] == teff.coords[ai][i]) ||
            error("to_performance_map: Rline grids differ between tables at " *
                  "Nc=$(Nc_axis[i])")
        Wc_line  = Wc_lines[i]
        PR_line  = tPR.values[ai][i]
        eff_line = teff.values[ai][i]
        for (j, wq) in enumerate(Wc_axis)
            PR_grid[i, j]  = _interp1_clamped(Wc_line, PR_line, wq)
            eta_grid[i, j] = _interp1_clamped(Wc_line, eff_line, wq)
        end
    end

    PerformanceMap(Nc_axis, Wc_axis, PR_grid, eta_grid)
end
