"""
Reader for NPSS turbomachinery map files (NEO `Subelement`/`Table` syntax).

Real NPSS map files wrap their tables in a subelement and may nest a Reynolds
subelement, e.g.

    Subelement CompressorRlineMap S_map {
        RlineMapDes = 2.212;  NcMapDes = 1.000;        // scalar design anchors
        Table TB_Wc(real alphaMap, real NcMap, real RlineMap) {
            alphaMap = 0.0 { NcMap = 0.5 { RlineMap = {…}  WcMap = {…} } … }
            NcMap.interp = "lagrange3";  RlineMap.extrap = "linear";  …
        }
        Table TB_PR(…) {…}   Table TB_eff(…) {…}
        Subelement CompressorReynoldsEffects S_Re {
            Table s_effRe(real RNI) { RNI = {…}  s_effRe = {…} }
            Table s_WcRe (real RNI) { RNI = {…}  s_WcRe  = {…} }
        }
    }

`read_npss_map(path)` parses this into a `ParsedSubelement` tree (subelements,
tables, scalars).  NPSS only cares about *hierarchy*, not names — within a table
the coordinate is the array named by the table's **last argument** and the value
is the other array — so the parser is fully positional and works with any axis
or leaf names (`ALPHA`/`SPED`/`R`/`FLOW`, `NcMap`/`WcMap`, …).

`compressor_map(path)` / `turbine_map(path)` build the scaled-once-removed
`CompressorMap`/`TurbineMap` the elements consume (call `scale_map` to place the
design point).  Table roles default to the common names but are overridable.
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
_scalar_value(t::AbstractString) = _is_number(t) ? parse(Float64, t) : String(t)

# ── Parse tree ────────────────────────────────────────────────────────────────

"""One `var = value { … }` level inside a table: numeric children or leaf arrays."""
struct _MapNode
    var::String
    val::Float64
    children::Vector{_MapNode}
    arrays::Vector{Pair{String,Vector{Float64}}}
end

"""A parsed `Table` block (data left as the raw nested-node tree + settings)."""
struct ParsedTable
    name::String
    argnames::Vector{String}
    root_children::Vector{_MapNode}                    # speed/alpha nesting
    root_arrays::Vector{Pair{String,Vector{Float64}}}  # 1-D tables: arrays directly
    settings::Dict{String,Any}                         # axis.interp, extrap flags, …
end

"""A parsed `Subelement` block (or the synthetic file root)."""
struct ParsedSubelement
    stype::String
    name::String
    scalars::Dict{String,Any}
    tables::Vector{ParsedTable}
    subelements::Vector{ParsedSubelement}
end

# ── Recursive-descent parser ──────────────────────────────────────────────────

# Parse a table body starting at toks[i] == "{"; returns (children, arrays, settings, next_i).
function _parse_table_body(toks::Vector{String}, i::Int)
    toks[i] == "{" || error("NPSS map: expected '{' at token $i, got '$(toks[i])'")
    i += 1
    children = _MapNode[]
    arrays   = Pair{String,Vector{Float64}}[]
    settings = Dict{String,Any}()

    while toks[i] != "}"
        name = toks[i]
        toks[i+1] == "=" || error("NPSS map: expected '=' after '$name'")
        i += 2
        if toks[i] == "{"                                   # leaf array
            i += 1
            vals = Float64[]
            while toks[i] != "}"
                toks[i] == "," || push!(vals, parse(Float64, toks[i]))
                i += 1
            end
            i += 1
            push!(arrays, name => vals)
        elseif _is_number(toks[i]) && toks[i+1] == "{"      # nested level
            v = parse(Float64, toks[i]); i += 1
            sub_c, sub_a, _, i = _parse_table_body(toks, i)
            push!(children, _MapNode(name, v, sub_c, sub_a))
        else                                                # scalar setting (interp/extrap/flag)
            settings[name] = _scalar_value(toks[i])
            i += 1
        end
    end
    (children, arrays, settings, i + 1)
end

# Parse a subelement body.  `is_root` parses to end-of-tokens; otherwise toks[i]
# is the opening "{" and parsing stops after the matching "}".
function _parse_subelement_body(toks::Vector{String}, i::Int, is_root::Bool)
    is_root || (i += 1)                                     # past "{"
    scalars = Dict{String,Any}()
    tables  = ParsedTable[]
    subs    = ParsedSubelement[]

    while i <= length(toks) && (is_root || toks[i] != "}")
        t = toks[i]
        if t == "Subelement"
            stype, sname = toks[i+1], toks[i+2]
            i += 3
            sc, tb, sb, i = _parse_subelement_body(toks, i, false)
            push!(subs, ParsedSubelement(stype, sname, sc, tb, sb))
        elseif t == "Table"
            tname = toks[i+1]
            toks[i+2] == "(" || error("NPSS map: expected '(' after Table $tname")
            i += 3
            argn = String[]
            while toks[i] != ")"
                tk = toks[i]
                (tk == "," || tk == "real" || tk == "int") || push!(argn, tk)
                i += 1
            end
            i += 1                                          # past ')'
            ch, arr, set, i = _parse_table_body(toks, i)
            push!(tables, ParsedTable(tname, argn, ch, arr, set))
        elseif i + 1 <= length(toks) && toks[i+1] == "="    # scalar: name = value
            scalars[t] = _scalar_value(toks[i+2])
            i += 3
        else
            i += 1                                          # stray token, skip
        end
    end
    is_root || (i += 1)                                     # past '}'
    (scalars, tables, subs, i)
end

"""
    read_npss_map(path) -> ParsedSubelement

Parse an NPSS map file into a subelement tree.  The returned node is a synthetic
file root whose `subelements`/`tables`/`scalars` are the file's top-level items.
"""
function read_npss_map(path::AbstractString)
    toks = _map_tokens(read(path, String))
    scalars, tables, subs, _ = _parse_subelement_body(toks, 1, true)
    (isempty(tables) && isempty(subs)) &&
        error("read_npss_map: no Table or Subelement blocks found in $path")
    ParsedSubelement("", "__root__", scalars, tables, subs)
end

# ── Tree navigation ───────────────────────────────────────────────────────────

_table_by_name(sub::ParsedSubelement, name) =
    (for t in sub.tables; t.name == name && return t; end; nothing)

"Depth-first search including `sub` itself."
function _find_subelement(sub::ParsedSubelement, pred)
    pred(sub) && return sub
    for s in sub.subelements
        r = _find_subelement(s, pred); r === nothing || return r
    end
    nothing
end

"Depth-first search of descendants only (not `sub`)."
function _find_descendant(sub::ParsedSubelement, pred)
    for s in sub.subelements
        pred(s) && return s
        r = _find_descendant(s, pred); r === nothing || return r
    end
    nothing
end

_getf(d::Dict, key, default) = haskey(d, key) ? Float64(d[key]) : Float64(default)

# ── ParsedTable → MapTable ────────────────────────────────────────────────────

function _coord_value(arrays, coord_name::AbstractString, tname)
    ci = findfirst(p -> first(p) == coord_name, arrays)
    isnothing(ci) && error("NPSS map: $tname: coordinate array '$coord_name' missing")
    vi = findfirst(p -> first(p) != coord_name, arrays)
    isnothing(vi) && error("NPSS map: $tname: value array missing")
    coord, vals = last(arrays[ci]), last(arrays[vi])
    length(coord) == length(vals) ||
        error("NPSS map: $tname: coordinate/value length mismatch " *
              "($(length(coord)) vs $(length(vals)))")
    (coord, vals)
end

function _axes_from_settings(argn::Vector{String}, settings::Dict)
    ax(name; di, de) = MapAxis(
        Symbol(get(settings, "$name.interp", String(di))),
        Symbol(get(settings, "$name.extrap", String(de))))
    if length(argn) >= 3
        (ax(argn[1]; di=:linear, de=:linear), ax(argn[2]; di=:lagrange3, de=:linear),
         ax(argn[3]; di=:lagrange3, de=:linear))
    else
        (MapAxis(:linear, :linear), ax(argn[1]; di=:lagrange3, de=:linear),
         ax(argn[2]; di=:lagrange3, de=:linear))
    end
end

function _build_maptable(pt::ParsedTable)
    coord_name = pt.argnames[end]
    blocks = length(pt.argnames) >= 3 ?
        [(nd.val, nd.children) for nd in pt.root_children] :
        [(0.0, pt.root_children)]

    alphas = Float64[]
    speeds = Vector{Float64}[]
    coords = Vector{Vector{Float64}}[]
    values = Vector{Vector{Float64}}[]
    for (aval, speed_nodes) in blocks
        isempty(speed_nodes) &&
            error("NPSS map: $(pt.name): no speed lines (is this a 1-D table?)")
        push!(alphas, aval)
        sp, co, va = Float64[], Vector{Float64}[], Vector{Float64}[]
        for snode in speed_nodes
            c, v = _coord_value(snode.arrays, coord_name, pt.name)
            push!(sp, snode.val); push!(co, c); push!(va, v)
        end
        push!(speeds, sp); push!(coords, co); push!(values, va)
    end
    MapTable(alphas, speeds, coords, values, _axes_from_settings(pt.argnames, pt.settings))
end

# ── Reynolds subelement → ReynoldsTables ──────────────────────────────────────

_has_re_tables(s::ParsedSubelement) = any(startswith(t.name, "s_") for t in s.tables)

function _onedim(pt::ParsedTable)
    coord_name = pt.argnames[end]
    _coord_value(pt.root_arrays, coord_name, pt.name)
end

function _build_reynolds(sub::Union{ParsedSubelement,Nothing})
    sub === nothing && return nothing
    eff_i  = findfirst(t -> occursin("eff", lowercase(t.name)), sub.tables)
    flow_i = findfirst(t -> !occursin("eff", lowercase(t.name)) &&
                            startswith(t.name, "s_"), sub.tables)
    (eff_i === nothing || flow_i === nothing) && return nothing
    coord, s_eff  = _onedim(sub.tables[eff_i])
    _,     s_flow = _onedim(sub.tables[flow_i])
    ReynoldsTables(coord, s_eff, s_flow, MapAxis(:linear, :linear))
end

# ── Public builders ───────────────────────────────────────────────────────────

"""
    compressor_map(path; flow="TB_Wc", pr="TB_PR", eff="TB_eff") -> CompressorMap

Load an NPSS R-line compressor map.  Returns an **unscaled** `CompressorMap`
(scale factors = 1) carrying the file's design anchors and Reynolds tables; call
[`scale_map`](@ref) to place the cycle design point.  `flow`/`pr`/`eff` name the
three tables and are overridable for files using other names.
"""
function compressor_map(path::AbstractString;
                        flow::AbstractString = "TB_Wc",
                        pr::AbstractString   = "TB_PR",
                        eff::AbstractString  = "TB_eff")
    root   = read_npss_map(path)
    mapsub = _find_subelement(root, s -> _table_by_name(s, flow) !== nothing)
    mapsub === nothing &&
        error("compressor_map: no subelement contains table '$flow' in $path")
    for nm in (pr, eff)
        _table_by_name(mapsub, nm) === nothing &&
            error("compressor_map: table '$nm' not found alongside '$flow'")
    end
    sc = mapsub.scalars
    re = _build_reynolds(_find_descendant(mapsub, _has_re_tables))
    CompressorMap{Float64}(
        _build_maptable(_table_by_name(mapsub, flow)),
        _build_maptable(_table_by_name(mapsub, pr)),
        _build_maptable(_table_by_name(mapsub, eff)),
        1.0, 1.0, 1.0, 1.0,
        _getf(sc, "NcMapDes", 1.0), _getf(sc, "RlineMapDes", 2.0),
        _getf(sc, "alphaMapDes", 0.0), _getf(sc, "RlineStall", 1.0), re)
end

"""
    turbine_map(path; flow="TB_Wp", eff="TB_eff") -> TurbineMap

Load an NPSS PR-parameterized turbine map (flow `Wp` and efficiency are outputs;
PR is an input).  Returns an **unscaled** `TurbineMap`; call [`scale_map`](@ref)
to place the design point.
"""
function turbine_map(path::AbstractString;
                     flow::AbstractString = "TB_Wp",
                     eff::AbstractString  = "TB_eff")
    root   = read_npss_map(path)
    mapsub = _find_subelement(root, s -> _table_by_name(s, flow) !== nothing)
    mapsub === nothing &&
        error("turbine_map: no subelement contains table '$flow' in $path")
    _table_by_name(mapsub, eff) === nothing &&
        error("turbine_map: table '$eff' not found alongside '$flow'")
    sc = mapsub.scalars
    re = _build_reynolds(_find_descendant(mapsub, _has_re_tables))
    TurbineMap{Float64}(
        _build_maptable(_table_by_name(mapsub, flow)),
        _build_maptable(_table_by_name(mapsub, eff)),
        1.0, 1.0, 1.0, 1.0,
        _getf(sc, "NpMapDes", 1.0), _getf(sc, "PRmapDes", 1.5), re)
end
