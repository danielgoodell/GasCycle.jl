# ── Assemblies: named, nestable groups of elements ───────────────────────────
#
# An Assembly is built with the same calls as a FlowNetwork (it wraps one) and
# publishes connection points with `expose!`. `add!(net, asm)` merges it into
# the network: elements, shafts, electrical nodes, heat links, boundaries and
# solver pairs are copied over, and names are prefixed with the assembly name
# ("PCU1.TAC.Comp"). The solver never sees assemblies — flattening costs
# nothing at run time — but names, lookups and reports keep the grouping.

"""
    Assembly(name)

A named group of elements — a turbine-alternator-compressor, a power
conversion unit, a heat rejection loop — built like a network:

    tac = Assembly("TAC")
    add!(tac, alt, comp, turb)
    add_shaft!(tac, shaft; drives = [comp, alt], driven_by = turb)
    expose!(tac, :inlet, alt)                # gas in at the alternator
    expose!(tac, :comp_out, comp, :outlet)
    expose!(tac, :elec, alt)                 # electrical source

Every network-building call works on an assembly (`add!`, `connect!`,
`connect_port!`, `add_hx_pair!`, `add_shaft!`, `set_boundary!`,
`connect_electrical!`, `add_conductance!`/`add_radiation!`/`add_loss_heat!`,
`add_free!`/`add_target!`, `add_speed_control!`), including adding another
assembly (nesting). Assemblies go anywhere an element goes in `connect!`,
`connect_port!`, `add_hx_pair!`, `set_state!` and `connect_electrical!`: an
assembly in a chain connects through its exposed `:inlet` / `:outlet` (and
`:elec` for electrical chains). `asm["Comp"]` looks a part up by its local name.
"""
struct Assembly
    name::String
    net::FlowNetwork
    ports::Dict{Symbol,Tuple{Any,Symbol}}      # exposed fluid ports → (element or Assembly, port)
    eports::Dict{Symbol,Any}                   # exposed electrical nodes
    parts::Dict{String,Any}                    # local name → element / node / sub-assembly
end
Assembly(name::AbstractString) =
    Assembly(String(name), FlowNetwork(), Dict{Symbol,Tuple{Any,Symbol}}(), Dict{Symbol,Any}(),
             Dict{String,Any}())

Base.show(io::IO, a::Assembly) =
    print(io, "Assembly(\"$(a.name)\", $(length(a.parts)) parts, ports: ",
          join(sort(string.(collect(keys(a.ports)) ∪ collect(keys(a.eports)))), ", "), ")")

function Base.getindex(a::Assembly, name::AbstractString)
    haskey(a.parts, name) && return a.parts[name]
    # "TAC.Comp": walk into sub-assemblies
    parts = split(name, '.'; limit = 2)
    if length(parts) == 2 && haskey(a.parts, parts[1]) && a.parts[parts[1]] isa Assembly
        return a.parts[parts[1]][String(parts[2])]
    end
    error("assembly \"$(a.name)\" has no part \"$name\"")
end

_register!(a::Assembly, x) = (a.parts[x.name] = x; x)

# ── Exposed ports ────────────────────────────────────────────────────────────

"""
    expose!(asm, name, el, port = :outlet / :inlet)
    expose!(asm, name, node)          # electrical node (or alternator) for :elec

Publish a connection point of an assembly. Fluid: `name` maps to
`(el, port)`; when `port` is omitted it defaults to `name` if that is a port
name (`:inlet`, `:outlet`, `:hot_inlet` …). An electrical exposure maps
`name` to a node; `:elec` is the one used when the assembly appears in a
`connect_electrical!` chain. `el` may itself be a sub-assembly (its exposed
port `port`).
"""
function expose!(a::Assembly, name::Symbol, el, port::Union{Symbol,Nothing} = nothing)
    if isnothing(port) && !(el isa AbstractElement || el isa Assembly)
        a.eports[name] = el                  # electrical node
    elseif isnothing(port) && el isa Alternator && name === :elec
        a.eports[name] = el
    else
        a.ports[name] = (el, something(port, name))
    end
    a
end

# Resolve (thing, port) through assembly exposures down to (element, port).
function _resolve(x, port::Symbol)
    x isa Assembly || return (x, port)
    haskey(x.ports, port) ||
        error("assembly \"$(x.name)\" exposes no port :$port (exposed: $(join(keys(x.ports), ", ")))")
    el, p = x.ports[port]
    _resolve(el, p)
end

function _resolve_e(x, which::Symbol = :elec)
    x isa Assembly || return x
    haskey(x.eports, which) ||
        error("assembly \"$(x.name)\" exposes no electrical node :$which")
    _resolve_e(x.eports[which], which)
end

# ── Building calls on assemblies (forward to the wrapped network) ────────────

function add!(a::Assembly, items...)
    for x in items
        _register!(a, x)          # local name, before any prefixing
        add!(a.net, x)
    end
    a
end

add_shaft!(a::Assembly, shaft::Shaft; kw...) =
    (_register!(a, shaft); add_shaft!(a.net, shaft; kw...); a)

for f in (:set_boundary!, :add_conductance!, :add_radiation!, :add_volume!,
          :add_loss_heat!, :add_free!, :add_target!, :add_speed_control!)
    @eval $f(a::Assembly, args...; kw...) = ($f(a.net, args...; kw...); a)
end

function _register_e!(a::Assembly, chain)
    for x in _flatten_chain(chain)
        x isa Assembly || haskey(a.parts, x.name) || _register!(a, x)
    end
end

connect!(a::Assembly, chain)  = (connect!(a.net, chain); a)
connect_port!(a::Assembly, args...; kw...) = (connect_port!(a.net, args...; kw...); a)
add_hx_pair!(a::Assembly, hx; kw...) = (add_hx_pair!(a.net, hx; kw...); a)
connect_electrical!(a::Assembly, chain) = (_register_e!(a, chain); connect_electrical!(a.net, chain); a)

# ── Assembly-aware connection on networks ────────────────────────────────────
# These methods catch calls where an Assembly appears and resolve it.

const _AsmOrEl = Union{AbstractElement,Assembly}

function connect_port!(net::FlowNetwork, src::Assembly, src_port::Symbol,
                       dst::_AsmOrEl, dst_port::Symbol; back_edge::Bool = false)
    s, sp = _resolve(src, src_port); d, dp = _resolve(dst, dst_port)
    connect_port!(net, s, sp, d, dp; back_edge)
end
function connect_port!(net::FlowNetwork, src::AbstractElement, src_port::Symbol,
                       dst::Assembly, dst_port::Symbol; back_edge::Bool = false)
    d, dp = _resolve(dst, dst_port)
    connect_port!(net, src, src_port, d, dp; back_edge)
end

function set_state!(net::FlowNetwork, first::Assembly; kw...)
    el, port = _resolve(first, :inlet)
    port === :inlet || error("set_state!: the seed must be an element's :inlet (assembly \"$(first.name)\" exposes :inlet as $port)")
    set_state!(net, el; kw...)
end

# A chain containing assemblies: connect consecutive items through their
# exposed :outlet / :inlet (a closing repeat of the first item is dropped,
# as for element chains).
function _connect_chain_with_assemblies!(net::FlowNetwork, ordered)
    if length(ordered) > 1 && ordered[end] === ordered[1]
        ordered = ordered[1:end-1]
    end
    for (a, b) in zip(ordered[1:end-1], ordered[2:end])
        s, sp = _resolve(a, :outlet); d, dp = _resolve(b, :inlet)
        connect_port!(net, s, sp, d, dp)
    end
    net
end

# ── Merging an assembly into a network (or a parent assembly) ────────────────

_prefix!(x, p) = (x.name = "$p.$(x.name)"; x)

"""
    add!(net, asm::Assembly)

Merge an assembly into a network (or into another assembly, for nesting).
Names of its elements, shafts, electrical nodes, heat links and solver pairs
are prefixed with the assembly name. An element may belong to one assembly
only; environment nodes (`ThermalBoundary`) are shared and keep their names.
"""
function add!(net::FlowNetwork, a::Assembly)
    src = a.net
    for el in src.elements
        any(e -> e === el, net.elements) &&
            error("add!: element \"$(el.name)\" of assembly \"$(a.name)\" is already in the network")
    end
    off = length(net.elements)
    for el in src.elements
        push!(net.elements, _prefix!(el, a.name))
    end
    for e in src.edges
        push!(net.edges, PortEdge(e.src + off, e.src_port, e.dst + off, e.dst_port, e.back_edge))
    end
    for sh in src.shafts
        push!(net.shafts, _prefix!(sh, a.name))
    end
    for (i, p, st) in src.boundaries
        push!(net.boundaries, (i + off, p, st))
    end
    if !isnothing(src.seed_state)
        isnothing(net.seed_state) ||
            error("add!: assembly \"$(a.name)\" sets a seed state but the network already has one — use set_boundary! inside assemblies")
        net.seed_el, net.seed_state = src.seed_el + off, src.seed_state
    end
    # Electrical nodes: reuse a node already present (shared across assemblies).
    idx = Int[]
    for n in src.electrical
        k = findfirst(x -> x === n, net.electrical)
        if isnothing(k)
            n isa Alternator || n isa Pump || _prefix!(n, a.name)   # fluid elements already prefixed
            push!(net.electrical, n)
            k = length(net.electrical)
        end
        push!(idx, k)
    end
    for (s, d) in src.e_edges
        push!(net.e_edges, (idx[s], idx[d]))
    end
    for l in src.heat_links
        push!(net.heat_links, _prefix!(l, a.name))
    end
    append!(net.volumes, src.volumes)          # volumes refer to the (shared) element objects
    for fv in src.free_vars
        push!(net.free_vars, FreeVar("$(a.name).$(fv.name)", fv.get, fv.set!, fv.ref, fv.lower, fv.upper))
    end
    for t in src.targets
        push!(net.targets, Target("$(a.name).$(t.name)", t.f, t.value, t.scale))
    end
    push!(net.assemblies, a.name)
    append!(net.assemblies, ("$(a.name).$p" for p in src.assemblies))
    _invalidate_plan!(net)
    net
end

"""    assemblies(net) -> Vector{String} — full paths of all merged assemblies."""
assemblies(net::FlowNetwork) = copy(net.assemblies)

"""
    members(net, path) -> Vector

Elements, shafts and electrical nodes whose names lie under assembly `path`
(e.g. `"PCU1"` or `"PCU1.TAC"`).
"""
function members(net::FlowNetwork, path::AbstractString)
    pre = path * "."
    [x for x in Iterators.flatten((net.elements, net.shafts, net.electrical))
     if startswith(x.name, pre)] |> unique
end
