# ── Port-graph topology helpers ──────────────────────────────────────────────
# Shared by inventory (network layer) and accounting (output layer).

# Concrete outlet port names (HX `:outlet` is an alias of `:cold_outlet`).
_outlet_names(el::AbstractElement) = network_outlets(el)
_outlet_names(::HeatExchanger)     = [:cold_outlet, :hot_outlet]
_canon_outlet(el::AbstractElement, p::Symbol) = p
_canon_outlet(::HeatExchanger, p::Symbol)     = p == :outlet ? :cold_outlet : p
_canon_inlet(el::AbstractElement, p::Symbol)  = p
_canon_inlet(::HeatExchanger, p::Symbol)      = p == :inlet ? :cold_inlet : p

function _inlet_port(el::AbstractElement, port::Symbol)
    el isa HeatExchanger && return port === :hot_inlet ? el.hot_inlet : el.cold_inlet
    if el isa Mixer
        i = port === :inlet ? 1 : port === :bleed_inlet ? 2 :
            parse(Int, string(port)[7:end])
        return el.inlets[i]
    end
    el.inlet
end

# Ports joined inside an element (same stream): (inlet, outlet) pairs.
_internal_links(el::AbstractElement) = [(:inlet, :outlet)]
_internal_links(::HeatExchanger) = [(:hot_inlet, :hot_outlet), (:cold_inlet, :cold_outlet)]
_internal_links(el::Splitter) = [(:inlet, p) for p in network_outlets(el)]
_internal_links(el::Mixer) = [(p, :outlet) for p in network_required_inlets(el)]

# Connected components of the port graph: Dict (element index, port) → root.
function _stream_components(net::FlowNetwork)
    parent = Dict{Tuple{Int,Symbol},Tuple{Int,Symbol}}()
    find(x) = (haskey(parent, x) || (parent[x] = x);
               parent[x] == x ? x : (parent[x] = find(parent[x])))
    unite(a, b) = (ra = find(a); rb = find(b); ra == rb || (parent[ra] = rb))
    for (i, el) in enumerate(net.elements), (a, b) in _internal_links(el)
        unite((i, a), (i, b))
    end
    for e in net.edges
        unite((e.src, _canon_outlet(net.elements[e.src], e.src_port)),
              (e.dst, _canon_inlet(net.elements[e.dst], e.dst_port)))
    end
    Dict(k => find(k) for k in collect(keys(parent)))
end


# State that a back-edge closing onto a seeded inlet must honour: the network
# seed (if the edge closes onto the seed element's inlet) or a boundary state
# set on the destination port. `nothing` otherwise.
function _seeded_state(net::FlowNetwork, edge::PortEdge)
    el = net.elements[edge.dst]
    dp = _canon_inlet(el, edge.dst_port)
    for (i, p, st) in net.boundaries
        i == edge.dst && _canon_inlet(el, p) == dp && return st
    end
    if edge.dst == net.seed_el && dp == _canon_inlet(el, :inlet) && !isnothing(net.seed_state)
        return net.seed_state
    end
    nothing
end

# ── Mass-flow propagation ────────────────────────────────────────────────────
# Mass flow is pure bookkeeping (no element changes it except by splitting or
# mixing), so every port's flow follows from the seed and boundary flows.
# `one_pass!` uses this to give each back-edge the flow implied by the current
# seeds, which keeps the residual a pure function of the unknowns when a seed
# flow is itself an unknown (inventory mode, flow targets).

# (outlet, inlets, combine) rules per element: outlet flow from inlet flows.
_flow_rules(el::AbstractElement) = [(:outlet, [:inlet], ws -> ws[1])]
_flow_rules(::HeatExchanger) = [(:hot_outlet, [:hot_inlet], ws -> ws[1]),
                                (:cold_outlet, [:cold_inlet], ws -> ws[1])]
_flow_rules(el::Splitter) = [(p, [:inlet], ws -> ws[1] * el.fracs[k])
                             for (k, p) in enumerate(network_outlets(el))]
_flow_rules(el::Mixer) = [(:outlet, network_required_inlets(el), ws -> sum(ws))]

"""
Flows at every reachable port, Dict (element index, canonical port) → W,
propagated from the seed and boundary states. Back-edges that close onto a
seeded port do not overwrite it.
"""
function _propagate_flows(net::FlowNetwork)
    W = Dict{Tuple{Int,Symbol},Any}()
    seeded = Set{Tuple{Int,Symbol}}()
    if net.seed_el > 0 && !isnothing(net.seed_state)
        k = (net.seed_el, _canon_inlet(net.elements[net.seed_el], :inlet))
        W[k] = net.seed_state.W; push!(seeded, k)
    end
    for (i, p, st) in net.boundaries
        k = (i, _canon_inlet(net.elements[i], p))
        W[k] = st.W; push!(seeded, k)
    end
    edges = [((e.src, _canon_outlet(net.elements[e.src], e.src_port)),
              (e.dst, _canon_inlet(net.elements[e.dst], e.dst_port))) for e in net.edges]
    changed = true
    while changed
        changed = false
        for (i, el) in enumerate(net.elements), (out, ins, f) in _flow_rules(el)
            haskey(W, (i, out)) && continue
            all(p -> haskey(W, (i, p)), ins) || continue
            W[(i, out)] = f([W[(i, p)] for p in ins]); changed = true
        end
        for (s, d) in edges
            (haskey(W, s) && !haskey(W, d) && !(d in seeded)) || continue
            W[d] = W[s]; changed = true
        end
    end
    W
end
