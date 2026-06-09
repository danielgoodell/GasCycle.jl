"""
FlowNetwork — assembles cycle elements into a solvable system.

The network is a directed port graph.  Each element exposes named inlet and
outlet ports; edges connect one element's outlet port to another's inlet port.
`one_pass!` traverses the graph in topological order, computing each element
once all its required inlets are available.

Circular dependencies (recuperators, closed loops) are handled by marking
certain edges as back-edges: they are seeded with the source element's stored
outlet from the *previous* pass, enabling fixed-point iteration to converge.

Usage:
  net = FlowNetwork()
  add!(net, comp, recup, reactor, turb)
  connect!(net, comp => recup => reactor => turb => comp)   # serial chain
  add_shaft!(net, shaft; drives=comp, driven_by=turb)
  add_hx_pair!(net, recup; hot=turb)                        # back-edge for hot side
  set_state!(net, comp; Pt=500e3, Tt=400.0, W=10.0, fluid=fluid)

  # Bleed/bypass branches use explicit port connections:
  connect_port!(net, splitter, :bleed_outlet, bleed_duct, :inlet)
  connect_port!(net, bleed_duct, :outlet, mixer, :bleed_inlet)

  sol = solve!(net)
"""

# ── Port edge ─────────────────────────────────────────────────────────────────

struct PortEdge
    src::Int          # source element index
    src_port::Symbol  # outlet port name on src
    dst::Int          # destination element index
    dst_port::Symbol  # inlet port name on dst
    back_edge::Bool   # true → seed from stored outlet; don't wait for src to run
end

# ── Network struct ────────────────────────────────────────────────────────────

mutable struct FlowNetwork
    elements::Vector{AbstractElement}
    edges::Vector{PortEdge}
    shafts::Vector{Shaft}
    seed_el::Int                       # element whose primary inlet is seeded
    seed_state::Union{FluidState,Nothing}
end

FlowNetwork() = FlowNetwork(AbstractElement[], PortEdge[], Shaft[], 0, nothing)

# ── Element lookup ────────────────────────────────────────────────────────────

function _el_idx(net::FlowNetwork, el::AbstractElement)
    idx = findfirst(e -> e === el, net.elements)
    isnothing(idx) && error("Element \"$(el.name)\" not found in network — did you call add!?")
    idx
end

# ── Port helpers ──────────────────────────────────────────────────────────────

"""Return the Port currently stored on an outlet port, or nothing."""
function _get_outlet(el::AbstractElement, port::Symbol)::Union{Port,Nothing}
    port == :outlet     && return hasproperty(el, :outlet)      ? el.outlet      :
                                  (el isa HeatExchanger         ? el.cold_outlet : nothing)
    port == :cold_outlet && el isa HeatExchanger && return el.cold_outlet
    port == :hot_outlet  && el isa HeatExchanger && return el.hot_outlet
    if el isa Splitter
        port == :bleed_outlet && length(el.outlets) >= 2 && return el.outlets[2]
        m = match(r"^outlet_(\d+)$", string(port))
        !isnothing(m) && return el.outlets[parse(Int, m.captures[1])]
    end
    nothing
end

"""Set an inlet port on an element."""
function _set_inlet!(el::AbstractElement, port::Symbol, p::Port)
    if port == :inlet
        if el isa HeatExchanger;  el.cold_inlet = p
        elseif el isa Mixer;      el.inlets[1]  = p
        else;                     el.inlet       = p
        end
    elseif port == :cold_inlet  && el isa HeatExchanger; el.cold_inlet = p
    elseif port == :hot_inlet   && el isa HeatExchanger; el.hot_inlet  = p
    elseif port == :bleed_inlet && el isa Mixer;         el.inlets[2]  = p
    else
        m = match(r"^inlet_(\d+)$", string(port))
        !isnothing(m) && el isa Mixer && (el.inlets[parse(Int, m.captures[1])] = p; return)
        error("Unknown inlet port :$port on element \"$(el.name)\"")
    end
end

"""Return the stored value of an inlet port, or nothing."""
function _get_stored_inlet(el::AbstractElement, port::Symbol)::Union{Port,Nothing}
    port == :inlet      && return hasproperty(el, :inlet) ? el.inlet :
                                  (el isa HeatExchanger   ? el.cold_inlet : nothing)
    port == :cold_inlet && el isa HeatExchanger && return el.cold_inlet
    port == :hot_inlet  && el isa HeatExchanger && return el.hot_inlet
    port == :bleed_inlet && el isa Mixer && length(el.inlets) >= 2 && return el.inlets[2]
    nothing
end

"""Inlet port names that must be satisfied before element can be computed."""
function _required_inlets(el::AbstractElement)::Vector{Symbol}
    el isa HeatExchanger && return [:cold_inlet, :hot_inlet]
    el isa Mixer         && return [i == 1 ? :inlet :
                                    i == 2 ? :bleed_inlet : Symbol("inlet_$i")
                                    for i in 1:el.n_inlets]
    [:inlet]
end

# ── Compute one element and store its outlets in avail ────────────────────────

const _AvailMap = Dict{Tuple{Int,Symbol}, Port}

function _compute_element!(net::FlowNetwork, el_idx::Int, avail::_AvailMap)
    el = net.elements[el_idx]

    if el isa HeatExchanger
        el.cold_inlet = avail[(el_idx, :cold_inlet)]
        el.hot_inlet  = avail[(el_idx, :hot_inlet)]
        compute_hx!(el)
        avail[(el_idx, :cold_outlet)] = el.cold_outlet
        avail[(el_idx, :hot_outlet)]  = el.hot_outlet
        avail[(el_idx, :outlet)]      = el.cold_outlet   # alias for serial chain

    elseif el isa Mixer
        for i in 1:el.n_inlets
            p = i == 1 ? :inlet : (i == 2 ? :bleed_inlet : Symbol("inlet_$i"))
            el.inlets[i] = avail[(el_idx, p)]
        end
        compute!(el)
        avail[(el_idx, :outlet)] = el.outlet

    elseif el isa Splitter
        compute!(el, avail[(el_idx, :inlet)])
        avail[(el_idx, :outlet)] = el.outlets[1]
        if length(el.outlets) >= 2
            avail[(el_idx, :bleed_outlet)] = el.outlets[2]
        end
        for i in 3:length(el.outlets)
            avail[(el_idx, Symbol("outlet_$i"))] = el.outlets[i]
        end

    else
        inlet = avail[(el_idx, :inlet)]
        out   = compute!(el, inlet)
        avail[(el_idx, :outlet)] = out
    end
end

# ── Public API ────────────────────────────────────────────────────────────────

"""Add elements to the network."""
function add!(net::FlowNetwork, elements::AbstractElement...)
    for el in elements
        isa(el, Shaft) ? push!(net.shafts, el) : push!(net.elements, el)
    end
end

"""
    connect!(net, el1 => el2 => el3 => ...)

Wire a serial flow chain.  Each element's `:outlet` feeds the next element's
`:inlet`.  For closed loops, write `... => last_el => first_el`; the closing
back-edge is stripped (the loop is closed by `set_state!`).
"""
function connect!(net::FlowNetwork, chain)
    ordered = _flatten_chain(chain)
    if length(ordered) > 1 && ordered[end] === ordered[1]
        ordered = ordered[1:end-1]
    end
    idxs = [_el_idx(net, el) for el in ordered]
    for i in 1:length(idxs)-1
        push!(net.edges, PortEdge(idxs[i], :outlet, idxs[i+1], :inlet, false))
    end
end

_flatten_chain(p::Pair) = [_flatten_chain(p.first)..., _flatten_chain(p.second)...]
_flatten_chain(el)      = [el]

"""
    connect_port!(net, src, src_port, dst, dst_port)

Add an explicit port-to-port edge.  Used for bleed and bypass branches where
the default `:outlet` → `:inlet` mapping doesn't apply.
"""
function connect_port!(net::FlowNetwork,
                       src::AbstractElement, src_port::Symbol,
                       dst::AbstractElement, dst_port::Symbol)
    push!(net.edges,
          PortEdge(_el_idx(net, src), src_port, _el_idx(net, dst), dst_port, false))
end

"""Add a shaft and link it to its producers and consumers."""
function add_shaft!(net::FlowNetwork, shaft::Shaft;
                    drives,
                    driven_by)
    push!(net.shafts, shaft)
    producers = isa(driven_by, AbstractVector) ? driven_by : [driven_by]
    consumers = isa(drives,    AbstractVector) ? drives    : [drives]
    link!(shaft, producers, consumers)
end

"""
    add_hx_pair!(net, hx; hot[, cold])

Register the back-edge from `hot`'s outlet to `hx`'s hot inlet.  The cold-side
path must already be wired by `connect!`.  The optional `cold` keyword is
accepted for backward compatibility but is no longer used.
"""
function add_hx_pair!(net::FlowNetwork, hx::HeatExchanger;
                      hot::AbstractElement, cold::Union{AbstractElement,Nothing}=nothing)
    push!(net.edges,
          PortEdge(_el_idx(net, hot), :outlet, _el_idx(net, hx), :hot_inlet, true))
end

"""Set the thermodynamic inlet state at the first element in the loop."""
function set_state!(net::FlowNetwork, first_element::AbstractElement;
                    Pt::Float64, Tt::Float64, W::Float64, fluid::FluidProperties)
    net.seed_el    = _el_idx(net, first_element)
    net.seed_state = FluidState(Pt, Tt, W, fluid)
end

# ── one_pass! ─────────────────────────────────────────────────────────────────

"""
    one_pass!(net)

Propagate thermodynamic state through all elements once, in topological order.

Back-edge inlets (e.g., recuperator hot side) are pre-seeded with the source
element's stored outlet from the *previous* pass, enabling fixed-point iteration
to resolve circular dependencies.
"""
function one_pass!(net::FlowNetwork)
    isnothing(net.seed_state) && error("call set_state! before solving")

    avail = _AvailMap()

    # Seed the initial element's inlet
    avail[(net.seed_el, :inlet)] = Port(net.seed_state)

    # Pre-seed back-edge inlets from the source element's stored outlet
    for edge in net.edges
        edge.back_edge || continue
        src_out = _get_outlet(net.elements[edge.src], edge.src_port)
        if isnothing(src_out)
            # First pass: no stored value yet — fall back to seed state
            avail[(edge.dst, edge.dst_port)] = Port(net.seed_state)
        else
            avail[(edge.dst, edge.dst_port)] = src_out
        end
    end

    # Topological traversal: repeatedly find elements whose inlets are all
    # available, compute them, then propagate their outlets downstream.
    processed = falses(length(net.elements))
    changed   = true
    while changed
        changed = false
        for el_idx in 1:length(net.elements)
            processed[el_idx] && continue
            net.elements[el_idx] isa Shaft && (processed[el_idx] = true; continue)

            # Check all required inlets are available
            all(p -> haskey(avail, (el_idx, p)),
                _required_inlets(net.elements[el_idx])) || continue

            _compute_element!(net, el_idx, avail)
            processed[el_idx] = true
            changed = true

            # Propagate this element's outlets to downstream inlets via forward edges
            for edge in net.edges
                (edge.src == el_idx && !edge.back_edge) || continue
                out = get(avail, (el_idx, edge.src_port), nothing)
                isnothing(out) && continue
                avail[(edge.dst, edge.dst_port)] = out
                # Normalize :inlet → :cold_inlet for HeatExchanger destinations
                if edge.dst_port == :inlet && net.elements[edge.dst] isa HeatExchanger
                    avail[(edge.dst, :cold_inlet)] = out
                end
            end
        end
    end

    # Shaft speed broadcast
    for sh in net.shafts
        _broadcast_speed!(sh)
    end

    nothing
end
