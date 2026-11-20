# ── Electrical network: wiring and power-flow solve ──────────────────────────

_e_name(el) = el.name

function _e_idx!(net::FlowNetwork, el)
    k = findfirst(e -> e === el, net.electrical)
    isnothing(k) || return k
    hasmethod(e_kind_out, Tuple{typeof(el)}) ||
        error("\"$(_e_name(el))\" is not an electrical node")
    push!(net.electrical, el)
    length(net.electrical)
end

"""
    connect_electrical!(net, a => b => c ...)

Wire electrical nodes in the order of power flow, from an `Alternator`
through cables, converters and buses to loads. Nodes are registered on first
use. A bus may feed several branches — connect each from the bus:

    connect_electrical!(net, alt => feeder => rect => dclink)
    connect_electrical!(net, dclink => plr)
    connect_electrical!(net, dclink => dcdc => bus => user)

The network must be a tree (each node fed once, sources are alternators) and
each connection must join matching kinds (:ac/:dc).
"""
function connect_electrical!(net::FlowNetwork, chain)
    items = _flatten_chain(chain)
    for (a, b) in zip(items[1:end-1], items[2:end])
        # Assemblies connect through exposed nodes: :elec as a source,
        # :elec_in (or :elec) as a destination.
        a = _resolve_e(a, :elec)
        b = b isa Assembly && haskey(b.eports, :elec_in) ? _resolve_e(b, :elec_in) : _resolve_e(b, :elec)
        ka, kb = e_kind_out(a), e_kind_in(b)
        isnothing(ka) && error("connect_electrical!: \"$(_e_name(a))\" is a load and cannot feed \"$(_e_name(b))\"")
        isnothing(kb) && error("connect_electrical!: \"$(_e_name(b))\" is a source and cannot be fed")
        ka === kb || error("connect_electrical!: \"$(_e_name(a))\" outputs :$ka but \"$(_e_name(b))\" takes :$kb")
        ia, ib = _e_idx!(net, a), _e_idx!(net, b)
        any(e -> e[2] == ib, net.e_edges) &&
            error("connect_electrical!: \"$(_e_name(b))\" is already fed (the electrical network must be a tree)")
        push!(net.e_edges, (ia, ib))
    end
    net
end

# Sweep order: sources first, then breadth-first down the tree.
function _e_order(net::FlowNetwork)
    n = length(net.electrical)
    children = [Int[] for _ in 1:n]
    fed = falses(n)
    for (a, b) in net.e_edges
        push!(children[a], b)
        fed[b] = true
    end
    order = Int[]
    for i in 1:n
        fed[i] && continue
        isnothing(e_kind_in(net.electrical[i])) ||
            error("electrical node \"$(_e_name(net.electrical[i]))\" has no feed")
        push!(order, i)
    end
    head = 1
    while head <= length(order)
        append!(order, children[order[head]])
        head += 1
    end
    length(order) == n || error("electrical network has a cycle")
    parent = zeros(Int, n)
    for (a, b) in net.e_edges
        parent[b] = a
    end
    (order, children, parent)
end

"""
Solve the electrical power flow: voltages down the tree, power demand up it,
repeated until resistive drops settle. Called from `one_pass!`.
"""
function _solve_electrical!(net::FlowNetwork; maxiter::Int = 50, rtol = 1e-14)
    isempty(net.electrical) && return nothing
    order, children, parent = _e_order(net)
    els = net.electrical
    V_prev = fill(NaN, length(els))
    for _ in 1:maxiter
        for i in order
            el = els[i]
            if parent[i] == 0
                _e_forward!(el, nothing, nothing)
            else
                p = els[parent[i]].e
                _e_forward!(el, p.V_out, p.f)
            end
        end
        for i in Iterators.reverse(order)
            P_out = sum(j -> els[j].e.P_in, children[i]; init = 0.0)
            _e_backward!(els[i], P_out)
        end
        V = [_primal(els[i].e.V_out) for i in eachindex(els)]
        done = all(i -> abs(V[i] - V_prev[i]) <= rtol * max(abs(V[i]), 1.0), eachindex(V))
        V_prev = V
        done && break
    end
    nothing
end

# ── Convenience queries ──────────────────────────────────────────────────────

_loads(net, parasitic::Bool) =
    (el for el in net.electrical if el isa ElectricalLoad && el.parasitic == parasitic)

"""    electrical_power(net_or_sol) -> Σ power delivered to non-parasitic loads [W]."""
electrical_power(net::FlowNetwork) = sum(el -> el.e.P_in, _loads(net, false); init = 0.0)

"""    parasitic_power(net_or_sol) -> Σ power dumped in parasitic loads [W]."""
parasitic_power(net::FlowNetwork) = sum(el -> el.e.P_in, _loads(net, true); init = 0.0)

"""
    add_speed_control!(net, plr, shaft; N)

Parasitic-load speed control: the solver sets the parasitic load's power
(`plr.P ≥ 0`) so the shaft runs at speed `N` [rpm]. Use with an off-design
shaft (speed unknown, power-balance residual).
"""
function add_speed_control!(net::FlowNetwork, plr::ElectricalLoad, shaft::Shaft; N::Real,
                            ref::Real = NaN)
    plr.parasitic || @warn "add_speed_control!: \"$(plr.name)\" is not a parasitic load"
    add_free!(net, plr, :P; lower = 0.0, ref = ref, name = "$(plr.name).P")
    add_target!(net, n -> shaft.N, N; name = "$(shaft.name).N")
end

"""
    electrical_report([io,] net_or_sol)

Table of the electrical network: per node, voltage, frequency, current, power
in/out and loss.
"""
electrical_report(net::FlowNetwork) = electrical_report(stdout, net)

function electrical_report(io::IO, net::FlowNetwork)
    isempty(net.electrical) && return println(io, "  (no electrical network)")
    order, _, parent = _e_order(net)
    f(x) = round(_primal(x), sigdigits = 5)
    w = maximum(length(_e_name(el)) for el in net.electrical)
    println(io, "  ", rpad("node", w), "  ", rpad("type", 15),
            lpad("V_out [V]", 11), lpad("f [Hz]", 9), lpad("I [A]", 10),
            lpad("P_in [kW]", 11), lpad("P_out [kW]", 11), lpad("loss [kW]", 11))
    for i in order
        el = net.electrical[i]; e = el.e
        typ = el isa ElectricalLoad ? (el.parasitic ? "ParasiticLoad" : "ElectricalLoad") :
              string(nameof(typeof(el)))
        println(io, "  ", rpad(_e_name(el), w), "  ", rpad(typ, 15),
                lpad(f(e.V_out), 11), lpad(f(e.f), 9), lpad(f(e.I), 10),
                lpad(f(e.P_in / 1e3), 11), lpad(f(e.P_out / 1e3), 11),
                lpad(f(e.P_loss / 1e3), 11))
    end
    nothing
end
