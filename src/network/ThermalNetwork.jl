# ── Heat links: wiring, per-pass application, environment coupling ───────────

function _check_stream(net::FlowNetwork, d, who)
    d isa Duct || error("$who: stream endpoints must be Ducts (heat ports), got $(typeof(d))")
    any(e -> e === d, net.elements) || error("$who: duct \"$(d.name)\" is not in the network")
    d.heated = true
end

function _check_endpoint(net, d, who)
    d isa ThermalBoundary && return
    _check_stream(net, d, who)
end

function _add_link!(net, kind, a, b, G, name)
    push!(net.heat_links, HeatLink(String(name), kind, a, b, G, 0.0, 1.0))
    net.heat_links[end]
end

"""
    add_conductance!(net, a, b; UA, name)

Linear heat link Q = UA·(T_a − T_b) [W], positive from `a` to `b`. Endpoints
are stream heat ports (`Duct`s in the network) or `ThermalBoundary` nodes:
wall conduction between two streams (crossover duct, shared wall), or a
convective/conductive leak from a stream to the surroundings.
"""
function add_conductance!(net::FlowNetwork, a, b; UA::Real,
                          name::AbstractString = "$(a.name)→$(b.name)")
    _check_endpoint(net, a, "add_conductance!"); _check_endpoint(net, b, "add_conductance!")
    _add_link!(net, :conductance, a, b, UA, name)
end

"""
    add_radiation!(net, a, b; area, emissivity = 1.0, view_factor = 1.0, name)

Radiative heat link Q = σ·ε·A·F·(T_a⁴ − T_b⁴) [W], typically from a stream's
duct or casing to an environment node. `emissivity·view_factor` is the
effective gray-body exchange factor — give the effective value for two-surface
exchange.
"""
function add_radiation!(net::FlowNetwork, a, b; area::Real, emissivity::Real = 1.0,
                        view_factor::Real = 1.0, name::AbstractString = "$(a.name)→$(b.name)")
    _check_endpoint(net, a, "add_radiation!"); _check_endpoint(net, b, "add_radiation!")
    _add_link!(net, :radiation, a, b, emissivity * area * view_factor, name)
end

"""
    add_loss_heat!(net, component, duct; fraction = 1.0, name)

Route `fraction` of an electrical component's dissipation (cable, converter,
…) into a stream — e.g. a PCU cold plate on the coolant loop. Losses not
routed anywhere go to the surroundings.
"""
function add_loss_heat!(net::FlowNetwork, comp, duct; fraction::Real = 1.0,
                        name::AbstractString = "$(comp.name)→$(duct.name)")
    any(e -> e === comp, net.electrical) && !(comp isa Union{ElectricalLoad,Alternator,Pump}) ||
        error("add_loss_heat!: \"$(comp.name)\" must be a cable, converter or bus in the electrical network")
    _check_stream(net, duct, "add_loss_heat!")
    _add_link!(net, :loss, comp, duct, fraction, name)
end

# Net heat into each heated duct from the links' current Q values.
function _apply_heat_links!(net::FlowNetwork)
    isempty(net.heat_links) && return nothing
    for el in net.elements
        el isa Duct && el.heated && (el.Q_ext = 0.0)
    end
    for l in net.heat_links
        l.a isa Duct && (l.a.Q_ext -= l.Q)
        l.b isa Duct && (l.b.Q_ext += l.Q)
    end
    nothing
end

# Before a solve: start each coupling unknown at its model value (after a pass
# has computed the states) and set its residual scale.
function _init_couplings!(net::FlowNetwork)
    for l in net.heat_links
        l.Q = heat_model(l)
        l.scale = max(abs(Float64(_primal(l.Q))), 1.0)
    end
    for el in net.elements
        if el isa Pump && !isnothing(el.inlet)
            el.P_elec = pump_power(el)
            el.scale = max(abs(Float64(_primal(el.P_elec))), 1.0)
        end
    end
    nothing
end

# ── Environment nodes ────────────────────────────────────────────────────────

"""    environment_nodes(net) -> Vector{ThermalBoundary} referenced by links or radiators."""
function environment_nodes(net::FlowNetwork)
    nodes = ThermalBoundary[]
    add(n) = n isa ThermalBoundary && !any(x -> x === n, nodes) && push!(nodes, n)
    for l in net.heat_links
        add(l.a); add(l.b)
    end
    for el in net.elements
        el isa Radiator && !isnothing(el.sink) && add(el.sink)
    end
    nodes
end

"""
    set_environment!(net, temps)

Set environment-node temperatures by name, e.g. from an external spacecraft
thermal model: `set_environment!(net, Dict("PanelA" => 245.0, "Deep space" => 4.0))`.
"""
function set_environment!(net::FlowNetwork, temps)
    nodes = environment_nodes(net)
    for (name, T) in temps
        k = findfirst(n -> n.name == name, nodes)
        isnothing(k) && error("set_environment!: no environment node \"$name\"")
        nodes[k].T = T
    end
    net
end

"""
    environment_heat(net_or_sol) -> Vector{Pair{String,Real}}

Heat into each environment node [W] (links and radiators combined) — the heat
loads to hand back to an external spacecraft thermal model.
"""
function environment_heat(net::FlowNetwork)
    nodes = environment_nodes(net)
    Q = Dict{String,Real}(n.name => 0.0 for n in nodes)
    for l in net.heat_links
        l.b isa ThermalBoundary && (Q[l.b.name] += l.Q)
        l.a isa ThermalBoundary && (Q[l.a.name] -= l.Q)
    end
    for el in net.elements
        el isa Radiator && !isnothing(el.sink) && !isnothing(el.inlet) &&
            (Q[el.sink.name] += Q_rejected(el))
    end
    [n.name => Q[n.name] for n in nodes]
end

"""    heat_link_report([io,] net_or_sol)   — table of heat links and their flows."""
heat_link_report(net::FlowNetwork) = heat_link_report(stdout, net)
function heat_link_report(io::IO, net::FlowNetwork)
    isempty(net.heat_links) && return println(io, "  (no heat links)")
    w = maximum(length(l.name) for l in net.heat_links)
    println(io, "  ", rpad("link", w), "  ", rpad("kind", 12), lpad("T_a [K]", 10),
            lpad("T_b [K]", 10), lpad("Q [W]", 12))
    Ts(d) = d isa Union{Duct,ThermalBoundary} ? string(round(_primal(_endpoint_T(d)), digits = 2)) : "—"
    for l in net.heat_links
        println(io, "  ", rpad(l.name, w), "  ", rpad(string(l.kind), 12), lpad(Ts(l.a), 10),
                lpad(Ts(l.b), 10), lpad(string(round(_primal(l.Q), sigdigits = 6)), 12))
    end
    nothing
end
