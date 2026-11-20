# ── Energy and loss accounting ───────────────────────────────────────────────
#
# First law: where energy enters and leaves the system (heat sources, radiators,
# open streams, shaft work) and whether it balances.
#
# Second law: exergy (availability) relative to a dead state T0. The heat
# sources supply exergy to the working fluid; it leaves as net shaft power,
# is destroyed by irreversibilities (T0·S_gen: turbomachinery, pressure losses,
# finite-ΔT heat transfer, mixing), or is rejected (radiators, coolant and
# other open streams). Gouy–Stodola makes this an exact identity, so the
# ledger's closure residual is a correctness check on the whole model.
#
# Element interface (extend for new element types):
#
#     element_inlets(el)  -> Vector{Port}     ports carrying flow in
#     element_outlets(el) -> Vector{Port}     ports carrying flow out
#     element_role(el)    -> :adiabatic | :work_in | :work_out | :heat_in | :heat_out
#     element_losses(el, T0) -> Vector{Loss}  irreversibilities / rejection
#
# The default `element_losses` reports the whole entropy generation of an
# adiabatic element under one kind (`loss_kind(el)`); heat sources/sinks and
# heat exchangers split thermal and pressure-loss parts.

"""
    Loss(element, kind, exergy)

One entry of the loss ledger: `exergy` [W] destroyed (T0·S_gen) or rejected at
`element`. `kind` is one of `:aero` (turbomachinery irreversibility),
`:pressure` (total-pressure loss), `:heat_transfer` (finite-ΔT heat exchange),
`:mixing`, `:heat_rejection` (exergy leaving with rejected heat or an open
stream), `:electrical` (alternator, cable and converter losses),
`:parasitic` (power dumped in parasitic loads), `:heat_leak` (exergy lost
through heat links to environment nodes), `:pumping` (pump irreversibility),
`:imbalance` (numerical energy imbalance of an element).
"""
struct Loss
    element::String
    kind::Symbol
    exergy::Real
end

# ── Element interface ────────────────────────────────────────────────────────

_ports(ps...) = Port[p for p in ps if !isnothing(p)]

element_inlets(el::AbstractElement)  = _ports(el.inlet)
element_outlets(el::AbstractElement) = _ports(el.outlet)
element_inlets(el::HeatExchanger)    = _ports(el.hot_inlet, el.cold_inlet)
element_outlets(el::HeatExchanger)   = _ports(el.hot_outlet, el.cold_outlet)
element_outlets(el::Splitter)        = _ports(el.outlets...)
element_inlets(el::Mixer)            = _ports(el.inlets...)

element_role(::AbstractElement) = :adiabatic
element_role(::Compressor)      = :work_in
element_role(::Turbine)         = :work_out
element_role(::HeatSource)      = :heat_in
element_role(::Radiator)        = :heat_out
element_role(::Alternator)      = :alternator   # shaft in, electricity out, loss heats the gas
element_role(::Pump)            = :pump         # electricity in, all of it heats the coolant

# Heat an adiabatic element receives from heat links [W].
external_heat(::AbstractElement) = 0.0
external_heat(el::Duct)          = el.heated ? el.Q_ext : 0.0

loss_kind(::AbstractElement) = :pressure
loss_kind(::Compressor)      = :aero
loss_kind(::Turbine)         = :aero
loss_kind(::Mixer)           = :mixing
loss_kind(::Splitter)        = :mixing
loss_kind(::Pump)            = :pumping

_H(p::Port) = (s = p[]; s.W * enthalpy(s))
_S(p::Port) = (s = p[]; s.W * entropy(s))
_Hsum(ps) = sum(_H, ps; init = 0.0)
_Ssum(ps) = sum(_S, ps; init = 0.0)

# Entropy rise of a stream from its inlet pressure to its outlet pressure at
# the outlet temperature — the pressure-loss part of a heat-transfer element's
# entropy change (path: heat at constant inlet pressure, then lose pressure).
function _dS_pressure(pin::Port, pout::Port)
    si, so = pin[], pout[]
    so.W * (entropy(so.fluid, so.Tt, so.Pt) - entropy(so.fluid, so.Tt, si.Pt))
end

function element_losses(el::AbstractElement, T0)
    S_gen = _Ssum(element_outlets(el)) - _Ssum(element_inlets(el))
    [Loss(el.name, loss_kind(el), T0 * S_gen)]
end

# Heat source: only the pressure loss is an internal irreversibility; the
# heating itself is the exergy input (see `_exergy_in`).
function element_losses(el::HeatSource, T0)
    [Loss(el.name, :pressure, T0 * _dS_pressure(el.inlet, el.outlet))]
end

# Radiator: pressure loss, plus all exergy of the heat it rejects.
function element_losses(el::Radiator, T0)
    dSp   = _dS_pressure(el.inlet, el.outlet)
    dH    = _H(el.outlet) - _H(el.inlet)                       # < 0
    dS_th = _S(el.outlet) - _S(el.inlet) - dSp                 # < 0
    [Loss(el.name, :pressure, T0 * dSp),
     Loss(el.name, :heat_rejection, -(dH - T0 * dS_th))]
end

function element_losses(el::HeatExchanger, T0)
    dSp = _dS_pressure(el.hot_inlet, el.hot_outlet) +
          _dS_pressure(el.cold_inlet, el.cold_outlet)
    S_gen = _Ssum(element_outlets(el)) - _Ssum(element_inlets(el))
    [Loss(el.name, :heat_transfer, T0 * (S_gen - dSp)),
     Loss(el.name, :pressure, T0 * dSp)]
end

# Alternator: the loss heats the gas; the entropy that generates is the
# destroyed exergy (the rest re-enters the cycle as gas exergy).
function element_losses(el::Alternator, T0)
    dSp = _dS_pressure(el.inlet, el.outlet)
    S_gen = _S(el.outlet) - _S(el.inlet)
    [Loss(el.name, :electrical, T0 * (S_gen - dSp)),
     Loss(el.name, :pressure, T0 * dSp)]
end

# Duct: only its pressure loss here; heat through its port is booked per
# heat link (`_link_losses`).
element_losses(el::Duct, T0) = [Loss(el.name, :pressure, T0 * _dS_pressure(el.inlet, el.outlet))]

# Electrical side: per node (P_in, P_loss) bookkeeping.
_is_conversion(el) = !(el isa Union{Alternator,ElectricalLoad,Pump})
_elec_losses(net) = sum(el -> el.e.P_loss, (el for el in net.electrical if _is_conversion(el)); init = 0.0)

_wired(net, el) = any(x -> x === el, net.electrical)

# Exergy per joule of heat entering duct d's stream: exact from the duct's
# states (ΔEx_th / Q_ext), with the 1 − T0/T limit when no net heat enters.
function _exergy_per_heat(d::Duct, T0)
    Q = d.Q_ext
    if abs(_primal(Q)) > 0
        dS_th = _S(d.outlet) - _S(d.inlet) - _dS_pressure(d.inlet, d.outlet)
        1 - T0 * dS_th / Q
    else
        1 - T0 / d.inlet[].Tt
    end
end

# Exergy destroyed or lost by each heat link. Summed over the links of a duct
# this equals the exergy its stream loses through its heat port, exactly.
function _link_losses(net, T0)
    out = Loss[]
    for l in net.heat_links
        ea = l.a isa Duct ? _exergy_per_heat(l.a, T0) : 0.0
        eb = l.b isa Duct ? _exergy_per_heat(l.b, T0) : 0.0
        if l.kind === :loss
            # Credit the electrical component with the exergy its heat gives the stream.
            push!(out, Loss(l.a.name, :electrical, -l.Q * eb))
        else
            kind = (l.a isa Duct && l.b isa Duct) ? :heat_transfer : :heat_leak
            push!(out, Loss(l.name, kind, l.Q * (ea - eb)))
        end
    end
    out
end

# Heat into environment nodes through links, and electrical loss heat into streams.
_Q_env(net) = sum(l -> (l.b isa ThermalBoundary ? l.Q : 0.0) - (l.a isa ThermalBoundary ? l.Q : 0.0),
                  net.heat_links; init = 0.0)
_Q_loss_to_streams(net) = sum(l -> l.kind === :loss ? l.Q : 0.0, net.heat_links; init = 0.0)

# Exergy supplied to the fluid by a heat source: enthalpy rise minus T0 times
# the entropy rise of the constant-pressure heating part.
function _exergy_in(el::HeatSource, T0)
    dH = _H(el.outlet) - _H(el.inlet)
    dS_th = _S(el.outlet) - _S(el.inlet) - _dS_pressure(el.inlet, el.outlet)
    dH - T0 * dS_th
end

# ── System boundary: open streams ────────────────────────────────────────────

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

"""
Open streams of a solved network: `entrances` (boundary streams, and the seed
inlet when nothing feeds it — an open-loop model) and `exits` (outlets with no
downstream element). Each entry is `(label, Port, (element index, port))`.
"""
function _open_streams(net::FlowNetwork)
    fed = Set{Tuple{Int,Symbol}}()
    srcs = Set{Tuple{Int,Symbol}}()
    for e in net.edges
        push!(fed,  (e.dst, _canon_inlet(net.elements[e.dst], e.dst_port)))
        push!(srcs, (e.src, _canon_outlet(net.elements[e.src], e.src_port)))
    end
    entrances = Tuple{String,Port,Tuple{Int,Symbol}}[]
    for (i, port, _) in net.boundaries
        el = net.elements[i]
        cp_ = _canon_inlet(el, port)
        (i, cp_) in fed && continue                       # overridden by a loop closure
        p = _inlet_port(el, cp_)
        isnothing(p) || push!(entrances, ("$(el.name).$cp_", p, (i, cp_)))
    end
    seed_el = net.elements[net.seed_el]
    seed_port = _canon_inlet(seed_el, :inlet)
    if !((net.seed_el, seed_port) in fed)
        p = _inlet_port(seed_el, seed_port)
        isnothing(p) || push!(entrances, ("$(seed_el.name).$seed_port (seed)", p, (net.seed_el, seed_port)))
    end
    exits = Tuple{String,Port,Tuple{Int,Symbol}}[]
    for (i, el) in enumerate(net.elements)
        for port in _outlet_names(el)
            (i, port) in srcs && continue
            p = network_get_outlet(el, port)
            isnothing(p) || push!(exits, ("$(el.name).$port", p, (i, port)))
        end
    end
    (entrances = entrances, exits = exits)
end

# A heat source with negative duty (a cooler modeled as HeatSource) rejects
# heat. Classification only — the accounting math is the same either way.
_heat_role(el, dH) = (r = element_role(el);
                      r === :heat_in && _primal(dH) < 0 ? :heat_out : r)

# ── First law ────────────────────────────────────────────────────────────────

"""
    energy_balance(net_or_sol) -> NamedTuple

First-law accounting of a solved network [W]:

- `Q_in`        — heat added by `HeatSource`s (those with positive duty)
- `Q_out`       — heat removed by `Radiator`s and by `HeatSource`s with
                  negative duty (e.g. a bleed cooler)
- `Q_open`      — net enthalpy carried out by open streams (coolant loops,
                  or the implicit cooler of an open-loop model)
- `P_shaft`     — net shaft power (turbines − compressors)
- `P_alternator`— shaft power drawn by alternators
- `shaft_surplus` — `P_shaft − P_alternator` (≈ 0 on balanced shafts that
                  drive alternators; equals `P_shaft` without alternators)
- `P_electric`  — alternator electrical output
- `P_delivered`, `P_parasitic`, `P_elec_losses` — where it goes: useful
                  loads, parasitic loads, cables and converters
- `Q_env`       — net heat leaked to environment nodes through heat links
- `Q_loss_to_streams` — electrical loss heat routed into streams
- `P_pumps`, `P_pumps_external` — pump power, and the part not drawn from
                  the cycle's own electrical network
- `imbalance`   — Σ enthalpy imbalance of adiabatic elements (should be ≈ 0)
- `residual`    — `Q_in + P_pumps_external + Q_loss_to_streams − Q_out −
                  Q_open − Q_env − shaft_surplus − P_delivered − P_parasitic −
                  P_elec_losses` (≈ 0; equals `imbalance`)
- `by_element`  — `(name, role, value)` per element (heat/work [W], or the
                  imbalance for adiabatic elements)
"""
energy_balance(r::SolveResult) = energy_balance(r.net)

function energy_balance(net::FlowNetwork)
    Q_in = Q_out = W_in = W_out = imb = P_alt = P_el = P_pump = P_pump_ext = 0.0
    by_el = Tuple{String,Symbol,Real}[]
    for el in net.elements
        dH = _Hsum(element_outlets(el)) - _Hsum(element_inlets(el))
        role = _heat_role(el, dH)
        if role === :heat_in
            Q_in += dH
        elseif role === :heat_out
            Q_out -= dH
        elseif role === :work_in
            W_in += dH
        elseif role === :work_out
            W_out -= dH
        elseif role === :alternator
            P_alt += el.P_mech
            P_el  += el.e.P_out
            imb   += dH - (el.P_mech - el.e.P_out)
        elseif role === :pump
            P_pump += el.P_elec
            _wired(net, el) || (P_pump_ext += el.P_elec)
            imb += dH - el.P_elec
        else
            imb += dH - external_heat(el)
        end
        push!(by_el, (el.name, role, role in (:heat_out, :work_out) ? -dH :
                                     role === :alternator ? el.P_mech :
                                     role === :pump ? el.P_elec : dH))
    end
    os = _open_streams(net)
    Q_open = sum(p -> _H(p[2]), os.exits; init = 0.0) -
             sum(p -> _H(p[2]), os.entrances; init = 0.0)
    P = W_out - W_in
    surplus = P - P_alt
    P_del, P_par, P_eloss = electrical_power(net), parasitic_power(net), _elec_losses(net)
    Q_env, Q_l2s = _Q_env(net), _Q_loss_to_streams(net)
    (Q_in = Q_in, Q_out = Q_out, Q_open = Q_open, P_shaft = P,
     P_alternator = P_alt, shaft_surplus = surplus, P_electric = P_el,
     P_delivered = P_del, P_parasitic = P_par, P_elec_losses = P_eloss,
     Q_env = Q_env, Q_loss_to_streams = Q_l2s, P_pumps = P_pump,
     P_pumps_external = P_pump_ext, imbalance = imb,
     residual = Q_in + P_pump_ext + Q_l2s - Q_out - Q_open - Q_env - surplus -
                P_del - P_par - P_eloss,
     by_element = by_el)
end

# ── Second law ───────────────────────────────────────────────────────────────

"""
    LossLedger

Result of [`loss_ledger`](@ref). Fields: `T0` (dead state [K]), `exergy_in`
(`(name, value)` per heat source [W]), `Q_in` [W], `P_shaft` [W] (turbines −
compressors), `P_delivered` [W] (useful electrical loads), `P_out` [W] (net
useful output: delivered electricity plus any unbalanced shaft surplus; equals
`P_shaft` without an electrical network), `losses` (`Vector{Loss}`),
`closure` (exergy in − P_out − Σ losses; ≈ 0), `electrical` (whether the
model has an electrical network).
"""
struct LossLedger
    T0::Real
    Q_in::Real
    exergy_in::Vector{Tuple{String,Real}}
    P_shaft::Real
    P_delivered::Real
    P_out::Real
    losses::Vector{Loss}
    closure::Real
    electrical::Bool
end

"""
    default_dead_state(net) -> T0 [K]

Lowest available heat-sink temperature: the minimum of radiator sink
temperatures and the inlet temperatures of the system's real inflows
(boundary streams, and the seed inlet of an open-loop model); for a closed
loop with neither, the seed-element inlet temperature.
"""
function default_dead_state(net::FlowNetwork)
    Ts = Float64[]
    for el in net.elements
        el isa Radiator && push!(Ts, _primal(_sink_T(el)))
    end
    for (_, p) in _open_streams(net).entrances      # real boundary inflows (and an unfed seed)
        push!(Ts, _primal(p[].Tt))
    end
    if isempty(Ts)
        p = _inlet_port(net.elements[net.seed_el], :inlet)
        isnothing(p) && error("default_dead_state: solve the network first")
        push!(Ts, _primal(p[].Tt))
    end
    minimum(Ts)
end

"""
    loss_ledger(net_or_sol; T0 = default_dead_state(net)) -> LossLedger

Second-law (exergy) accounting of a solved network relative to dead state
`T0` [K]. The exergy the heat sources add to the working fluid equals the net
useful output (net shaft power, or delivered electrical power when the model
has an electrical network) plus every loss in the ledger — destroyed (T0·S_gen
of each irreversibility, alternator/cable/converter losses, parasitic loads)
or rejected (radiators, coolant and other open streams).
`closure` is the residual of that identity and should be ≈ 0.

Displaying the ledger prints the waterfall from exergy input to net power,
sorted by element; `losses_by_kind(ledger)` groups it by kind.
"""
loss_ledger(r::SolveResult; kw...) = loss_ledger(r.net; kw...)

function loss_ledger(net::FlowNetwork; T0::Real = default_dead_state(net))
    losses = Loss[]
    ex_in = Tuple{String,Real}[]
    Q_in = 0.0
    W_in = W_out = P_alt = 0.0
    for el in net.elements
        dH = _Hsum(element_outlets(el)) - _Hsum(element_inlets(el))
        role = _heat_role(el, dH)
        if role === :heat_in
            push!(ex_in, (el.name, _exergy_in(el, T0)))
            Q_in += dH
        elseif el isa HeatSource                       # negative duty: rejection
            push!(losses, Loss(el.name, :heat_rejection, -_exergy_in(el, T0)))
        elseif role === :work_in
            W_in += dH
        elseif role === :work_out
            W_out -= dH
        elseif role === :alternator
            P_alt += el.P_mech
            imb = dH - (el.P_mech - el.e.P_out)
            abs(_primal(imb)) > 0 && push!(losses, Loss(el.name, :imbalance, -imb))
        elseif role === :pump
            # Power from outside the model is an exergy input; drawn from the
            # cycle's own bus it is already inside the electrical bookkeeping.
            _wired(net, el) || push!(ex_in, ("$(el.name) (external power)", el.P_elec))
            imb = dH - el.P_elec
            abs(_primal(imb)) > 0 && push!(losses, Loss(el.name, :imbalance, -imb))
        elseif role === :adiabatic
            # Energy an adiabatic element creates/destroys is a model error,
            # but it must be booked for the exergy identity to close.
            imb = dH - external_heat(el)
            abs(_primal(imb)) > 0 && push!(losses, Loss(el.name, :imbalance, -imb))
        end
        append!(losses, element_losses(el, T0))
    end
    # Open streams: net exergy carried out, grouped per physical stream (mass
    # in = mass out per stream, so the arbitrary property reference states
    # cancel). Streams are the connected components of the port graph.
    os = _open_streams(net)
    comp = _stream_components(net)
    key(lp) = comp[lp[3]]
    for k in unique(key.([os.exits; os.entrances]))
        outs = [p for p in os.exits if key(p) == k]
        ins  = [p for p in os.entrances if key(p) == k]
        ex = sum(p -> _H(p[2]) - T0 * _S(p[2]), outs; init = 0.0) -
             sum(p -> _H(p[2]) - T0 * _S(p[2]), ins;  init = 0.0)
        label = isempty(outs) ? join(first.(ins), ", ") : join(first.(outs), ", ")
        push!(losses, Loss(label, :heat_rejection, ex))
    end
    # Electrical side: conversion losses and parasitic loads are destroyed
    # exergy (dissipated as heat to the surroundings).
    append!(losses, _link_losses(net, T0))
    for el in net.electrical
        el isa Union{Alternator,Pump} && continue
        if el isa ElectricalLoad
            el.parasitic && push!(losses, Loss(el.name, :parasitic, el.e.P_in))
        else
            push!(losses, Loss(el.name, :electrical, el.e.P_loss))
        end
    end
    P = W_out - W_in
    P_del = electrical_power(net)
    P_out = P - P_alt + P_del
    Ex = sum(last, ex_in; init = 0.0)
    closure = Ex - P_out - sum(l -> l.exergy, losses; init = 0.0)
    LossLedger(T0, Q_in, ex_in, P, P_del, P_out, losses, closure, !isempty(net.electrical))
end

"""
    losses_by_kind(ledger) -> Vector{Pair{Symbol,Real}}

Total exergy loss per kind, largest first.
"""
function losses_by_kind(L::LossLedger)
    kinds = unique(l.kind for l in L.losses)
    tot = [k => sum(l.exergy for l in L.losses if l.kind === k) for k in kinds]
    sort(tot; by = p -> -_primal(p.second))
end

"""
    ledger_rows(ledger) -> Vector{Tuple{String,Symbol,Float64}}

Losses aggregated per (element, kind) in ledger order, negligible entries
dropped — the rows of the printed waterfall and of `losswaterfall`.
"""
function ledger_rows(L::LossLedger)
    Ex = _primal(sum(last, L.exergy_in; init = 0.0))
    agg = Tuple{String,Symbol,Float64}[]
    for l in L.losses
        k = findfirst(a -> a[1] == l.element && a[2] === l.kind, agg)
        isnothing(k) ? push!(agg, (l.element, l.kind, _primal(l.exergy))) :
                       (agg[k] = (agg[k][1], agg[k][2], agg[k][3] + _primal(l.exergy)))
    end
    filter(a -> abs(a[3]) >= 1e-6 * max(1.0, abs(Ex)), agg)
end

function Base.show(io::IO, ::MIME"text/plain", L::LossLedger)
    kW(x) = lpad(string(round(_primal(x) / 1e3, digits = 3)), 11)
    Ex = sum(last, L.exergy_in; init = 0.0)
    pct(x) = Ex == 0 ? "" : lpad(string(round(100 * _primal(x) / _primal(Ex), digits = 2)), 8) * " %"
    println(io, "Loss ledger (exergy, dead state T0 = $(round(_primal(L.T0), digits = 2)) K)")
    println(io, "  heat input Q = ", strip(kW(L.Q_in)), " kW")
    for (name, ex) in L.exergy_in
        println(io, "  ", rpad("+ exergy in  $name", 40), kW(ex), " kW", pct(ex))
    end
    for (el, kind, ex) in ledger_rows(L)
        println(io, "  ", rpad("− $el  ($kind)", 40), kW(ex), " kW", pct(ex))
    end
    if L.electrical
        surplus = L.P_out - L.P_delivered
        abs(_primal(surplus)) >= 1e-6 * max(1.0, abs(_primal(Ex))) &&
            println(io, "  ", rpad("  (unbalanced shaft surplus)", 40), kW(surplus), " kW", pct(surplus))
        println(io, "  ", rpad("= net output (delivered electric)", 40), kW(L.P_out), " kW", pct(L.P_out))
    else
        println(io, "  ", rpad("= net shaft power", 40), kW(L.P_out), " kW", pct(L.P_out))
    end
    print(io, "  closure residual: ", round(_primal(L.closure), sigdigits = 3), " W")
end

Base.show(io::IO, L::LossLedger) =
    print(io, "LossLedger(T0=$(round(_primal(L.T0), digits=2)) K, $(length(L.losses)) losses, ",
          "P_out=$(round(_primal(L.P_out)/1e3, digits=3)) kW)")

"""
    losses_by_assembly(ledger, net; depth = 1) -> Vector{Pair{String,Real}}

Total exergy loss per assembly (paths up to `depth` levels, e.g. `"PCU1"` at
depth 1, `"PCU1.TAC"` at depth 2), largest first. Losses outside any
assembly are grouped under `"(top level)"`.
"""
function losses_by_assembly(L::LossLedger, net::FlowNetwork; depth::Int = 1)
    paths = filter(p -> count(==('.'), p) < depth, net.assemblies)
    sort!(paths; by = length, rev = true)                 # longest (deepest) match first
    tot = Dict{String,Real}()
    for l in L.losses
        k = findfirst(p -> startswith(l.element, p * "."), paths)
        key = isnothing(k) ? "(top level)" : paths[k]
        tot[key] = get(tot, key, 0.0) + l.exergy
    end
    sort([k => v for (k, v) in tot]; by = p -> -_primal(p.second))
end
