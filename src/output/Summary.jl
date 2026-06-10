"""
NPSS-style solution summary.

`summary(sol)` prints a station table (Tt/Pt/W/h/s at every port, in physical
flow order), a component table (work, heat, efficiency, pressure ratio per
element), and cycle totals.  `stations(sol)` returns the station list as data
for plotting and post-processing.
"""

using Printf
import ForwardDiff

# ── Dual-safe scalar extraction ───────────────────────────────────────────────

_scalar(x::AbstractFloat) = Float64(x)
_scalar(x::Real) = x isa ForwardDiff.Dual ? _scalar(ForwardDiff.value(x)) : Float64(x)

# ── Flow-order station traversal ──────────────────────────────────────────────

"""Outlet port reached when the flow enters `el` through `entry`."""
_exit_port(el::AbstractElement, entry::Symbol) = :outlet
_exit_port(el::HeatExchanger, entry::Symbol) =
    entry == :hot_inlet ? :hot_outlet : :cold_outlet

"""Does an edge labeled `edge_port` carry the flow that left through `exitp`?"""
function _port_matches(el::AbstractElement, edge_port::Symbol, exitp::Symbol)
    edge_port == exitp && return true
    el isa HeatExchanger && exitp == :cold_outlet && edge_port == :outlet
end

_port_suffix(port::Symbol) = port == :outlet ? "out" :
                             port == :inlet  ? "in"  :
                             replace(string(port), "outlet" => "out", "inlet" => "in")

"""Canonical name for an outlet port (resolves the HX `:outlet` alias)."""
_normal_out(el::AbstractElement, port::Symbol) = port
_normal_out(el::HeatExchanger, port::Symbol) =
    port == :outlet ? :cold_outlet : port

"""
    stations(net::FlowNetwork; branches=true) -> Vector{Pair{String,FluidState}}

Ordered list of `"Element.port" => FluidState` stations.  The main flow path
is walked first (following back-edges through heat-exchanger hot sides, so a
closed loop reads in physical order); branch sub-paths (bleed flows) follow
unless `branches=false`.
"""
function stations(net::FlowNetwork; branches::Bool = true)
    sts, _ = _station_paths(net; branches)
    sts
end

stations(r::SolveResult; kwargs...) = stations(r.net; kwargs...)

"""
Walk the port graph.  Returns `(stations, closed)` where `closed` is true if
the main path returned to the seed element (a closed cycle).
"""
function _station_paths(net::FlowNetwork; branches::Bool = true)
    isnothing(net.seed_state) && error("stations: call set_state! and solve! first")

    sts     = Pair{String,Any}[]
    emitted = Set{Tuple{Int,Symbol}}()
    visited = Set{Int}()
    closed  = false

    seed = net.elements[net.seed_el]
    push!(sts, "$(seed.name).in" => net.seed_state)

    # Walk one flow path, emitting each element's exit station, until the loop
    # closes, the path ends, or an already-visited element is reached.
    function walk!(el_idx::Int, entry::Symbol, main::Bool)
        for _ in 1:length(net.edges) + 1   # guard against malformed graphs
            el    = net.elements[el_idx]
            exitp = _exit_port(el, entry)
            push!(visited, el_idx)
            (el_idx, exitp) in emitted && return
            out = _get_outlet(el, exitp)
            isnothing(out) && return
            push!(sts, "$(el.name).$(_port_suffix(exitp))" => out[])
            push!(emitted, (el_idx, exitp))

            i = findfirst(e -> e.src == el_idx &&
                               _port_matches(el, e.src_port, exitp), net.edges)
            isnothing(i) && return
            edge = net.edges[i]
            if edge.dst == net.seed_el && edge.dst_port == :inlet
                main && (closed = true)
                return
            end
            el_idx, entry = edge.dst, edge.dst_port
        end
    end

    walk!(net.seed_el, :inlet, true)

    if branches
        # Branch sub-paths: any outlet port on a visited element that the
        # main walk didn't pass through (e.g. a splitter's bleed outlet).
        for edge in net.edges
            edge.src in visited || continue
            el = net.elements[edge.src]
            np = _normal_out(el, edge.src_port)
            (edge.src, np) in emitted && continue
            out = _get_outlet(el, np)
            isnothing(out) && continue
            push!(sts, "$(el.name).$(_port_suffix(np))" => out[])
            push!(emitted, (edge.src, np))
            edge.dst in visited || walk!(edge.dst, edge.dst_port, false)
        end

        # Boundary streams (set_boundary!) and their dangling outlet ports —
        # e.g. the coolant side of a heat-rejection exchanger.
        for (el_idx, port, state) in net.boundaries
            el = net.elements[el_idx]
            push!(sts, "$(el.name).$(_port_suffix(port))" => state)
        end
        for el_idx in sort(collect(visited))
            el = net.elements[el_idx]
            for port in network_outlets(el)
                np = _normal_out(el, port)
                np == port || continue        # skip aliases (HX :outlet)
                (el_idx, np) in emitted && continue
                out = _get_outlet(el, np)
                isnothing(out) && continue
                push!(sts, "$(el.name).$(_port_suffix(np))" => out[])
                push!(emitted, (el_idx, np))
            end
        end
    end

    (sts, closed)
end

# ── Station table ─────────────────────────────────────────────────────────────

function _print_stations(io::IO, net::FlowNetwork)
    sts  = stations(net)
    wide = maximum(length(first(p)) for p in sts; init=7)

    @printf(io, "  %-*s %10s %10s %8s %10s %10s\n", wide, "Station",
            "Tt [K]", "Pt [kPa]", "W [kg/s]", "h [kJ/kg]", "s [kJ/kgK]")
    for (label, s) in sts
        @printf(io, "  %-*s %10.2f %10.2f %8.4f %10.2f %10.4f\n", wide, label,
                _scalar(s.Tt), _scalar(s.Pt) / 1e3, _scalar(s.W),
                _scalar(enthalpy(s)) / 1e3, _scalar(entropy(s)) / 1e3)
    end
end

# ── Component table ───────────────────────────────────────────────────────────

_describe(el::AbstractElement) = ""
_describe(el::Duct)       = @sprintf("dP/P=%.3f", _scalar(el.dPqP))
_describe(el::Splitter)   = "fracs=" * string(el.fracs)
_describe(el::Mixer)      = ""

function _describe(el::Union{Compressor,Turbine})
    isnothing(el.inlet) && return "(not computed)"
    P_kW  = _scalar(specific_work(el) * el.inlet[].W) / 1e3
    η_lbl = el.η_type == :isentropic ? "ηad" : "ηp"
    str   = @sprintf("PR=%.4f  %s=%.4f  P=%.2f kW",
                     _scalar(pressure_ratio(el)), η_lbl, _scalar(el.η_poly), P_kW)
    if el.mode == :off_design
        s = el.inlet[]
        str *= @sprintf("  Wc=%.4f  Nc=%.1f", el.Wc_map,
                        corrected_speed(el.N_shaft, _scalar(s.Tt)))
    end
    str
end

function _describe(el::HeatExchanger)
    str = @sprintf("ε=%.4f", _scalar(el.ε))
    isnothing(el.hot_inlet) ||
        (str *= @sprintf("  Q=%.2f kW", _scalar(Q_transferred(el)) / 1e3))
    str
end

function _describe(el::HeatSource)
    if !isnothing(el.inlet) && !isnothing(el.outlet)
        Q = (enthalpy(el.outlet[]) - enthalpy(el.inlet[])) * el.inlet[].W
        @sprintf("Q=%.2f kW  dP/P=%.3f", _scalar(Q) / 1e3, _scalar(el.dPqP))
    else
        @sprintf("Q=%.2f kW  dP/P=%.3f", _scalar(el.Q) / 1e3, _scalar(el.dPqP))
    end
end

_describe(sh::Shaft) = begin
    str = @sprintf("N=%.0f rpm  balance=%.2f kW", sh.N, _scalar(power_balance(sh)) / 1e3)
    sh.P_load == 0.0 || (str *= @sprintf("  P_load=%.2f kW", sh.P_load / 1e3))
    str
end

function _print_components(io::IO, net::FlowNetwork)
    rows = [(el.name, string(nameof(typeof(el))), _describe(el))
            for el in vcat(net.elements, net.shafts)]
    wn = maximum(length(r[1]) for r in rows; init=4)
    wt = maximum(length(r[2]) for r in rows; init=4)
    for (name, typ, desc) in rows
        @printf(io, "  %-*s %-*s %s\n", wn, name, wt, typ, desc)
    end
end

# ── Cycle totals + entry point ────────────────────────────────────────────────

"""
    summary([io,] sol::SolveResult)

Print an NPSS-style report of the solved cycle: solver status, station table
in flow order, component table, and cycle totals.
"""
function Base.summary(io::IO, r::SolveResult)
    println(io, "GasCycle solution — ", r.status, " (", r.iterations,
            " iters, |F|=", @sprintf("%.3g", _scalar(r.residual_norm)), ")")
    println(io)
    _print_stations(io, r.net)
    println(io)
    _print_components(io, r.net)
    println(io)

    P_net = _scalar(net_power(r))
    η     = _scalar(cycle_efficiency(r))
    @printf(io, "  Net shaft power:  %10.2f kW\n", P_net / 1e3)
    η == 0.0 || @printf(io, "  Heat input:       %10.2f kW\n", P_net / η / 1e3)
    η == 0.0 || @printf(io, "  Cycle efficiency: %10.2f %%\n", η * 100)
    nothing
end

Base.summary(r::SolveResult) = summary(stdout, r)
