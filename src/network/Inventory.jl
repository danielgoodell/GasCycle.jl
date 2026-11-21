# ── Volumes and loop inventory ───────────────────────────────────────────────
#
# Each component volume holds m = ∫ρ dV. Along a volume the temperature and
# pressure are taken to vary linearly from inlet to outlet state, and the
# density is integrated by Gauss–Legendre quadrature — exact enough for the
# strong axial gradient of a recuperator (for an ideal gas with a linear T
# profile the exact answer is P·V/(R·T_lm); 5 points match it to ~1e-7 at a
# 2:1 temperature ratio). Works for any fluid backend and is differentiable.

"""A gas volume attached to a component (one side of a heat exchanger)."""
mutable struct VolumeSpec
    el::Any          # element
    side::Symbol     # :main, :hot, :cold
    V::Real          # [m³]
    n::Int           # quadrature points
end

# Gauss–Legendre nodes/weights on [0, 1].
const _GL = Dict(
    1 => ([0.5], [1.0]),
    2 => ([0.5 - sqrt(3)/6, 0.5 + sqrt(3)/6], [0.5, 0.5]),
    3 => ([0.5 - sqrt(15)/10, 0.5, 0.5 + sqrt(15)/10], [5/18, 8/18, 5/18]),
    4 => (0.5 .+ 0.5 .* [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526],
          0.5 .* [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538]),
    5 => (0.5 .+ 0.5 .* [-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640],
          0.5 .* [0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891]),
    6 => (0.5 .+ 0.5 .* [-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
                          0.2386191860831969, 0.6612093864662645, 0.9324695142031521],
          0.5 .* [0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
                  0.4679139345726910, 0.3607615730481386, 0.1713244923791704]),
)

"""
    add_volume!(net, el; V, side = :main, n = 5)

Attach a gas volume `V` [m³] to element `el` (for a `HeatExchanger`, give
`side = :hot` or `:cold`; add both sides separately). Its mass is the density
integrated along linear inlet→outlet temperature and pressure profiles with
`n` Gauss points (1–6). Volumes make the loop inventory a model quantity —
see [`loop_mass`](@ref) and [`set_inventory!`](@ref). Where two components
join directly, put the connecting pipe's volume on either of them or on a
zero-loss `Duct`.
"""
function add_volume!(net::FlowNetwork, el; V::Real, side::Symbol = :main, n::Int = 5)
    any(e -> e === el, net.elements) || error("add_volume!: \"$(el.name)\" is not in the network")
    haskey(_GL, n) || error("add_volume!: n must be 1–6")
    if el isa HeatExchanger
        side in (:hot, :cold) || error("add_volume!: give side = :hot or :cold for heat exchanger \"$(el.name)\"")
    else
        side === :main || error("add_volume!: side applies to heat exchangers only")
    end
    any(v -> v.el === el && v.side === side, net.volumes) &&
        error("add_volume!: \"$(el.name)\" ($side) already has a volume")
    push!(net.volumes, VolumeSpec(el, side, V, n))
    net
end

# Inlet and outlet ports of a volume.
function _volume_ports(v::VolumeSpec)
    el = v.el
    v.side === :hot  && return (el.hot_inlet, el.hot_outlet)
    v.side === :cold && return (el.cold_inlet, el.cold_outlet)
    el isa Mixer    && return (el.outlet, el.outlet)
    el isa Splitter && return (el.inlet, el.inlet)
    (el.inlet, el.outlet)
end

"""    volume_mass(v::VolumeSpec) -> gas mass in the volume [kg]."""
function volume_mass(v::VolumeSpec)
    pin, pout = _volume_ports(v)
    (isnothing(pin) || isnothing(pout)) && error("volume of \"$(v.el.name)\": solve the network first")
    a, b = pin[], pout[]
    x, w = _GL[v.n]
    ρ̄ = sum(w[k] * density(a.fluid, a.Tt + (b.Tt - a.Tt) * x[k], a.Pt + (b.Pt - a.Pt) * x[k])
            for k in eachindex(x))
    v.V * ρ̄
end

# Stream-component key of a volume.
function _volume_key(net::FlowNetwork, comps, v::VolumeSpec)
    i = findfirst(e -> e === v.el, net.elements)
    port = v.side === :hot ? :hot_inlet : v.side === :cold ? :cold_inlet :
           v.el isa Mixer ? :outlet : :inlet
    get(comps, (i, port), (i, port))
end

function _el_key(net::FlowNetwork, comps, el; side = :main)
    i = findfirst(e -> e === el, net.elements)
    isnothing(i) && error("\"$(el.name)\" is not in the network")
    port = el isa HeatExchanger ? (side === :hot ? :hot_inlet : :cold_inlet) :
           el isa Mixer ? :outlet : :inlet
    get(comps, (i, port), (i, port))
end

"""
    loop_mass(net, el; side = :main) -> Real

Total gas mass [kg] in the volumes of the loop (connected stream) that passes
through element `el` (`side = :hot`/`:cold` picks a heat-exchanger side).
"""
function loop_mass(net::FlowNetwork, el; side::Symbol = :main)
    comps = _stream_components(net)
    k = _el_key(net, comps, el; side)
    sum(volume_mass(v) for v in net.volumes if _volume_key(net, comps, v) == k; init = 0.0)
end

"""    total_mass(net) -> Σ mass in all volumes [kg]."""
total_mass(net::FlowNetwork) = sum(volume_mass, net.volumes; init = 0.0)

# The seed (network seed or boundary guess) of the loop through `el`, with
# get/set closures for one of its fields.
function _loop_seed(net::FlowNetwork, el, field::Symbol; side = :main)
    comps = _stream_components(net)
    k = _el_key(net, comps, el; side)
    if net.seed_el > 0 && !isnothing(net.seed_state)
        s_el = net.elements[net.seed_el]
        if get(comps, (net.seed_el, _canon_inlet(s_el, :inlet)), nothing) == k
            return (() -> getfield(net.seed_state, field),
                    x -> (net.seed_state = update(net.seed_state; (field => x,)...); nothing))
        end
    end
    for (j, (i, p, st)) in enumerate(net.boundaries)
        if get(comps, (i, _canon_inlet(net.elements[i], p)), nothing) == k
            return (() -> getfield(net.boundaries[j][3], field),
                    x -> (b = net.boundaries[j];
                          net.boundaries[j] = (b[1], b[2], update(b[3]; (field => x,)...)); nothing))
        end
    end
    error("no seed or boundary state found for the loop through \"$(el.name)\"")
end

"""
    set_inventory!(net, el; mass, vary = :W, side = :main, name)

Inventory mode for the loop through `el`: the solver varies that loop's
circulating mass flow (`vary = :W`, the loop's seed or boundary-guess state)
so the gas mass in the loop's volumes equals `mass` [kg]. Use it where the
pressure level floats — off-design with maps, where the hardware fixes the
corrected flows — to make charge (inventory) the control input instead of
flow, as on a real closed Brayton unit. `vary = :Pt` varies the seed pressure
instead (for models where flow is specified and pressure must follow).
Needs volumes (`add_volume!`) on the loop's components.
"""
function set_inventory!(net::FlowNetwork, el; mass::Real, vary::Symbol = :W,
                        side::Symbol = :main, name::AbstractString = "$(el.name) loop inventory")
    vary in (:W, :Pt) || error("set_inventory!: vary must be :W or :Pt")
    get, set! = _loop_seed(net, el, vary; side)
    add_free!(net, "$(el.name) loop $vary"; get = get, set! = set!, lower = 0.0)
    add_target!(net, n -> loop_mass(n, el; side), mass; name = name)
end

"""    inventory_report([io,] net) — per-volume V, mean density and mass, with loop totals."""
inventory_report(net::FlowNetwork) = inventory_report(stdout, net)
function inventory_report(io::IO, net::FlowNetwork)
    isempty(net.volumes) && return println(io, "  (no volumes)")
    comps = _stream_components(net)
    keys_ = unique(_volume_key(net, comps, v) for v in net.volumes)
    f(x, d) = round(_primal(x), digits = d)
    for (li, k) in enumerate(keys_)
        vs = [v for v in net.volumes if _volume_key(net, comps, v) == k]
        println(io, "  loop $li")
        tot = 0.0
        for v in vs
            m = volume_mass(v)
            tot += _primal(m)
            lbl = v.side === :main ? v.el.name : "$(v.el.name) ($(v.side))"
            println(io, "    ", rpad(lbl, 28), lpad("V=$(f(v.V, 5)) m³", 16),
                    lpad("ρ̄=$(f(m / v.V, 4)) kg/m³", 22), lpad("m=$(f(m, 5)) kg", 16))
        end
        println(io, "    ", rpad("total", 28), lpad("", 38), lpad("m=$(round(tot, digits = 5)) kg", 16))
    end
    nothing
end
