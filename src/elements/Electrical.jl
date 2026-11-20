"""
Electrical power-flow elements.

The electrical side is a radial tree of nodes connected with
`connect_electrical!`: a source (an `Alternator`) at the root, then any chain
of cables, converters and buses, ending in loads. It is solved inside every
`one_pass!` before the gas path (the alternator's losses heat the gas), so it
joins the cycle's Newton solve and stays differentiable.

Modeling level: steady power flow. DC nodes carry (V, I, P); AC nodes carry
three-phase line-to-line RMS voltage, frequency and real power at unity power
factor (no reactive power or harmonics). Voltages propagate downstream from the
source and converter set points; power demand propagates upstream from the
loads; resistive voltage drops couple the two and are converged by a short
fixed-point iteration (a couple of sweeps without cables).

Node interface (implemented below and by `Alternator`):

    e_kind_in(el)   :ac, :dc, or nothing (source)
    e_kind_out(el)  :ac, :dc, or nothing (load)
    _e_forward!(el, V_in, f_in) -> (V_out, f_out)      voltage sweep
    _e_backward!(el, P_out)     -> P_in                 power sweep

All numeric fields are `Real` so any of them can be a free variable and can
carry an outer gradient's Duals.
"""

abstract type AbstractElectrical end

"""Electrical operating point of a node, written by the electrical solve."""
mutable struct ElecPoint
    V_in::Real     # input voltage [V] (DC, or AC line-line RMS)
    V_out::Real    # output voltage [V]
    f::Real        # output frequency [Hz] (0 for DC)
    I::Real        # current in the element [A] (DC, or AC line current)
    P_in::Real     # power in [W]
    P_out::Real    # power out [W]
    P_loss::Real   # dissipated [W]
end
ElecPoint() = ElecPoint(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# Current carried at power P and voltage V.
_current(kind::Symbol, P, V) = kind === :dc ? P / V : P / (sqrt(3) * V)

# ── Cable ─────────────────────────────────────────────────────────────────────

"""
    Cable(name; kind = :dc, R = nothing, length, area, resistivity = 1.72e-8,
          temp_coeff = 0.00393, T = 293.15)

Resistive transmission line — the electrical "duct". `R` is the resistance of
one conductor [Ω]; give it directly, or give `length` [m] and conductor
cross-section `area` [m²] and it is computed from `resistivity` (copper at
20 °C by default) corrected to conductor temperature `T` [K]:
R = ρ₂₀·(1 + α·(T − 293.15))·length/area.

DC: two conductors, loss 2I²R, drop 2IR. AC: three-phase, loss 3I²R, line-line
drop √3·I·R (resistive, unity power factor). Exact for the given power: the
current solves P_out = V_in·I − loss(I). Demanding more than the line can
deliver gives NaN (the solver treats it as infeasible).
"""
mutable struct Cable <: AbstractElectrical
    name::String
    kind::Symbol
    R::Real
    e::ElecPoint
end

function Cable(name::String; kind::Symbol = :dc, R = nothing, length = nothing,
               area = nothing, resistivity = 1.72e-8, temp_coeff = 0.00393, T = 293.15)
    kind in (:ac, :dc) || error("Cable \"$name\": kind must be :ac or :dc")
    if isnothing(R)
        (isnothing(length) || isnothing(area)) &&
            error("Cable \"$name\": give R, or length and area")
        R = resistivity * (1 + temp_coeff * (T - 293.15)) * length / area
    end
    Cable(name, kind, R, ElecPoint())
end

e_kind_in(el::Cable)  = el.kind
e_kind_out(el::Cable) = el.kind

function _e_forward!(el::Cable, V_in, f_in)
    drop = el.kind === :dc ? 2 * el.e.I * el.R : sqrt(3) * el.e.I * el.R
    el.e.V_in, el.e.V_out, el.e.f = V_in, V_in - drop, f_in
    (el.e.V_out, f_in)
end

function _e_backward!(el::Cable, P_out)
    V, R = el.e.V_in, el.R
    # DC:  P_out = V I − 2R I²        AC:  P_out = √3 V I − 3R I²
    a, b = el.kind === :dc ? (2R, V) : (3R, sqrt(3) * V)
    I = if iszero(R)
        P_out / b
    else
        disc = b^2 - 4a * P_out
        _primal(disc) < 0 ? oftype(disc, NaN) : (b - sqrt(disc)) / (2a)
    end
    loss = a * I^2
    el.e.I, el.e.P_out, el.e.P_loss, el.e.P_in = I, P_out, loss, P_out + loss
    el.e.P_in
end

# ── Converters (rectifier, DC-DC, transformer, generic box) ──────────────────

"""
    Converter(name; input, output, V_out = nothing, ratio = nothing,
              η = nothing, P_noload = 0, k_lin = 0, k_quad = 0)

Power-conversion stage. Output voltage is either regulated (`V_out` [V]) or a
fixed `ratio` of the input voltage. Losses as a function of output power:

    P_loss = P_noload + k_lin·P_out + k_quad·P_out²

(standby + switching + conduction). `η` is shorthand for a constant
efficiency: k_lin = 1/η − 1. AC→AC stages keep the input frequency.

Use the named constructors: [`Rectifier`](@ref), [`DCDCConverter`](@ref),
[`Transformer`](@ref), [`PowerConverter`](@ref).
"""
mutable struct Converter <: AbstractElectrical
    name::String
    kind_in::Symbol
    kind_out::Symbol
    V_set::Real       # regulated output voltage, or NaN
    ratio::Real       # V_out / V_in when not regulated, or NaN
    P_noload::Real
    k_lin::Real
    k_quad::Real
    e::ElecPoint
end

function Converter(name::String; input::Symbol, output::Symbol, V_out = nothing,
                   ratio = nothing, η = nothing, P_noload = 0.0, k_lin = 0.0, k_quad = 0.0)
    input in (:ac, :dc) && output in (:ac, :dc) ||
        error("Converter \"$name\": input/output must be :ac or :dc")
    isnothing(V_out) == isnothing(ratio) &&
        error("Converter \"$name\": give exactly one of V_out (regulated) or ratio")
    if !isnothing(η)
        0 < η <= 1 || error("Converter \"$name\": η must be in (0, 1]")
        k_lin = 1 / η - 1
    end
    Converter(name, input, output, something(V_out, NaN), something(ratio, NaN),
              P_noload, k_lin, k_quad, ElecPoint())
end

"""
    Rectifier(name; V_out = nothing, ratio = 1.35, kw...)

AC→DC. Default is a passive three-phase diode bridge, V_dc = 1.35·V_LL; give
`V_out` for an active (regulated) rectifier. Loss keywords as in `Converter`.
"""
Rectifier(name::String; V_out = nothing, ratio = nothing, kw...) =
    Converter(name; input = :ac, output = :dc, V_out = V_out,
              ratio = isnothing(V_out) ? something(ratio, 1.35) : nothing, kw...)

"""    DCDCConverter(name; V_out, kw...)   — regulated DC→DC stage."""
DCDCConverter(name::String; V_out, kw...) =
    Converter(name; input = :dc, output = :dc, V_out = V_out, kw...)

"""    Transformer(name; ratio, kw...)   — AC→AC, V_out = ratio·V_in, same frequency."""
Transformer(name::String; ratio, kw...) =
    Converter(name; input = :ac, output = :ac, ratio = ratio, kw...)

"""
    PowerConverter(name; input = :ac, output = :dc, V_out, kw...)

Generic regulated conversion box — e.g. "AC in, DC out at bus voltage" when the
internal topology (rectifier, DC link, DC-DC) is not being modeled.
"""
PowerConverter(name::String; input::Symbol = :ac, output::Symbol = :dc, V_out, kw...) =
    Converter(name; input = input, output = output, V_out = V_out, kw...)

e_kind_in(el::Converter)  = el.kind_in
e_kind_out(el::Converter) = el.kind_out

function _e_forward!(el::Converter, V_in, f_in)
    V_out = isnan(_primal(el.V_set)) ? el.ratio * V_in : el.V_set
    f_out = el.kind_out === :ac ? f_in : zero(f_in)
    el.e.V_in, el.e.V_out, el.e.f = V_in, V_out, f_out
    (V_out, f_out)
end

function _e_backward!(el::Converter, P_out)
    loss = el.P_noload + el.k_lin * P_out + el.k_quad * P_out^2
    el.e.P_out, el.e.P_loss, el.e.P_in = P_out, loss, P_out + loss
    el.e.I = _current(el.kind_in, el.e.P_in, el.e.V_in)
    el.e.P_in
end

# ── Bus ───────────────────────────────────────────────────────────────────────

"""
    ElectricalBus(name; kind = :dc)

Lossless junction (a DC link or distribution bus): one feed, any number of
branches at the same voltage. Its power out is the sum of its branches.
"""
mutable struct ElectricalBus <: AbstractElectrical
    name::String
    kind::Symbol
    e::ElecPoint
end
ElectricalBus(name::String; kind::Symbol = :dc) = ElectricalBus(name, kind, ElecPoint())

e_kind_in(el::ElectricalBus)  = el.kind
e_kind_out(el::ElectricalBus) = el.kind

function _e_forward!(el::ElectricalBus, V_in, f_in)
    el.e.V_in, el.e.V_out, el.e.f = V_in, V_in, f_in
    (V_in, f_in)
end

function _e_backward!(el::ElectricalBus, P_out)
    el.e.P_out, el.e.P_in, el.e.P_loss = P_out, P_out, zero(P_out)
    el.e.I = _current(el.kind, P_out, el.e.V_in)
    P_out
end

# ── Loads ─────────────────────────────────────────────────────────────────────

"""
    ElectricalLoad(name; P, kind = :dc, parasitic = false)

Power sink drawing `P` [W]. A parasitic load is power dumped on purpose — see
[`ParasiticLoad`](@ref). `P` can be a free variable (`add_free!(net, load, :P)`).
"""
mutable struct ElectricalLoad <: AbstractElectrical
    name::String
    kind::Symbol
    P::Real
    parasitic::Bool
    e::ElecPoint
end
ElectricalLoad(name::String; P, kind::Symbol = :dc, parasitic::Bool = false) =
    ElectricalLoad(name, kind, P, parasitic, ElecPoint())

"""
    ParasiticLoad(name; P = 0.0, kind = :dc)

Parasitic load radiator (PLR): dumps surplus electrical power as heat to space,
the standard speed control of a space Brayton unit. Pair it with a speed
target — [`add_speed_control!`](@ref) — so the solver sets `P`.
"""
ParasiticLoad(name::String; P = 0.0, kind::Symbol = :dc) =
    ElectricalLoad(name; P = P, kind = kind, parasitic = true)

e_kind_in(el::ElectricalLoad)  = el.kind
e_kind_out(::ElectricalLoad)   = nothing

function _e_forward!(el::ElectricalLoad, V_in, f_in)
    el.e.V_in, el.e.V_out, el.e.f = V_in, V_in, f_in
    (V_in, f_in)
end

function _e_backward!(el::ElectricalLoad, P_out)
    el.e.P_in, el.e.P_out, el.e.P_loss = el.P, zero(el.P), zero(el.P)
    el.e.I = _current(el.kind, el.P, el.e.V_in)
    el.P
end
