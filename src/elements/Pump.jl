"""
    Pump(name; ΔP, η = 0.5, kind = :dc)

Liquid coolant pump with a canned motor: raises total pressure by `ΔP` [Pa]
and draws electrical power

    P_elec = W·ΔP / (ρ·η)

where `η` is the overall wire-to-water efficiency and ρ the inlet density.
All of that power ends up in the coolant (hydraulic work + motor and
impeller losses), so h_out = h_in + P_elec / W.

Wire it to a DC (or AC, `kind = :ac`) bus with `connect_electrical!(net, bus
=> pump)` to draw from the cycle's electrical output; unwired, its power is
counted as an external input. `P_elec` is a solver unknown (the electrical
network is solved before the coolant flow is known in each pass), so a pump
puts the model on the coupled Newton path.
"""
mutable struct Pump <: AbstractElement
    name::String
    ΔP::Real
    η::Real
    kind::Symbol
    P_elec::Real     # electrical draw [W] — solver unknown
    scale::Float64   # residual normalization [W]
    e::ElecPoint
    inlet::Union{Port,Nothing}
    outlet::Union{Port,Nothing}
end

function Pump(name::String; ΔP, η = 0.5, kind::Symbol = :dc)
    0 < η <= 1 || error("Pump \"$name\": η must be in (0, 1]")
    Pump(name, ΔP, η, kind, 0.0, 1.0, ElecPoint(), nothing, nothing)
end

"""Electrical power the pump needs at its current inlet state [W]."""
function pump_power(el::Pump)
    s = el.inlet[]
    s.W * el.ΔP / (density(s) * el.η)
end

function compute!(el::Pump, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]
    Pt_out = s.Pt + el.ΔP
    h_out  = enthalpy(s) + el.P_elec / s.W
    el.outlet = Port(update(s; Pt = Pt_out, Tt = T_from_h(s.fluid, h_out, Pt_out; T_guess = s.Tt)))
    el.outlet
end

n_residuals(::Pump) = 1
residuals(el::Pump) = [(el.P_elec - pump_power(el)) / el.scale]
indep_vars(el::Pump) = [el.P_elec]
set_indep_vars!(el::Pump, x::AbstractVector) = (el.P_elec = x[1]; nothing)

# Electrical load interface
e_kind_in(el::Pump)  = el.kind
e_kind_out(::Pump)   = nothing
function _e_forward!(el::Pump, V_in, f_in)
    el.e.V_in, el.e.V_out, el.e.f = V_in, V_in, f_in
    (V_in, f_in)
end
function _e_backward!(el::Pump, _P_out)
    el.e.P_in, el.e.P_out, el.e.P_loss = el.P_elec, zero(el.P_elec), zero(el.P_elec)
    el.e.I = _current(el.kind, el.P_elec, el.e.V_in)
    el.P_elec
end
