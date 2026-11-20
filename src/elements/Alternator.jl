"""
    Alternator(name; V_rated, N_rated, poles = 2, η = nothing,
               R_phase = 0.0, P_core_rated = 0.0, core_exp = 1.5,
               P_windage_rated = 0.0, dPqP = 0.0)

Shaft-driven permanent-magnet alternator, cooled by the main gas flow: it sits
in the gas path (`connect!(net, comp => alt => recup …)`) and every watt of
loss heats the gas passing through it. It is also the source of the
electrical network (`connect_electrical!(net, alt => cable => …)`) and a load
on its shaft (`add_shaft!(net, shaft; drives = [comp, alt], driven_by = turb)`).

Electrical output: three-phase AC at unity power factor.
- frequency  f = poles·N/120 [Hz]
- open-circuit line-line voltage ∝ speed: V_oc = V_rated·N/N_rated
- terminal voltage V = V_oc − √3·I·R_phase

Losses (all to the gas), either a constant efficiency `η`
(P_loss = P_elec·(1/η − 1)) or the loss model
- copper   3·I²·R_phase
- core     P_core_rated·(f/f_rated)^core_exp
- windage  P_windage_rated·(N/N_rated)³

Shaft power drawn: P_mech = P_elec + P_loss. `dPqP` is the gas-side pressure
loss.
"""
mutable struct Alternator <: AbstractElement
    name::String
    η::Real                  # constant efficiency, or NaN for the loss model
    V_rated::Real            # line-line RMS open-circuit voltage at N_rated [V]
    N_rated::Real            # [rpm]
    poles::Int
    R_phase::Real            # [Ω]
    P_core_rated::Real       # [W] at rated frequency
    core_exp::Real
    P_windage_rated::Real    # [W] at N_rated
    dPqP::Real
    N_shaft::Real            # [rpm] — set by Shaft
    P_mech::Real             # shaft power drawn [W]
    e::ElecPoint
    inlet::Union{Port,Nothing}
    outlet::Union{Port,Nothing}
end

function Alternator(name::String; V_rated, N_rated, poles::Int = 2, η = nothing,
                    R_phase = 0.0, P_core_rated = 0.0, core_exp = 1.5,
                    P_windage_rated = 0.0, dPqP = 0.0)
    isnothing(η) || 0 < η <= 1 || error("Alternator \"$name\": η must be in (0, 1]")
    Alternator(name, something(η, NaN), V_rated, N_rated, poles, R_phase, P_core_rated,
               core_exp, P_windage_rated, dPqP, N_rated, 0.0, ElecPoint(),
               nothing, nothing)
end

# ── Electrical source ────────────────────────────────────────────────────────

e_kind_in(::Alternator)  = nothing
e_kind_out(::Alternator) = :ac

electrical_frequency(el::Alternator) = el.poles * el.N_shaft / 120

function _e_forward!(el::Alternator, _V, _f)
    V_oc = el.V_rated * el.N_shaft / el.N_rated
    V = V_oc - sqrt(3) * el.e.I * el.R_phase
    el.e.V_in, el.e.V_out, el.e.f = V_oc, V, electrical_frequency(el)
    (V, el.e.f)
end

function _e_backward!(el::Alternator, P_out)
    I = _current(:ac, P_out, el.e.V_out)
    loss = if !isnan(_primal(el.η))
        P_out * (1 / el.η - 1)
    else
        f_rated = el.poles * el.N_rated / 120
        n = el.N_shaft / el.N_rated
        3 * I^2 * el.R_phase + el.P_core_rated * (el.e.f / f_rated)^el.core_exp +
            el.P_windage_rated * n^3
    end
    el.e.I, el.e.P_out, el.e.P_loss = I, P_out, loss
    el.e.P_in = el.P_mech = P_out + loss
    el.P_mech
end

# ── Gas side: the loss heats the coolant gas ─────────────────────────────────

function compute!(el::Alternator, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]
    Pt_out = s.Pt * (1 - el.dPqP)
    h_out  = enthalpy(s) + el.e.P_loss / s.W
    Tt_out = T_from_h(s.fluid, h_out, Pt_out; T_guess = s.Tt)
    el.outlet = Port(update(s; Pt = Pt_out, Tt = Tt_out))
    el.outlet
end

n_residuals(::Alternator) = 0
residuals(::Alternator)   = Float64[]
indep_vars(::Alternator)  = Float64[]
set_indep_vars!(::Alternator, x::AbstractVector) = nothing

"""Shaft power drawn by the alternator [W]."""
shaft_power(el::Alternator) = el.P_mech
