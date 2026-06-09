"""
HeatSource element — external heat addition (reactor core, solar receiver, etc.).

Adds heat Q [W] to the working fluid.  No change in composition or mass flow.
A small total pressure loss dPqP models flow friction through the heat source.

Operating modes:
  :fixed_Q      — Q is given; Tt_out is computed from energy balance.
  :fixed_TtExit — target Tt_out is given; Q is computed from energy balance.
"""
mutable struct HeatSource{T<:Real} <: AbstractElement
    name::String
    Q::T         # heat addition [W]  (positive = heat in)
    TtExit::T    # target exit total temperature [K]  (:fixed_TtExit mode)
    dPqP::T      # fractional total pressure loss
    mode::Symbol # :fixed_Q or :fixed_TtExit

    inlet::Union{Port, Nothing}
    outlet::Union{Port, Nothing}
end

function HeatSource(name::String;
                    Q       = 0.0,
                    TtExit  = 1200.0,
                    dPqP    = 0.02,
                    mode::Symbol = :fixed_TtExit)
    T = promote_type(typeof(Q), typeof(TtExit), typeof(dPqP))
    HeatSource{T}(name, T(Q), T(TtExit), T(dPqP), mode, nothing, nothing)
end

function compute!(el::HeatSource, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]
    fp = s.fluid

    Pt_out = s.Pt * (1 - el.dPqP)

    if el.mode == :fixed_TtExit
        Tt_out = el.TtExit
        h_out  = enthalpy(fp, Tt_out, Pt_out)
        h_in   = enthalpy(fp, s.Tt,   s.Pt)
        Q_val  = (h_out - h_in) * s.W
        # Only store if type matches; in AD context Q_val may be Dual while el.Q is Float64
        Q_val isa typeof(el.Q) && (el.Q = Q_val)
    else  # :fixed_Q
        h_in   = enthalpy(fp, s.Tt, s.Pt)
        h_out  = h_in + el.Q / s.W
        Tt_out = T_from_h(fp, h_out, Pt_out; T_guess = s.Tt * 1.5)
    end

    outlet_state = update(s; Pt = Pt_out, Tt = Tt_out)
    el.outlet = Port(outlet_state)
    el.outlet
end

n_residuals(el::HeatSource)          = 0
residuals(el::HeatSource)            = Float64[]
indep_vars(el::HeatSource)           = Float64[]
set_indep_vars!(el::HeatSource, x::AbstractVector) = nothing
