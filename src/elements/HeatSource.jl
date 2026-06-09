"""
HeatSource element — external heat addition (reactor core, solar receiver, etc.).

Adds heat Q [W] to the working fluid.  No change in composition or mass flow.
A small total pressure loss dPqP models flow friction through the heat source.

Operating modes:
  :fixed_Q      — Q is given; Tt_out is computed from energy balance.
  :fixed_TtExit — target Tt_out is given; Q is the independent variable
                  (residual: Tt_out_computed - TtExit_target = 0).
"""
mutable struct HeatSource <: AbstractElement
    name::String
    Q::Float64           # heat addition [W]  (positive = heat in)
    TtExit::Float64      # target exit total temperature [K]  (:fixed_TtExit mode)
    dPqP::Float64        # fractional total pressure loss
    mode::Symbol         # :fixed_Q or :fixed_TtExit

    inlet::Union{Port, Nothing}
    outlet::Union{Port, Nothing}
end

function HeatSource(name::String;
                    Q::Float64       = 0.0,
                    TtExit::Float64  = 1200.0,
                    dPqP::Float64    = 0.02,
                    mode::Symbol     = :fixed_TtExit)
    HeatSource(name, Q, TtExit, dPqP, mode, nothing, nothing)
end

function compute!(el::HeatSource, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]
    fp = s.fluid

    Pt_out = s.Pt * (1.0 - el.dPqP)

    if el.mode == :fixed_TtExit
        Tt_out = el.TtExit
        h_out  = enthalpy(fp, Tt_out, Pt_out)
        h_in   = enthalpy(fp, s.Tt,   s.Pt)
        el.Q   = (h_out - h_in) * s.W
    else  # :fixed_Q
        h_in   = enthalpy(fp, s.Tt, s.Pt)
        h_out  = h_in + el.Q / s.W
        Tt_out = T_from_h(fp, h_out, Pt_out; T_guess = s.Tt * 1.5)
    end

    outlet_state = update(s; Pt = Pt_out, Tt = Tt_out)
    el.outlet = Port(outlet_state)
    el.outlet
end

# compute! handles both modes directly — no solver participation needed.
n_residuals(el::HeatSource)          = 0
residuals(el::HeatSource)            = Float64[]
indep_vars(el::HeatSource)           = Float64[]
set_indep_vars!(el::HeatSource, x::AbstractVector) = nothing
