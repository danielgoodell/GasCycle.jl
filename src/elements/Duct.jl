"""
Duct element — constant-enthalpy flow with a total pressure loss.

  Pt_out = Pt_in * (1 - dPqP)
  Tt_out = Tt_in            (isenthalpic, no heat transfer)
  W_out  = W_in

No residual equations; no independent variables.
"""
mutable struct Duct <: AbstractElement
    name::String
    dPqP::Float64   # fractional total pressure loss (ΔP/P)
    inlet::Union{Port, Nothing}
    outlet::Union{Port, Nothing}
end

Duct(name::String; dPqP::Float64 = 0.02) = Duct(name, dPqP, nothing, nothing)

function compute!(el::Duct, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]
    outlet_state = update(s; Pt = s.Pt * (1.0 - el.dPqP))
    el.outlet = Port(outlet_state)
    el.outlet
end

n_residuals(el::Duct)  = 0
residuals(el::Duct)    = Float64[]
indep_vars(el::Duct)   = Float64[]
set_indep_vars!(el::Duct, x) = nothing
