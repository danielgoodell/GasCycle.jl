"""
Duct element — flow with a total pressure loss, and a heat port.

  Pt_out = Pt_in * (1 - dPqP)
  h_out  = h_in + Q_ext / W   (Q_ext = net heat from attached heat links;
                               0 → isenthalpic, Tt_out = Tt_in)
  W_out  = W_in

A duct is also the stream endpoint of heat links (`add_conductance!`,
`add_radiation!`, `add_loss_heat!`): a crossover duct, a wall segment, or a
zero-loss segment after a component whose casing leaks heat.

No residual equations; no independent variables.
"""
mutable struct Duct{T<:Real} <: AbstractElement
    name::String
    dPqP::T  # fractional total pressure loss (ΔP/P)
    heated::Bool     # true once a heat link attaches
    Q_ext::Real      # net heat into the stream from heat links [W]
    inlet::Union{Port, Nothing}
    outlet::Union{Port, Nothing}
end

Duct(name::String; dPqP = 0.02) =
    Duct{typeof(dPqP)}(name, dPqP, false, 0.0, nothing, nothing)

function compute!(el::Duct, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]
    Pt_out = s.Pt * (1 - el.dPqP)
    outlet_state = if el.heated
        h_out = enthalpy(s) + el.Q_ext / s.W
        update(s; Pt = Pt_out, Tt = T_from_h(s.fluid, h_out, Pt_out; T_guess = s.Tt))
    else
        update(s; Pt = Pt_out)
    end
    el.outlet = Port(outlet_state)
    el.outlet
end

n_residuals(el::Duct)  = 0
residuals(el::Duct)    = Float64[]
indep_vars(el::Duct)   = Float64[]
set_indep_vars!(el::Duct, x::AbstractVector) = nothing
