"""
    Valve(name; A, Cd = 0.8, opening = 1.0)

Throttling valve: an orifice of flow area `A·opening` [m²] with discharge
coefficient `Cd`. Total-pressure loss

    ΔPt = W² / (2·ρ·(Cd·A·opening)²)

with ρ at the inlet (incompressible orifice — valid while ΔPt ≪ Pt), at
constant total enthalpy. `opening` (0–1] is `Real`, so it can be a free
variable: a bypass valve around the turbine (Splitter → Valve → Mixer)
controlled to hold speed, or a throttle in a coolant loop. A demand the valve
cannot pass (ΔPt ≥ Pt) gives NaN, which the solver treats as infeasible.
"""
mutable struct Valve <: AbstractElement
    name::String
    A::Real
    Cd::Real
    opening::Real
    inlet::Union{Port,Nothing}
    outlet::Union{Port,Nothing}
end

Valve(name::String; A, Cd = 0.8, opening = 1.0) = Valve(name, A, Cd, opening, nothing, nothing)

"""Total-pressure drop across the valve at its current inlet state [Pa]."""
function valve_dP(el::Valve, s::FluidState)
    A_eff = el.Cd * el.A * el.opening
    s.W^2 / (2 * density(s) * A_eff^2)
end

function compute!(el::Valve, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]
    Pt_out = s.Pt - valve_dP(el, s)
    _primal(Pt_out) > 0 || (Pt_out = oftype(Pt_out, NaN))
    el.outlet = Port(update(s; Pt = Pt_out))
    el.outlet
end

n_residuals(::Valve) = 0
residuals(::Valve)   = Float64[]
indep_vars(::Valve)  = Float64[]
set_indep_vars!(::Valve, x::AbstractVector) = nothing
