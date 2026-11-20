"""
Thermal links: heat exchanged between streams, environment nodes and
electrical components outside the main flow elements.

Endpoints
- a **stream**: a `Duct` (its heat port). The stream presents the mean of the
  duct's inlet and outlet temperatures (lumped wall) — accurate when the link
  shifts the stream temperature only slightly; split the duct into segments,
  or use a `HeatExchanger` in `:UA` mode, when it does not.
- an **environment node**: [`ThermalBoundary`](@ref), a named equivalent sink
  temperature (e.g. one per spacecraft-model node, fed from a thermal model).
- an **electrical component**: its loss heat (`add_loss_heat!`).

Each link's heat flow `Q` [W] (positive from `a` to `b`) is an unknown of the
cycle's Newton solve with the residual Q − Q_model(T_a, T_b) = 0, so links can
couple any parts of the loop, stay differentiable, and can be targeted.
"""

"""
    ThermalBoundary(name; T)

Environment node at a fixed (equivalent sink) temperature `T` [K] — space,
structure, a cabin, or a node of an external spacecraft thermal model. Use
several for different surroundings. `T` is `Real`, so it can be a free
variable or carry Duals; update many at once with [`set_environment!`](@ref).
Radiators can reject to one: `Radiator(...; sink = node)`.
"""
mutable struct ThermalBoundary
    name::String
    T::Real
end
ThermalBoundary(name::String; T) = ThermalBoundary(name, T)

const _σ_SB = 5.670374419e-8

"""
A heat link between two endpoints (see `add_conductance!`, `add_radiation!`,
`add_loss_heat!`). `Q` [W] flows from `a` to `b`.
"""
mutable struct HeatLink
    name::String
    kind::Symbol          # :conductance, :radiation, :loss
    a::Any                # Duct, ThermalBoundary, or electrical node (:loss)
    b::Any                # Duct or ThermalBoundary
    G::Real               # UA [W/K] (:conductance), ε·A·F [m²] (:radiation), fraction (:loss)
    Q::Real               # heat flow a → b [W] — solver unknown
    scale::Float64        # residual normalization [W]
end

_endpoint_T(d) = d isa ThermalBoundary ? d.T :
    (isnothing(d.inlet) || isnothing(d.outlet)) ? error("heat link endpoint \"$(d.name)\" not computed") :
    (d.inlet[].Tt + d.outlet[].Tt) / 2

function heat_model(l::HeatLink)
    if l.kind === :conductance
        l.G * (_endpoint_T(l.a) - _endpoint_T(l.b))
    elseif l.kind === :radiation
        _σ_SB * l.G * (_endpoint_T(l.a)^4 - _endpoint_T(l.b)^4)
    else  # :loss — a share of an electrical component's dissipation
        l.G * l.a.e.P_loss
    end
end

# Newton interface (collected alongside element unknowns).
n_residuals(::HeatLink) = 1
residuals(l::HeatLink)  = [(l.Q - heat_model(l)) / l.scale]
indep_vars(l::HeatLink) = [l.Q]
set_indep_vars!(l::HeatLink, x::AbstractVector) = (l.Q = x[1]; nothing)
