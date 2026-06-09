"""
Immutable snapshot of thermodynamic state at a cycle station.

Derived properties (cp, h, s, γ, ρ) are computed on demand from the
fluid backend — they are never stored, so the struct stays small and
ForwardDiff-friendly.

The type parameter T allows Dual numbers to flow through: use
FluidState{Float64} for normal operation, FluidState{Dual{...}} when
differentiating cycle performance with respect to design parameters.
"""
struct FluidState{T<:Real}
    Pt::T              # total pressure     [Pa]
    Tt::T              # total temperature  [K]
    W::T               # mass flow rate     [kg/s]
    fluid::FluidProperties  # thermodynamic backend (FPTFluid, IdealGasFluid, …)

    # Inner constructor: all three numeric fields must have the same type T
    FluidState{T}(Pt::T, Tt::T, W::T, fluid::FluidProperties) where T<:Real =
        new{T}(Pt, Tt, W, fluid)
end

# Outer promoting constructor: accepts mixed numeric types, promotes to a common T
function FluidState(Pt, Tt, W, fluid::FluidProperties)
    T = promote_type(typeof(Pt), typeof(Tt), typeof(W))
    FluidState{T}(T(Pt), T(Tt), T(W), fluid)
end

# Derived property accessors — delegate to the fluid backend
cp(s::FluidState)       = cp(s.fluid, s.Tt, s.Pt)
enthalpy(s::FluidState) = enthalpy(s.fluid, s.Tt, s.Pt)
entropy(s::FluidState)  = entropy(s.fluid, s.Tt, s.Pt)
density(s::FluidState)  = density(s.fluid, s.Tt, s.Pt)
gamma(s::FluidState)    = gamma(s.fluid, s.Tt, s.Pt)

"""Return a new FluidState with updated fields; unspecified fields copied from s."""
function update(s::FluidState; Pt=s.Pt, Tt=s.Tt, W=s.W)
    FluidState(Pt, Tt, W, s.fluid)
end
