"""
Immutable snapshot of thermodynamic state at a cycle station.

Derived properties (cp, h, s, γ, ρ) are computed on demand from the
fluid backend — they are never stored, so the struct stays small and
ForwardDiff-friendly.
"""
struct FluidState
    Pt::Float64            # total pressure     [Pa]
    Tt::Float64            # total temperature  [K]
    W::Float64             # mass flow rate     [kg/s]
    fluid::FluidProperties # thermodynamic backend (FPTFluid, IdealGasFluid, …)

    # Inner constructor handles automatic promotion from other numeric types
    FluidState(Pt, Tt, W, fluid) = new(Float64(Pt), Float64(Tt), Float64(W), fluid)
end

# Derived property accessors — delegate to the fluid backend
cp(s::FluidState)      = cp(s.fluid, s.Tt, s.Pt)
enthalpy(s::FluidState) = enthalpy(s.fluid, s.Tt, s.Pt)
entropy(s::FluidState)  = entropy(s.fluid, s.Tt, s.Pt)
density(s::FluidState)  = density(s.fluid, s.Tt, s.Pt)
gamma(s::FluidState)    = gamma(s.fluid, s.Tt, s.Pt)

"""Return a new FluidState with updated total pressure and temperature."""
function update(s::FluidState; Pt=s.Pt, Tt=s.Tt, W=s.W)
    FluidState(Pt, Tt, W, s.fluid)
end
