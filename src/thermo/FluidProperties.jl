"""
Abstract interface for fluid thermodynamic property backends.

All cycle components call these methods exclusively — they never
reference a specific backend. Swap FPTFluid for IdealGasFluid (or any
other implementation) without touching component code.

All inputs:  T [K],  P [Pa]
All outputs: SI units (J, kg, K, Pa)

Methods accept any numeric type (Real subtypes) so that ForwardDiff Dual
numbers propagate through compute! for automatic differentiation.
Backends that override T_from_h / T_from_s with exact closed-form
inversions (like IdealGasFluid) fully support AD; backends using the
default bisection fallback return the primal root but drop derivatives.
"""
abstract type FluidProperties end

"""
    cp(fp, T, P) -> specific heat [J/(kg·K)]
"""
function cp(fp::FluidProperties, T, P)
    error("cp not implemented for $(typeof(fp))")
end

"""
    enthalpy(fp, T, P) -> specific enthalpy [J/kg], referenced to 0 K
"""
function enthalpy(fp::FluidProperties, T, P)
    error("enthalpy not implemented for $(typeof(fp))")
end

"""
    entropy(fp, T, P) -> specific entropy [J/(kg·K)]
"""
function entropy(fp::FluidProperties, T, P)
    error("entropy not implemented for $(typeof(fp))")
end

"""
    density(fp, T, P) -> mass density [kg/m³]
"""
function density(fp::FluidProperties, T, P)
    error("density not implemented for $(typeof(fp))")
end

"""
    gamma(fp, T, P) -> isentropic exponent Cp/Cv [-]
"""
function gamma(fp::FluidProperties, T, P)
    error("gamma not implemented for $(typeof(fp))")
end

"""
    viscosity(fp, T, P) -> dynamic viscosity μ [Pa·s]
"""
function viscosity(fp::FluidProperties, T, P)
    error("viscosity not implemented for $(typeof(fp))")
end

"""
    conductivity(fp, T, P) -> thermal conductivity k [W/(m·K)]
"""
function conductivity(fp::FluidProperties, T, P)
    error("conductivity not implemented for $(typeof(fp))")
end

"""
    prandtl(fp, T, P) -> Prandtl number cp·μ/k [-]
"""
prandtl(fp::FluidProperties, T, P) =
    cp(fp, T, P) * viscosity(fp, T, P) / conductivity(fp, T, P)

"""
    T_from_h(fp, h_target, P; T_guess) -> T [K]

Invert enthalpy: find T such that enthalpy(fp, T, P) ≈ h_target.
Default implementation uses bisection; backends may override for speed and AD support.
"""
function T_from_h(fp::FluidProperties, h_target, P; T_guess=500.0)
    T_lo, T_hi = 100.0, 5000.0
    for _ in 1:60
        T_mid = 0.5 * (T_lo + T_hi)
        enthalpy(fp, T_mid, P) < h_target ? (T_lo = T_mid) : (T_hi = T_mid)
        (T_hi - T_lo) < 1e-6 && break
    end
    0.5 * (T_lo + T_hi)
end

"""
    T_from_s(fp, s_target, P; T_guess) -> T [K]

Invert entropy: find T such that entropy(fp, T, P) ≈ s_target.
Used for isentropic process calculations.
"""
function T_from_s(fp::FluidProperties, s_target, P; T_guess=500.0)
    T_lo, T_hi = 100.0, 5000.0
    for _ in 1:60
        T_mid = 0.5 * (T_lo + T_hi)
        entropy(fp, T_mid, P) < s_target ? (T_lo = T_mid) : (T_hi = T_mid)
        (T_hi - T_lo) < 1e-6 && break
    end
    0.5 * (T_lo + T_hi)
end
