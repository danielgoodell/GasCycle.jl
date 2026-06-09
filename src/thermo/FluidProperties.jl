"""
Abstract interface for fluid thermodynamic property backends.

All cycle components call these methods exclusively — they never
reference a specific backend. Swap FPTFluid for IdealGasFluid (or any
other implementation) without touching component code.

All inputs:  T [K],  P [Pa]
All outputs: SI units (J, kg, K, Pa)
"""
abstract type FluidProperties end

"""
    cp(fp, T, P) -> Float64

Specific heat at constant pressure [J/(kg·K)].
"""
function cp(fp::FluidProperties, T::Float64, P::Float64)::Float64
    error("cp not implemented for $(typeof(fp))")
end

"""
    enthalpy(fp, T, P) -> Float64

Specific enthalpy [J/kg], referenced to 0 K.
"""
function enthalpy(fp::FluidProperties, T::Float64, P::Float64)::Float64
    error("enthalpy not implemented for $(typeof(fp))")
end

"""
    entropy(fp, T, P) -> Float64

Specific entropy [J/(kg·K)].
"""
function entropy(fp::FluidProperties, T::Float64, P::Float64)::Float64
    error("entropy not implemented for $(typeof(fp))")
end

"""
    density(fp, T, P) -> Float64

Mass density [kg/m³].
"""
function density(fp::FluidProperties, T::Float64, P::Float64)::Float64
    error("density not implemented for $(typeof(fp))")
end

"""
    gamma(fp, T, P) -> Float64

Isentropic exponent Cp/Cv [-].
"""
function gamma(fp::FluidProperties, T::Float64, P::Float64)::Float64
    error("gamma not implemented for $(typeof(fp))")
end

"""
    T_from_h(fp, h_target, P; T_guess) -> Float64

Invert enthalpy: find T such that enthalpy(fp, T, P) ≈ h_target.
Default implementation uses bisection; backends may override for speed.
"""
function T_from_h(fp::FluidProperties, h_target::Float64, P::Float64;
                  T_guess::Float64=500.0)::Float64
    # Bisection over a wide range — backends may override with Newton
    T_lo, T_hi = 100.0, 5000.0
    for _ in 1:60
        T_mid = 0.5 * (T_lo + T_hi)
        if enthalpy(fp, T_mid, P) < h_target
            T_lo = T_mid
        else
            T_hi = T_mid
        end
        (T_hi - T_lo) < 1e-6 && break
    end
    0.5 * (T_lo + T_hi)
end

"""
    T_from_s(fp, s_target, P; T_guess) -> Float64

Invert entropy: find T such that entropy(fp, T, P) ≈ s_target.
Used for isentropic process calculations.
"""
function T_from_s(fp::FluidProperties, s_target::Float64, P::Float64;
                  T_guess::Float64=500.0)::Float64
    T_lo, T_hi = 100.0, 5000.0
    for _ in 1:60
        T_mid = 0.5 * (T_lo + T_hi)
        if entropy(fp, T_mid, P) < s_target
            T_lo = T_mid
        else
            T_hi = T_mid
        end
        (T_hi - T_lo) < 1e-6 && break
    end
    0.5 * (T_lo + T_hi)
end
