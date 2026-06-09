"""
Ideal-gas backend for monatomic fluid mixtures (e.g., He-Xe).

Cp is constant = (5/2) * R_specific for a monatomic ideal gas.
Equation of state: P = ρ R T.
Entropy includes both temperature and pressure terms.

Use this for unit tests and examples that run without an FPT file.
All methods have closed-form expressions and fully support ForwardDiff AD.
"""
struct IdealGasFluid <: FluidProperties
    R::Float64      # specific gas constant [J/(kg·K)]
    cp_val::Float64 # constant Cp = (5/2) R  [J/(kg·K)]
    h_ref::Float64  # enthalpy at T_ref=0 K (typically 0)

    function IdealGasFluid(; M_molar::Float64)
        R_universal = 8314.46261815324  # J/(kmol·K)
        R = R_universal / M_molar       # J/(kg·K)
        cp_val = 2.5 * R                # monatomic ideal gas
        new(R, cp_val, 0.0)
    end
end

"""
    HeXeIdealGas(x_He) -> IdealGasFluid

Construct an ideal-gas He-Xe mixture by helium mole fraction x_He.
"""
function HeXeIdealGas(x_He::Float64)
    M_He = 4.002602    # kg/kmol
    M_Xe = 131.293
    M_mix = x_He * M_He + (1.0 - x_He) * M_Xe
    IdealGasFluid(M_molar=M_mix)
end

cp(fp::IdealGasFluid, T, P)       = fp.cp_val
enthalpy(fp::IdealGasFluid, T, P) = fp.cp_val * T + fp.h_ref
entropy(fp::IdealGasFluid, T, P)  =
    fp.cp_val * log(T / 298.15) - fp.R * log(P / 101325.0)
density(fp::IdealGasFluid, T, P)  = P / (fp.R * T)
gamma(fp::IdealGasFluid, T, P)    = fp.cp_val / (fp.cp_val - fp.R)

# Exact closed-form inversions — fully AD-compatible
T_from_h(fp::IdealGasFluid, h_target, P; T_guess=500.0) =
    (h_target - fp.h_ref) / fp.cp_val

T_from_s(fp::IdealGasFluid, s_target, P; T_guess=500.0) =
    298.15 * exp((s_target + fp.R * log(P / 101325.0)) / fp.cp_val)
