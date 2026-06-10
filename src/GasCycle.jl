module GasCycle

# ── Units ─────────────────────────────────────────────────────────────────────
include("units.jl")

# ── Thermodynamics ────────────────────────────────────────────────────────────
include("thermo/FluidProperties.jl")
include("thermo/IdealGasFluid.jl")
include("thermo/FPTFluid.jl")

# ── Core abstractions ─────────────────────────────────────────────────────────
include("core/FluidState.jl")
include("core/Port.jl")
include("core/Element.jl")

# ── Performance maps ──────────────────────────────────────────────────────────
include("maps/PerformanceMap.jl")
include("maps/MapScaling.jl")

# ── Elements ──────────────────────────────────────────────────────────────────
include("elements/Compressor.jl")
include("elements/Turbine.jl")
include("elements/Duct.jl")
include("elements/Shaft.jl")
include("elements/HeatSource.jl")
include("elements/HeatExchanger.jl")
include("elements/Splitter.jl")
include("elements/Mixer.jl")

# ── Network & Solver ──────────────────────────────────────────────────────────
include("network/FlowNetwork.jl")
include("solver/Solver.jl")

# ── Public API ────────────────────────────────────────────────────────────────
export R_to_K, K_to_R, psia_to_Pa, Pa_to_psia
export lbm_to_kg, kg_to_lbm, lbps_to_kgps, kgps_to_lbps
export btulbm_to_Jkg, Jkg_to_btulbm, btulbmR_to_JkgK, JkgK_to_btulbmR
export lbmft3_to_kgm3, kgm3_to_lbmft3, rpm_to_radps, radps_to_rpm
export hp_to_W, W_to_hp

export FluidProperties, IdealGasFluid, HeXeIdealGas, FPTFluid
# Note: `cp` is not exported to avoid conflict with Base.Filesystem.cp in Julia ≥ 1.12.
# Use GasCycle.cp(...) or `import GasCycle: cp` to access it.
export enthalpy, entropy, density, gamma, T_from_h, T_from_s, h_from_s

export FluidState, Port, AbstractElement
export update

export PerformanceMap, scale_map, query, corrected_speed, corrected_flow

export Compressor, Turbine, Duct, Shaft, HeatSource, HeatExchanger, Splitter, Mixer
export compute!, compute_hx!, specific_work, pressure_ratio, Q_transferred, power_balance
export n_residuals, residuals, indep_vars, set_indep_vars!
export link!

export FlowNetwork, add!, connect!, connect_port!, add_shaft!, add_hx_pair!, set_state!, one_pass!
export solve!, cycle_efficiency, net_power, SolveResult

end # module GasCycle
