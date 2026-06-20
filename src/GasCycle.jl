module GasCycle

# ── Units ─────────────────────────────────────────────────────────────────────
include("units.jl")

# ── Thermodynamics ────────────────────────────────────────────────────────────
include("thermo/FluidProperties.jl")
include("thermo/IdealGasFluid.jl")
include("thermo/FPTFluid.jl")
include("thermo/ConstantPropertyLiquid.jl")
include("thermo/PolynomialLiquid.jl")
include("thermo/NobleGasMixture.jl")

# ── Core abstractions ─────────────────────────────────────────────────────────
include("core/FluidState.jl")
include("core/Port.jl")
include("core/Element.jl")

# ── Performance maps ──────────────────────────────────────────────────────────
include("maps/TurbomachineMap.jl")
include("maps/NPSSMapReader.jl")

# ── Elements ──────────────────────────────────────────────────────────────────
include("elements/Compressor.jl")
include("elements/Turbine.jl")
include("elements/Duct.jl")
include("elements/Shaft.jl")
include("elements/HeatSource.jl")
include("elements/Radiator.jl")
include("elements/HeatExchanger.jl")
include("elements/Splitter.jl")
include("elements/Mixer.jl")

# ── Network & Solver ──────────────────────────────────────────────────────────
include("network/FlowNetwork.jl")
include("solver/Solver.jl")

# ── Output ────────────────────────────────────────────────────────────────────
include("output/Summary.jl")
include("output/PlotRecipes.jl")

# ── Public API ────────────────────────────────────────────────────────────────
export R_to_K, K_to_R, psia_to_Pa, Pa_to_psia
export lbm_to_kg, kg_to_lbm, lbps_to_kgps, kgps_to_lbps
export btulbm_to_Jkg, Jkg_to_btulbm, btulbmR_to_JkgK, JkgK_to_btulbmR
export lbmft3_to_kgm3, kgm3_to_lbmft3, rpm_to_radps, radps_to_rpm
export hp_to_W, W_to_hp

export FluidProperties, IdealGasFluid, HeXeIdealGas, FPTFluid, ConstantPropertyLiquid, PolynomialLiquid
export NobleGas, NobleGasMixture, NobleGasFluid, HeXe
export HELIUM, NEON, ARGON, KRYPTON, XENON
# Note: `cp` is not exported to avoid conflict with Base.Filesystem.cp in Julia ≥ 1.12.
# Use GasCycle.cp(...) or `import GasCycle: cp` to access it.
export enthalpy, entropy, density, gamma, T_from_h, T_from_s, h_from_s
export viscosity, conductivity, prandtl

export FluidState, Port, AbstractElement
export update

export TurbomachineMap, CompressorMap, TurbineMap, FunctionMap
export eval_map, scale_map, corrected_speed, corrected_flow
export ReynoldsModel, ReDesIndex, RawRe, FunctionReynolds, re_coord, reynolds_index
export read_npss_map, compressor_map, turbine_map

export Compressor, Turbine, Duct, Shaft, HeatSource, HeatExchanger, Radiator, Splitter, Mixer
export compute!, compute_hx!, size_UA!, specific_work, pressure_ratio, Q_transferred, Q_rejected, power_balance
export n_residuals, residuals, indep_vars, set_indep_vars!
export link!

export FlowNetwork, add!, connect!, connect_port!, add_shaft!, add_hx_pair!, set_state!, set_boundary!, one_pass!
export solve!, cycle_efficiency, net_power, SolveResult
export stations   # summary(sol) extends Base.summary — no export needed
export tsdiagram, tsdiagram!, mapplot, mapplot!   # RecipesBase user plots

end # module GasCycle
