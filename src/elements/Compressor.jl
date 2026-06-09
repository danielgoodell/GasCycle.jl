"""
Compressor element.

Raises total pressure by PR and total temperature by the isentropic work
divided by polytropic efficiency.  Supports two operating modes:

  :design   — PR is a free parameter set by the user.  Residual is empty
              (0 equations); the design point is fully specified.
              When a PerformanceMap is attached, scale factors are solved
              so the map passes through this design point.

  :off_design — PR and efficiency come from the performance map interpolated
               at (Nc, Wc) from the current shaft speed and inlet conditions.
               Residual: map_Wc - actual_Wc = 0  (one equation).
"""
mutable struct Compressor <: AbstractElement
    name::String
    PR::Float64          # total pressure ratio
    η_poly::Float64      # polytropic efficiency (≈ adiabatic for small stages)
    map::Union{PerformanceMap, Nothing}
    mode::Symbol         # :design or :off_design
    N_shaft::Float64     # shaft speed [rpm] — set by Shaft each iteration

    # State (written by compute!, read by residuals)
    inlet::Union{Port, Nothing}
    outlet::Union{Port, Nothing}
    Wc_map::Float64      # corrected flow from map (off-design only)
end

function Compressor(name::String;
                    PR::Float64    = 2.0,
                    η_poly::Float64 = 0.87,
                    map::Union{PerformanceMap,Nothing} = nothing,
                    mode::Symbol   = :design)
    Compressor(name, PR, η_poly, map, mode, 0.0, nothing, nothing, 0.0)
end

function compute!(el::Compressor, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]

    if el.mode == :off_design && !isnothing(el.map)
        Nc = corrected_speed(el.N_shaft, s.Tt)
        Wc_act = corrected_flow(s.W, s.Tt, s.Pt)
        # Query map: get PR and η at this Nc; Wc_act is the actual corrected flow
        # We use the map's Nc axis at the stored Wc from last iteration
        PR_map, η_map = query(el.map, Nc, Wc_act)
        el.PR   = PR_map
        el.η_poly = η_map
        el.Wc_map = Wc_act
    end

    fp = s.fluid
    h_in  = enthalpy(fp, s.Tt, s.Pt)
    s_in  = entropy(fp, s.Tt, s.Pt)

    Pt_out = s.Pt * el.PR

    # Isentropic exit temperature
    Tt_is = T_from_s(fp, s_in, Pt_out; T_guess = s.Tt * el.PR^0.3)

    # Actual exit enthalpy via polytropic efficiency
    h_is_out = enthalpy(fp, Tt_is, Pt_out)
    h_out    = h_in + (h_is_out - h_in) / el.η_poly

    Tt_out = T_from_h(fp, h_out, Pt_out; T_guess = Tt_is * 1.1)

    outlet_state = update(s; Pt = Pt_out, Tt = Tt_out)
    el.outlet = Port(outlet_state)
    el.outlet
end

n_residuals(el::Compressor) = el.mode == :off_design ? 1 : 0

function residuals(el::Compressor)
    el.mode == :design && return Float64[]
    # Off-design: map corrected flow must equal actual corrected flow
    s = el.inlet[]
    Wc_act = corrected_flow(s.W, s.Tt, s.Pt)
    [el.Wc_map - Wc_act]
end

# In :design mode PR is a fixed user parameter, not a solver unknown.
# In :off_design mode shaft speed is the unknown (PR comes from the map).
indep_vars(el::Compressor) = el.mode == :off_design ? [el.N_shaft] : Float64[]
function set_indep_vars!(el::Compressor, x::AbstractVector)
    el.mode == :off_design && (el.N_shaft = x[1])
end

"""Work input per unit mass flow [W per kg/s = J/kg]"""
specific_work(el::Compressor) =
    enthalpy(el.outlet[]) - enthalpy(el.inlet[])
