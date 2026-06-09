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

The type parameter T allows design variables PR and η_poly to carry
ForwardDiff Dual numbers for gradient-based design optimization.
"""
mutable struct Compressor{T<:Real} <: AbstractElement
    name::String
    PR::T            # total pressure ratio
    η_poly::T        # polytropic efficiency (≈ adiabatic for small stages)
    map::Union{PerformanceMap, Nothing}
    mode::Symbol     # :design or :off_design
    N_shaft::Float64 # shaft speed [rpm] — set by Shaft each iteration

    # State (written by compute!, read by residuals)
    inlet::Union{Port, Nothing}
    outlet::Union{Port, Nothing}
    Wc_map::Float64  # corrected flow from map (off-design only)
end

function Compressor(name::String;
                    PR         = 2.0,
                    η_poly      = 0.87,
                    map::Union{PerformanceMap,Nothing} = nothing,
                    mode::Symbol = :design)
    T = promote_type(typeof(PR), typeof(η_poly))
    Compressor{T}(name, T(PR), T(η_poly), map, mode, 0.0, nothing, nothing, 0.0)
end

function compute!(el::Compressor, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]

    if el.mode == :off_design && !isnothing(el.map)
        Nc = corrected_speed(el.N_shaft, s.Tt)
        Wc_act = corrected_flow(s.W, s.Tt, s.Pt)
        PR_map, η_map = query(el.map, Nc, Wc_act)
        el.PR     = PR_map
        el.η_poly = η_map
        el.Wc_map = Wc_act
    end

    fp     = s.fluid
    Pt_out = s.Pt * el.PR
    Tt_out = _polytropic_outlet(fp, s.Tt, s.Pt, Pt_out, el.η_poly)

    outlet_state = update(s; Pt = Pt_out, Tt = Tt_out)
    el.outlet = Port(outlet_state)
    el.outlet
end

n_residuals(el::Compressor) = el.mode == :off_design ? 1 : 0

function residuals(el::Compressor)
    el.mode == :design && return Float64[]
    s = el.inlet[]
    Wc_act = corrected_flow(s.W, s.Tt, s.Pt)
    [el.Wc_map - Wc_act]
end

indep_vars(el::Compressor) = el.mode == :off_design ? [el.N_shaft] : Float64[]
function set_indep_vars!(el::Compressor, x::AbstractVector)
    el.mode == :off_design && (el.N_shaft = x[1])
end

"""Work input per unit mass flow [J/kg]"""
specific_work(el::Compressor) =
    enthalpy(el.outlet[]) - enthalpy(el.inlet[])
