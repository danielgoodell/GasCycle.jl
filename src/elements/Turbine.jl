"""
Turbine element — mirror of Compressor.

Expands flow from inlet Pt to outlet Pt = Pt_in / PR, extracting enthalpy
to drive the shaft.  Polytropic efficiency < 1 means less work is extracted
than the isentropic ideal.

Modes:
  :design           — PR is user-specified (fixed parameter).
  :pressure_closure — PR is auto-computed so the turbine exhausts to P_exit
                      (the compressor inlet pressure for a closed cycle).
                      Provide `P_exit` at construction.
  :off_design       — PR and η from the performance map; shaft speed is the
                      independent variable.

The type parameter T allows design variables PR, η_poly, and P_exit to carry
ForwardDiff Dual numbers for gradient-based design optimization.
"""
mutable struct Turbine{T<:Real} <: AbstractElement
    name::String
    PR::T            # expansion ratio (Pt_in / Pt_out)
    η_poly::T        # polytropic efficiency
    map::Union{PerformanceMap, Nothing}
    mode::Symbol
    N_shaft::Float64 # [rpm] — set by Shaft
    P_exit::T        # target exit pressure for :pressure_closure [Pa]

    inlet::Union{Port, Nothing}
    outlet::Union{Port, Nothing}
    Wc_map::Float64
end

function Turbine(name::String;
                 PR      = 2.0,
                 η_poly   = 0.90,
                 map::Union{PerformanceMap,Nothing} = nothing,
                 mode::Symbol = :design,
                 P_exit   = 101325.0)
    T = promote_type(typeof(PR), typeof(η_poly), typeof(P_exit))
    Turbine{T}(name, T(PR), T(η_poly), map, mode, 0.0, T(P_exit), nothing, nothing, 0.0)
end

function compute!(el::Turbine, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]

    # PR_eff is computed locally so Dual numbers can flow through pressure_closure
    # without trying to store a Dual value into el.PR (which may be Float64).
    PR_eff = if el.mode == :pressure_closure
        s.Pt / el.P_exit   # local only; el.PR left unchanged
    elseif el.mode == :off_design && !isnothing(el.map)
        Nc = corrected_speed(el.N_shaft, s.Tt)
        Wc_act = corrected_flow(s.W, s.Tt, s.Pt)
        PR_map, η_map = query(el.map, Nc, Wc_act)
        el.PR    = PR_map
        el.η_poly = η_map
        el.Wc_map = Wc_act
        el.PR
    else
        el.PR
    end

    fp     = s.fluid
    Pt_out = s.Pt / PR_eff
    Tt_out = _polytropic_outlet(fp, s.Tt, s.Pt, Pt_out, el.η_poly)

    outlet_state = update(s; Pt = Pt_out, Tt = Tt_out)
    el.outlet = Port(outlet_state)
    el.outlet
end

n_residuals(el::Turbine) = el.mode == :off_design ? 1 : 0

function residuals(el::Turbine)
    el.mode == :off_design || return Float64[]
    s = el.inlet[]
    Wc_act = corrected_flow(s.W, s.Tt, s.Pt)
    [el.Wc_map - Wc_act]
end

function indep_vars(el::Turbine)
    el.mode == :off_design && return [el.N_shaft]
    Float64[]
end

function set_indep_vars!(el::Turbine, x::AbstractVector)
    el.mode == :off_design && (el.N_shaft = x[1])
end

specific_work(el::Turbine) =
    enthalpy(el.inlet[]) - enthalpy(el.outlet[])
