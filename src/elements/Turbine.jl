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
"""
mutable struct Turbine <: AbstractElement
    name::String
    PR::Float64          # expansion ratio (Pt_in / Pt_out)
    η_poly::Float64      # polytropic efficiency
    map::Union{PerformanceMap, Nothing}
    mode::Symbol
    N_shaft::Float64     # [rpm] — set by Shaft
    P_exit::Float64      # target exit pressure for :pressure_closure [Pa]

    inlet::Union{Port, Nothing}
    outlet::Union{Port, Nothing}
    Wc_map::Float64
end

function Turbine(name::String;
                 PR::Float64     = 2.0,
                 η_poly::Float64 = 0.90,
                 map::Union{PerformanceMap,Nothing} = nothing,
                 mode::Symbol    = :design,
                 P_exit::Float64 = 101325.0)
    Turbine(name, PR, η_poly, map, mode, 0.0, P_exit, nothing, nothing, 0.0)
end

function compute!(el::Turbine, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]

    if el.mode == :pressure_closure
        # PR set so turbine exhausts exactly to the specified exit pressure
        el.PR = s.Pt / el.P_exit
    elseif el.mode == :off_design && !isnothing(el.map)
        Nc = corrected_speed(el.N_shaft, s.Tt)
        Wc_act = corrected_flow(s.W, s.Tt, s.Pt)
        PR_map, η_map = query(el.map, Nc, Wc_act)
        el.PR    = PR_map
        el.η_poly = η_map
        el.Wc_map = Wc_act
    end

    fp = s.fluid
    h_in = enthalpy(fp, s.Tt, s.Pt)
    s_in = entropy(fp, s.Tt, s.Pt)

    Pt_out = s.Pt / el.PR

    Tt_is  = T_from_s(fp, s_in, Pt_out; T_guess = s.Tt / el.PR^0.3)
    h_is   = enthalpy(fp, Tt_is, Pt_out)
    h_out  = h_in - (h_in - h_is) * el.η_poly

    Tt_out = T_from_h(fp, h_out, Pt_out; T_guess = Tt_is * 0.95)

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

# Only :off_design has solver unknowns (shaft speed from map).
# :design and :pressure_closure have fixed or auto-computed PR.
function indep_vars(el::Turbine)
    el.mode == :off_design && return [el.N_shaft]
    Float64[]
end

function set_indep_vars!(el::Turbine, x::AbstractVector)
    el.mode == :off_design && (el.N_shaft = x[1])
end

specific_work(el::Turbine) =
    enthalpy(el.inlet[]) - enthalpy(el.outlet[])
