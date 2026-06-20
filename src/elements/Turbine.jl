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
  :off_design       — the map is PR-parameterized (PR is its native *input*), so
                      the expansion ratio PR is the independent variable owned by
                      the solver.  The map evaluated forward at (Np, PR) returns
                      the corrected flow Wp it passes and the efficiency; the
                      residual is flow continuity Wp - actual_Wp = 0.  Pressure
                      closure of the loop is enforced by the network back-edge,
                      not here.  Shaft speed belongs to the Shaft element.

An optional `reynolds::ReynoldsModel` supplies the abstract index fed to the
map's Reynolds-correction tables; the default `nothing` disables the correction.

Efficiency semantics are selected by `η_type` (see Compressor): :polytropic
(default, small-stage integration) or :isentropic (adiabatic efficiency on
the full Δh_is, matching NPSS `effDes`).

The type parameter T allows design variables PR, η_poly, and P_exit to carry
ForwardDiff Dual numbers for gradient-based design optimization.
"""
mutable struct Turbine{T<:Real} <: AbstractElement
    name::String
    PR::T            # expansion ratio (Pt_in / Pt_out); off-design solver independent
    η_poly::T        # efficiency value; meaning set by η_type
    η_type::Symbol   # :polytropic or :isentropic
    map::Union{TurbomachineMap, Nothing}
    reynolds::Union{ReynoldsModel, Nothing}  # map Reynolds-index model (nothing ⇒ no correction)
    mode::Symbol
    N_shaft::Float64 # [rpm] — set by Shaft
    P_exit::T        # target exit pressure for :pressure_closure [Pa]

    inlet::Union{Port, Nothing}
    outlet::Union{Port, Nothing}
end

function Turbine(name::String;
                 PR      = 2.0,
                 η_poly   = 0.90,
                 η_type::Symbol = :polytropic,
                 map::Union{TurbomachineMap,Nothing} = nothing,
                 reynolds::Union{ReynoldsModel,Nothing} = nothing,
                 mode::Symbol = :design,
                 P_exit   = 101325.0)
    η_type in (:polytropic, :isentropic) ||
        error("Turbine \"$name\": η_type must be :polytropic or :isentropic, got :$η_type")
    T = promote_type(typeof(PR), typeof(η_poly), typeof(P_exit))
    Turbine{T}(name, T(PR), T(η_poly), η_type, map, reynolds, mode, 0.0, T(P_exit),
               nothing, nothing)
end

function compute!(el::Turbine, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]

    # PR_eff is computed locally so Dual numbers can flow through pressure_closure
    # without trying to store a Dual value into el.PR (which may be Float64).
    PR_eff = if el.mode == :pressure_closure
        s.Pt / el.P_exit   # local only; el.PR left unchanged
    elseif el.mode == :off_design && !isnothing(el.map)
        # PR is the solver independent (the map's native input); read η off the
        # map at the current (Np, PR).
        Np = corrected_speed(el.N_shaft, s.Tt)
        rc = re_coord(el.reynolds, s)
        el.η_poly = eval_map(el.map, Np, el.PR, rc).eff
        el.PR
    else
        el.PR
    end

    fp     = s.fluid
    Pt_out = s.Pt / PR_eff
    Tt_out = _efficiency_outlet(el.η_type, fp, s.Tt, s.Pt, Pt_out, el.η_poly)

    outlet_state = update(s; Pt = Pt_out, Tt = Tt_out)
    el.outlet = Port(outlet_state)
    el.outlet
end

n_residuals(el::Turbine) = el.mode == :off_design ? 1 : 0

# Flow continuity: the corrected flow the map passes at (Np, PR) must equal the
# actual corrected flow.  PR is the unknown that closes it; loop pressure closure
# is handled by the network back-edge.
function residuals(el::Turbine)
    el.mode == :off_design || return Float64[]
    s = el.inlet[]
    Np = corrected_speed(el.N_shaft, s.Tt)
    rc = re_coord(el.reynolds, s)
    Wp_pred = eval_map(el.map, Np, el.PR, rc).Wp
    Wp_act  = corrected_flow(s.W, s.Tt, s.Pt)
    [(Wp_pred - Wp_act) / Wp_pred]
end

indep_vars(el::Turbine) = el.mode == :off_design ? [el.PR] : Float64[]
function set_indep_vars!(el::Turbine, x::AbstractVector)
    el.mode == :off_design && (el.PR = x[1])
end

specific_work(el::Turbine) =
    enthalpy(el.inlet[]) - enthalpy(el.outlet[])

"""
    pressure_ratio(el::Turbine) -> Real

Actual expansion ratio Pt_in / Pt_out from the last computed state.  Unlike
the `el.PR` field, this is correct in :pressure_closure mode, where the
effective PR is computed locally and `el.PR` is left unchanged.  Falls back
to `el.PR` if the element has not been computed yet.
"""
pressure_ratio(el::Turbine) =
    (isnothing(el.inlet) || isnothing(el.outlet)) ? el.PR :
    el.inlet[].Pt / el.outlet[].Pt
