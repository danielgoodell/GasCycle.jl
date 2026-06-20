"""
Compressor element.

Raises total pressure by PR and total temperature by the isentropic work
divided by polytropic efficiency.  Supports two operating modes:

  :design   — PR is a free parameter set by the user.  Residual is empty
              (0 equations); the design point is fully specified.  A
              TurbomachineMap is scaled (via `scale_map`) to pass through this
              design point before being attached for off-design.

  :off_design — PR and efficiency come from the performance map evaluated
               FORWARD in its native coordinates at (Nc, Rline), where Nc
               follows from the shaft speed and Rline — the map R-line
               coordinate — is the independent variable owned by the solver.
               The map *returns* the corrected flow Wc it would pass; the
               residual is flow continuity Wc_pred - actual_Wc = 0 (one
               equation, one unknown).  No Wc inversion, so the formulation is
               well-conditioned near choke.  Shaft speed belongs to the Shaft.

An optional `reynolds::ReynoldsModel` supplies the abstract index fed to the
map's Reynolds-correction tables; the default `nothing` disables the correction.

Efficiency semantics are selected by `η_type`:
  :polytropic (default) — η_poly is the polytropic (small-stage) efficiency,
                          applied by stepwise integration.
  :isentropic           — η_poly is the adiabatic (isentropic) efficiency,
                          applied to the full Δh_is in one step.  This matches
                          NPSS `effDes`/map `eff` semantics.

The type parameter T allows design variables PR and η_poly to carry
ForwardDiff Dual numbers for gradient-based design optimization.
"""
mutable struct Compressor{T<:Real} <: AbstractElement
    name::String
    PR::T            # total pressure ratio
    η_poly::T        # efficiency value; meaning set by η_type
    η_type::Symbol   # :polytropic or :isentropic
    map::Union{TurbomachineMap, Nothing}
    reynolds::Union{ReynoldsModel, Nothing}  # map Reynolds-index model (nothing ⇒ no correction)
    mode::Symbol     # :design or :off_design
    N_shaft::Float64 # shaft speed [rpm] — set by Shaft each iteration

    # State (written by compute!, read by residuals)
    inlet::Union{Port, Nothing}
    outlet::Union{Port, Nothing}
    Rline::Float64   # map R-line coordinate — solver independent (off-design only)
end

function Compressor(name::String;
                    PR         = 2.0,
                    η_poly      = 0.87,
                    η_type::Symbol = :polytropic,
                    map::Union{TurbomachineMap,Nothing} = nothing,
                    reynolds::Union{ReynoldsModel,Nothing} = nothing,
                    mode::Symbol = :design)
    η_type in (:polytropic, :isentropic) ||
        error("Compressor \"$name\": η_type must be :polytropic or :isentropic, got :$η_type")
    T = promote_type(typeof(PR), typeof(η_poly))
    # Seed R-line from the map's design anchor (a valid mid-map coordinate);
    # 0.0 means "not yet seeded" and is filled on the first off-design pass.
    Compressor{T}(name, T(PR), T(η_poly), η_type, map, reynolds, mode, 0.0,
                  nothing, nothing, 0.0)
end

function compute!(el::Compressor, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]

    if el.mode == :off_design && !isnothing(el.map)
        Nc = corrected_speed(el.N_shaft, s.Tt)
        # Seed the R-line coordinate from the map design anchor on the first
        # pass, before the solver has taken ownership of it.
        el.Rline > 0.0 || (el.Rline = design_line(el.map))
        rc = re_coord(el.reynolds, s)
        out = eval_map(el.map, Nc, el.Rline, rc)
        el.PR     = out.PR
        el.η_poly = out.eff
    end

    fp     = s.fluid
    Pt_out = s.Pt * el.PR
    Tt_out = _efficiency_outlet(el.η_type, fp, s.Tt, s.Pt, Pt_out, el.η_poly)

    outlet_state = update(s; Pt = Pt_out, Tt = Tt_out)
    el.outlet = Port(outlet_state)
    el.outlet
end

n_residuals(el::Compressor) = el.mode == :off_design ? 1 : 0

# Flow continuity: the corrected flow the map passes at (Nc, Rline) must equal
# the actual corrected flow coming in.  Rline is the unknown that closes it.
function residuals(el::Compressor)
    el.mode == :off_design || return Float64[]
    s = el.inlet[]
    Nc = corrected_speed(el.N_shaft, s.Tt)
    rc = re_coord(el.reynolds, s)
    Wc_pred = eval_map(el.map, Nc, el.Rline, rc).Wc
    Wc_act  = corrected_flow(s.W, s.Tt, s.Pt)
    [(Wc_pred - Wc_act) / Wc_pred]
end

indep_vars(el::Compressor) = el.mode == :off_design ? [el.Rline] : Float64[]
function set_indep_vars!(el::Compressor, x::AbstractVector)
    el.mode == :off_design && (el.Rline = x[1])
end

"""Work input per unit mass flow [J/kg]"""
specific_work(el::Compressor) =
    enthalpy(el.outlet[]) - enthalpy(el.inlet[])

"""
    pressure_ratio(el::Compressor) -> Real

Actual pressure ratio Pt_out / Pt_in from the last computed state.
Falls back to `el.PR` if the element has not been computed yet.
"""
pressure_ratio(el::Compressor) =
    (isnothing(el.inlet) || isnothing(el.outlet)) ? el.PR :
    el.outlet[].Pt / el.inlet[].Pt
