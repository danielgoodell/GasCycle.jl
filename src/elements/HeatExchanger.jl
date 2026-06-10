"""
Counter-flow heat exchanger / recuperator.

Models a two-stream heat exchanger using the ε-NTU method.  The two streams
(hot and cold) must be declared as a pair and linked in the FlowNetwork via
`add_hx_pair!` (recuperator) or `set_boundary!` (coolant side of a
heat-rejection cooler).

  Q = ε * Q_max
  Q_max is computed from the limiting endpoint enthalpy change.

Effectiveness comes from one of two modes:

  :effectiveness — ε is a fixed parameter (default).
  :UA            — ε is computed each pass from the counter-flow ε-NTU
                   relation with NTU = UA / C_min, so effectiveness
                   responds to off-design flow and property changes.
                   Selected automatically when the `UA` keyword is given.

Both streams conserve mass flow; only Tt changes (pressure loss dPqP applied
to each stream independently).

Usage in FlowNetwork:
  hx = HeatExchanger("Recup"; ε=0.95)              # fixed effectiveness
  hx = HeatExchanger("Cooler"; UA=1200.0)          # ε from UA [W/K]
  add_hx_pair!(net, hx, hot_inlet_port, cold_inlet_port)
"""
mutable struct HeatExchanger{T<:Real} <: AbstractElement
    name::String
    ε::T           # effectiveness (0–1); recomputed each pass in :UA mode
    UA::T          # overall conductance [W/K] (:UA mode)
    mode::Symbol   # :effectiveness or :UA
    dPqP_hot::T
    dPqP_cold::T

    # Set by add_hx_pair! before first compute!
    hot_inlet::Union{Port, Nothing}
    cold_inlet::Union{Port, Nothing}
    hot_outlet::Union{Port, Nothing}
    cold_outlet::Union{Port, Nothing}
end

function HeatExchanger(name::String;
                       ε         = 0.90,
                       UA        = nothing,
                       dPqP_hot  = 0.01,
                       dPqP_cold = 0.01)
    mode   = isnothing(UA) ? :effectiveness : :UA
    UA_val = isnothing(UA) ? NaN : UA
    T = promote_type(typeof(ε), typeof(UA_val), typeof(dPqP_hot), typeof(dPqP_cold))
    HeatExchanger{T}(name, T(ε), T(UA_val), mode, T(dPqP_hot), T(dPqP_cold),
                     nothing, nothing, nothing, nothing)
end

"""
    _effectiveness_NTU_counterflow(NTU, Cr) -> ε

Counter-flow ε-NTU relation, with the balanced-stream limit ε = NTU/(1+NTU)
as Cr → 1.
"""
function _effectiveness_NTU_counterflow(NTU, Cr)
    if abs(1 - Cr) < 1e-9
        NTU / (1 + NTU)
    else
        e = exp(-NTU * (1 - Cr))
        (1 - e) / (1 - Cr * e)
    end
end

"""Effectiveness used for the current pass (fixed or computed from UA)."""
function _current_effectiveness(hx::HeatExchanger, sh::FluidState, sc::FluidState)
    hx.mode == :effectiveness && return hx.ε

    C_hot  = sh.W * cp(sh.fluid, sh.Tt, sh.Pt)
    C_cold = sc.W * cp(sc.fluid, sc.Tt, sc.Pt)
    C_min  = min(C_hot, C_cold)
    C_max  = max(C_hot, C_cold)

    ε_val = _effectiveness_NTU_counterflow(hx.UA / C_min, C_min / C_max)
    # Store for reporting; in AD context ε_val may be Dual while hx.ε is Float64
    ε_val isa typeof(hx.ε) && (hx.ε = ε_val)
    ε_val
end

function _qmax_enthalpy(sh::FluidState, sc::FluidState)
    fp_h = sh.fluid
    fp_c = sc.fluid

    h_h_in = enthalpy(fp_h, sh.Tt, sh.Pt)
    h_c_in = enthalpy(fp_c, sc.Tt, sc.Pt)
    h_h_at_cold = enthalpy(fp_h, sc.Tt, sh.Pt)
    h_c_at_hot = enthalpy(fp_c, sh.Tt, sc.Pt)

    Q_hot_limit = sh.W * (h_h_in - h_h_at_cold)
    Q_cold_limit = sc.W * (h_c_at_hot - h_c_in)

    # For reverse temperature ordering, both limits are negative and the
    # physically limiting magnitude is the value closer to zero.
    sh.Tt >= sc.Tt ? min(Q_hot_limit, Q_cold_limit) :
                     max(Q_hot_limit, Q_cold_limit)
end

"""
    compute_hx!(hx) -> (hot_outlet, cold_outlet)

Compute both outlet states.  Call this once per iteration after both inlet
ports have been written by their upstream elements.
"""
function compute_hx!(hx::HeatExchanger)
    isnothing(hx.hot_inlet)  && error("HeatExchanger $(hx.name): hot_inlet not set")
    isnothing(hx.cold_inlet) && error("HeatExchanger $(hx.name): cold_inlet not set")

    sh = hx.hot_inlet[]
    sc = hx.cold_inlet[]

    fp_h = sh.fluid
    fp_c = sc.fluid

    Q_max = _qmax_enthalpy(sh, sc)
    Q     = _current_effectiveness(hx, sh, sc) * Q_max   # total heat transferred [W]

    Pt_h_out = sh.Pt * (1 - hx.dPqP_hot)
    Pt_c_out = sc.Pt * (1 - hx.dPqP_cold)

    h_h_out = enthalpy(fp_h, sh.Tt, sh.Pt) - Q / sh.W
    h_c_out = enthalpy(fp_c, sc.Tt, sc.Pt) + Q / sc.W

    Tt_h_out = T_from_h(fp_h, h_h_out, Pt_h_out; T_guess = sh.Tt)
    Tt_c_out = T_from_h(fp_c, h_c_out, Pt_c_out; T_guess = sc.Tt)

    hx.hot_outlet  = Port(update(sh; Pt = Pt_h_out, Tt = Tt_h_out))
    hx.cold_outlet = Port(update(sc; Pt = Pt_c_out, Tt = Tt_c_out))

    (hx.hot_outlet, hx.cold_outlet)
end

# HeatExchanger is triggered through compute_hx!, not the standard compute! path.
# The FlowNetwork handles the two-stream topology specially.
compute!(el::HeatExchanger, inlet::Port) = error(
    "Use compute_hx!(hx) for HeatExchanger elements; call add_hx_pair! in FlowNetwork.")

n_residuals(el::HeatExchanger) = 0
residuals(el::HeatExchanger)   = Float64[]
indep_vars(el::HeatExchanger)  = Float64[]
set_indep_vars!(el::HeatExchanger, x::AbstractVector) = nothing

"""Heat transferred from hot to cold stream [W]."""
function Q_transferred(hx::HeatExchanger)
    isnothing(hx.hot_inlet)  && return 0.0
    if !isnothing(hx.hot_outlet)
        sh = hx.hot_inlet[]
        oh = hx.hot_outlet[]
        return sh.W * (enthalpy(sh) - enthalpy(oh))
    end

    sh = hx.hot_inlet[]
    sc = hx.cold_inlet[]
    _current_effectiveness(hx, sh, sc) * _qmax_enthalpy(sh, sc)
end
