"""
Counter-flow heat exchanger / recuperator.

Models a two-stream heat exchanger using the ε-NTU method.  The two streams
(hot and cold) must be declared as a pair and linked in the FlowNetwork via
`add_hx_pair!`.

  Q = ε * Q_max
  Q_max is computed from the limiting endpoint enthalpy change.

Both streams conserve mass flow; only Tt changes (pressure loss dPqP applied
to each stream independently).

Usage in FlowNetwork:
  hx = HeatExchanger("Recup"; ε=0.95)
  add_hx_pair!(net, hx, hot_inlet_port, cold_inlet_port)
"""
mutable struct HeatExchanger{T<:Real} <: AbstractElement
    name::String
    ε::T           # effectiveness (0–1)
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
                       dPqP_hot  = 0.01,
                       dPqP_cold = 0.01)
    T = promote_type(typeof(ε), typeof(dPqP_hot), typeof(dPqP_cold))
    HeatExchanger{T}(name, T(ε), T(dPqP_hot), T(dPqP_cold),
                     nothing, nothing, nothing, nothing)
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
    Q     = hx.ε * Q_max           # total heat transferred [W]

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
    hx.ε * _qmax_enthalpy(sh, sc)
end
