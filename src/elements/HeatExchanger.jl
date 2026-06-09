"""
Counter-flow heat exchanger / recuperator.

Models a two-stream heat exchanger using the ε-NTU method.  The two streams
(hot and cold) must be declared as a pair and linked in the FlowNetwork via
`add_hx_pair!`.

  Q = ε * Q_max
  Q_max = min(Ċ_hot, Ċ_cold) * (Tt_hot_in - Tt_cold_in)
  Ċ = ṁ * cp    [W/K]

Both streams conserve mass flow; only Tt changes (pressure loss dPqP applied
to each stream independently).

Usage in FlowNetwork:
  hx = HeatExchanger("Recup"; ε=0.95)
  add_hx_pair!(net, hx, hot_inlet_port, cold_inlet_port)
"""
mutable struct HeatExchanger <: AbstractElement
    name::String
    ε::Float64        # effectiveness (0–1)
    dPqP_hot::Float64
    dPqP_cold::Float64

    # Set by add_hx_pair! before first compute!
    hot_inlet::Union{Port, Nothing}
    cold_inlet::Union{Port, Nothing}
    hot_outlet::Union{Port, Nothing}
    cold_outlet::Union{Port, Nothing}
end

function HeatExchanger(name::String;
                       ε::Float64        = 0.90,
                       dPqP_hot::Float64  = 0.01,
                       dPqP_cold::Float64 = 0.01)
    HeatExchanger(name, ε, dPqP_hot, dPqP_cold,
                  nothing, nothing, nothing, nothing)
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

    Cp_h = cp(fp_h, sh.Tt, sh.Pt)
    Cp_c = cp(fp_c, sc.Tt, sc.Pt)

    Cdot_h = sh.W * Cp_h   # W/K
    Cdot_c = sc.W * Cp_c

    Q_max = min(Cdot_h, Cdot_c) * (sh.Tt - sc.Tt)
    Q     = hx.ε * Q_max           # total heat transferred [W]

    Tt_h_out = sh.Tt - Q / Cdot_h
    Tt_c_out = sc.Tt + Q / Cdot_c

    Pt_h_out = sh.Pt * (1.0 - hx.dPqP_hot)
    Pt_c_out = sc.Pt * (1.0 - hx.dPqP_cold)

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
set_indep_vars!(el::HeatExchanger, x) = nothing

"""Heat transferred from hot to cold stream [W]."""
function Q_transferred(hx::HeatExchanger)
    isnothing(hx.hot_inlet)  && return 0.0
    sh = hx.hot_inlet[]
    sc = hx.cold_inlet[]
    Cp_h  = cp(sh.fluid, sh.Tt, sh.Pt)
    Cp_c  = cp(sc.fluid, sc.Tt, sc.Pt)
    Q_max = min(sh.W * Cp_h, sc.W * Cp_c) * (sh.Tt - sc.Tt)
    hx.ε * Q_max
end
