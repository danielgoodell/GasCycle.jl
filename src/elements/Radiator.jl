"""
Radiator element — space heat rejection by thermal radiation.

The working fluid cools along the radiator as it radiates to an effective
sink temperature:

  dQ = σ · ε_r · dA · (T⁴ − T_sink⁴)

Because T falls along the flow path, the panel is integrated in `N_seg`
segments (Heun predictor-corrector per segment, so segment count converges
quickly).  A small total pressure loss dPqP is distributed evenly along the
length.

Operating modes:
  :fixed_area   — area A is given; Tt_out follows from the integration.
                  This is the off-design mode: compressor inlet temperature
                  responds to the operating point.
  :fixed_TtExit — target Tt_out is given; the required area A is computed
                  (design sizing) and stored, so a subsequent :fixed_area
                  run reproduces the design point.

Radiation is two-way: if the fluid enters below T_sink the flux reverses
and the stream warms toward the sink.
"""
mutable struct Radiator{T<:Real} <: AbstractElement
    name::String
    A::T            # total radiating area [m²] (computed in :fixed_TtExit mode)
    emissivity::T   # surface emissivity (0–1)
    T_sink::T       # effective radiation sink temperature [K]
    TtExit::T       # target exit total temperature [K] (:fixed_TtExit mode)
    dPqP::T         # fractional total pressure loss
    mode::Symbol    # :fixed_area or :fixed_TtExit
    N_seg::Int

    inlet::Union{Port, Nothing}
    outlet::Union{Port, Nothing}
end

const _σ_SB = 5.670374419e-8  # Stefan-Boltzmann [W/(m²·K⁴)]

function Radiator(name::String;
                  A          = 0.0,
                  emissivity = 0.85,
                  T_sink     = 200.0,
                  TtExit     = 400.0,
                  dPqP       = 0.01,
                  mode::Symbol = :fixed_area,
                  N_seg::Int   = 50)
    mode in (:fixed_area, :fixed_TtExit) ||
        error("Radiator $name: mode must be :fixed_area or :fixed_TtExit, got :$mode")
    mode == :fixed_area && A <= 0 &&
        error("Radiator $name: :fixed_area mode needs A > 0 (or use mode=:fixed_TtExit to size it)")
    T = promote_type(typeof(A), typeof(emissivity), typeof(T_sink),
                     typeof(TtExit), typeof(dPqP))
    Radiator{T}(name, T(A), T(emissivity), T(T_sink), T(TtExit), T(dPqP),
                mode, N_seg, nothing, nothing)
end

"""Segment-march the radiator at fixed area; returns (Tt_out, Pt_out)."""
function _radiate_fixed_area(el::Radiator, s::FluidState)
    fp    = s.fluid
    σεdA  = _σ_SB * el.emissivity * el.A / el.N_seg
    Ts4   = el.T_sink^4
    dPqPs = el.dPqP / el.N_seg

    T = s.Tt
    P = s.Pt
    h = enthalpy(fp, T, P)
    for _ in 1:el.N_seg
        P_next = P * (1 - dPqPs)
        # Heun: predict with inlet-T flux, correct with mean of T⁴ endpoints
        h_pred = h - σεdA * (T^4 - Ts4) / s.W
        T_pred = T_from_h(fp, h_pred, P_next; T_guess = T)
        h_next = h - σεdA * (0.5 * (T^4 + T_pred^4) - Ts4) / s.W
        T = T_from_h(fp, h_next, P_next; T_guess = T_pred)
        h = h_next
        P = P_next
    end
    (T, P)
end

"""Size the radiator area for a target exit temperature; returns (A, Pt_out)."""
function _radiate_size_area(el::Radiator, s::FluidState)
    fp = s.fluid
    Pt_out = s.Pt * (1 - el.dPqP)
    Tt_out = el.TtExit

    Tt_out > el.T_sink || error(
        "Radiator $(el.name): TtExit=$(el.TtExit) K must exceed T_sink=$(el.T_sink) K " *
        "(radiative rejection cannot cool below the sink)")
    s.Tt > Tt_out || error(
        "Radiator $(el.name): inlet Tt=$(s.Tt) K must exceed TtExit=$(el.TtExit) K")

    # Trapezoid in enthalpy: dA = W dh / (σ ε (T⁴ − Tsink⁴)), with segment
    # pressure ramped linearly for the property evaluations.
    h_in  = enthalpy(fp, s.Tt, s.Pt)
    h_out = enthalpy(fp, Tt_out, Pt_out)
    dh    = (h_out - h_in) / el.N_seg
    σε    = _σ_SB * el.emissivity
    Ts4   = el.T_sink^4

    A = zero(h_in)
    T_prev = s.Tt
    for i in 1:el.N_seg
        P_i = s.Pt + (Pt_out - s.Pt) * (i / el.N_seg)
        T_i = T_from_h(fp, h_in + i * dh, P_i; T_guess = T_prev)
        A  += s.W * (-dh) / (σε * (0.5 * (T_prev^4 + T_i^4) - Ts4))
        T_prev = T_i
    end
    (A, Pt_out)
end

function compute!(el::Radiator, inlet::Port)::Port
    el.inlet = inlet
    s = inlet[]

    if el.mode == :fixed_TtExit
        A_val, Pt_out = _radiate_size_area(el, s)
        # Only store if type matches; in AD context A_val may be Dual while el.A is Float64
        A_val isa typeof(el.A) && (el.A = A_val)
        Tt_out = el.TtExit
    else  # :fixed_area
        Tt_out, Pt_out = _radiate_fixed_area(el, s)
    end

    el.outlet = Port(update(s; Pt = Pt_out, Tt = Tt_out))
    el.outlet
end

n_residuals(el::Radiator)  = 0
residuals(el::Radiator)    = Float64[]
indep_vars(el::Radiator)   = Float64[]
set_indep_vars!(el::Radiator, x::AbstractVector) = nothing

"""Heat rejected to the sink [W] (positive = heat out of the fluid)."""
function Q_rejected(el::Radiator)
    (isnothing(el.inlet) || isnothing(el.outlet)) && return 0.0
    s = el.inlet[]
    s.W * (enthalpy(s) - enthalpy(el.outlet[]))
end
