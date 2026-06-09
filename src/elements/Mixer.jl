"""
Mixer — combines N inlet streams into one outlet by enthalpy-mass balance.

  W_out  = Σ Wᵢ
  h_out  = Σ(Wᵢ hᵢ) / W_out          (enthalpy balance, not Cp·ΔT)
  Pt_out = min(Ptᵢ)                   (conservative; real mixing has irreversibility)
  Tt_out = T_from_h(fluid, h_out, Pt_out)

Port naming:
  inlet        — primary inlet (stream 1, used by default serial connect!)
  bleed_inlet  — secondary inlet (stream 2, 2-way mix only)
  inlet_3, inlet_4, ... — additional inlets for N-way mixes
  outlet       — single outlet
"""
mutable struct Mixer <: AbstractElement
    name::String
    n_inlets::Int
    inlets::Vector{Union{Port,Nothing}}
    outlet::Union{Port,Nothing}
end

function Mixer(name::String; n_inlets::Int = 2)
    n_inlets >= 2 || error("Mixer \"$name\": n_inlets must be ≥ 2")
    Mixer(name, n_inlets, Vector{Union{Port,Nothing}}(fill(nothing, n_inlets)), nothing)
end

function compute!(mx::Mixer)
    all(!isnothing, mx.inlets) ||
        error("Mixer \"$(mx.name)\": not all inlets have been set")
    states  = [p[] for p in mx.inlets]
    W_tot   = sum(s.W for s in states)
    fluid   = states[1].fluid
    h_mix   = sum(s.W * enthalpy(s.fluid, s.Tt, s.Pt) for s in states) / W_tot
    Pt_mix  = minimum(s.Pt for s in states)
    Tt_mix  = T_from_h(fluid, h_mix, Pt_mix)
    mx.outlet = Port(FluidState(Pt_mix, Tt_mix, W_tot, fluid))
end

n_residuals(::Mixer)                       = 0
residuals(::Mixer)                         = Float64[]
indep_vars(::Mixer)                        = Float64[]
set_indep_vars!(::Mixer, ::AbstractVector) = nothing
