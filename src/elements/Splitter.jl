"""
Splitter — divides one inlet stream into N outlet streams.

All outlets share the same total pressure and temperature as the inlet;
only mass flow is divided according to `fracs` (which must sum to 1.0).

Port naming:
  inlet      — single inlet
  outlet     — primary outlet (fracs[1], used by default serial connect!)
  bleed_outlet — secondary outlet (fracs[2], 2-way split only)
  outlet_3, outlet_4, ... — additional outlets for N-way splits
"""
mutable struct Splitter <: AbstractElement
    name::String
    fracs::Vector{Float64}
    inlet::Union{Port,Nothing}
    outlets::Vector{Union{Port,Nothing}}
end

function Splitter(name::String; fracs::Vector{Float64})
    abs(sum(fracs) - 1.0) > 1e-8 &&
        error("Splitter \"$name\": fracs must sum to 1.0, got $(sum(fracs))")
    Splitter(name, fracs, nothing, fill(nothing, length(fracs)))
end

function compute!(sp::Splitter, inlet::Port)
    s = inlet[]
    sp.inlet   = inlet
    sp.outlets = [Port(FluidState(s.Pt, s.Tt, s.W * f, s.fluid)) for f in sp.fracs]
    sp.outlets[1]
end

n_residuals(::Splitter)                        = 0
residuals(::Splitter)                          = Float64[]
indep_vars(::Splitter)                         = Float64[]
set_indep_vars!(::Splitter, ::AbstractVector)  = nothing
