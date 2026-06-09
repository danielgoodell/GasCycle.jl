"""
Shaft element — tracks the mechanical power balance.

Links turbines (power producers) to compressors (power consumers) on a
common shaft.

In :design mode for a CLOSED POWER CYCLE:
  - The turbine PR is set by pressure closure (not shaft balance).
  - Net shaft power = W_turbine - W_compressor > 0 is the generator output.
  - The shaft has NO solver residual in design mode.

In :off_design mode:
  - Shaft speed N is the independent variable.
  - The power balance (W_turbine - W_compressor = 0) is enforced as a residual.
    This is appropriate when turbine PR comes from a performance map and shaft
    speed is adjusted until the power balance closes.
"""
mutable struct Shaft <: AbstractElement
    name::String
    N::Float64                       # shaft speed [rpm]
    mode::Symbol
    producers::Vector{AbstractElement}   # turbines
    consumers::Vector{AbstractElement}   # compressors
end

function Shaft(name::String; N::Float64 = 10000.0, mode::Symbol = :design)
    Shaft(name, N, mode, AbstractElement[], AbstractElement[])
end

"""Called by FlowNetwork after connect! to wire turbines and compressors."""
function link!(shaft::Shaft, producers, consumers)
    append!(shaft.producers, producers)
    append!(shaft.consumers, consumers)
    _broadcast_speed!(shaft)
end

function _broadcast_speed!(shaft::Shaft)
    for el in Iterators.flatten((shaft.producers, shaft.consumers))
        el.N_shaft = shaft.N
    end
end

# Shaft has no flow path — compute! is a no-op
compute!(el::Shaft, inlet::Port) = inlet

# In :design mode, shaft balance is an output (not a solver constraint).
# In :off_design mode, power balance is the residual; N is the indep var.
n_residuals(el::Shaft) = el.mode == :off_design ? 1 : 0

function residuals(el::Shaft)
    el.mode == :off_design || return Float64[]
    _broadcast_speed!(el)
    W_prod  = sum(specific_work(t) * t.inlet[].W for t in el.producers; init=0.0)
    W_cons  = sum(specific_work(c) * c.inlet[].W for c in el.consumers; init=0.0)
    W_scale = max(abs(W_prod), abs(W_cons), 1.0)
    [(W_prod - W_cons) / W_scale]
end

indep_vars(el::Shaft) = el.mode == :off_design ? [el.N] : Float64[]
function set_indep_vars!(el::Shaft, x::AbstractVector)
    if el.mode == :off_design
        el.N = x[1]
        _broadcast_speed!(el)
    end
end

"""
    power_balance(shaft) -> Float64

Net shaft power [W] = turbine work - compressor work.
Positive value goes to the generator in a closed power cycle.
"""
function power_balance(shaft::Shaft)
    W_prod = sum(specific_work(t) * t.inlet[].W for t in shaft.producers; init=0.0)
    W_cons = sum(specific_work(c) * c.inlet[].W for c in shaft.consumers; init=0.0)
    W_prod - W_cons
end
