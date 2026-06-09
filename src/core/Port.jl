"""
A mutable connector that holds a FluidState reference.

Components write to their outlet Port; the next component reads from it
as its inlet Port. Since Ports are shared by reference, writing the outlet
of one element automatically updates the inlet of the downstream element —
no explicit copy step is needed after `connect!`.

The type parameter T mirrors FluidState{T}: Port{Float64} for normal use,
Port{Dual{...}} when automatic differentiation is active.
"""
mutable struct Port{T<:Real}
    state::FluidState{T}
    name::String
end

Port(state::FluidState{T}) where {T<:Real} = Port{T}(state, "")

"""Read the current FluidState from a port."""
Base.getindex(p::Port) = p.state

"""Write a new FluidState to a port (T must match port's parameter)."""
Base.setindex!(p::Port{T}, s::FluidState{T}) where T = (p.state = s; p)

# Forwarding accessors so element code can write `port.Pt` instead of `port[].Pt`
Base.getproperty(p::Port, sym::Symbol) =
    sym in (:state, :name) ? getfield(p, sym) : getproperty(p.state, sym)
