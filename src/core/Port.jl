"""
A mutable connector that holds a FluidState reference.

Components write to their outlet Port; the next component reads from it
as its inlet Port. Since Ports are shared by reference, writing the outlet
of one element automatically updates the inlet of the downstream element —
no explicit copy step is needed after `connect!`.
"""
mutable struct Port
    state::FluidState
    name::String
end

Port(state::FluidState) = Port(state, "")

"""Read the current FluidState from a port."""
Base.getindex(p::Port) = p.state

"""Write a new FluidState to a port."""
Base.setindex!(p::Port, s::FluidState) = (p.state = s; p)

# Forwarding accessors so element code can write `port.Pt` instead of `port[].Pt`
Base.getproperty(p::Port, sym::Symbol) =
    sym in (:state, :name) ? getfield(p, sym) : getproperty(p.state, sym)
