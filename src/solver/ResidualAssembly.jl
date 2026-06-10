# ── Helper: gather/scatter independent variables ──────────────────────────────

function _collect_indeps(net::FlowNetwork)
    x      = Float64[]
    slices = UnitRange{Int}[]
    for el in net.elements
        v = indep_vars(el)
        push!(slices, length(x)+1 : length(x)+length(v))
        append!(x, v)
    end
    for sh in net.shafts
        v = indep_vars(sh)
        push!(slices, length(x)+1 : length(x)+length(v))
        append!(x, v)
    end
    (x, slices)
end

function _scatter_indeps!(net::FlowNetwork, x::AbstractVector, slices)
    k = 1
    for el in net.elements
        set_indep_vars!(el, view(x, slices[k]))
        k += 1
    end
    for sh in net.shafts
        set_indep_vars!(sh, view(x, slices[k]))
        k += 1
    end
end

function _collect_residuals(net::FlowNetwork)
    F = Float64[]
    for el in net.elements
        append!(F, residuals(el))
    end
    for sh in net.shafts
        append!(F, residuals(sh))
    end
    F
end

"""Collect all element outlet Tt values for convergence tracking."""
function _outlet_Tt_vec(net::FlowNetwork)
    vals = Real[]
    for el in net.elements
        if el isa HeatExchanger
            isnothing(el.cold_outlet) || push!(vals, el.cold_outlet[].Tt)
            isnothing(el.hot_outlet)  || push!(vals, el.hot_outlet[].Tt)
        elseif el isa Splitter
            for out in el.outlets
                isnothing(out) || push!(vals, out[].Tt)
            end
        elseif hasproperty(el, :outlet) && !isnothing(el.outlet)
            push!(vals, el.outlet[].Tt)
        end
    end
    vals
end

"""
Build the initial guess for the back-edge Newton.
Returns [Tt₁, Pt₁, Tt₂, Pt₂, …] from stored element outlet states,
falling back to seed_state when no outlet is stored yet.
The returned vector has a concrete element type (Float64 or Dual) inferred
from the outlet states — NonlinearSolve requires a concrete-typed u0.
"""
function _back_edge_z0(net::FlowNetwork, back_edges)
    # Collect pairs [Tt, Pt] per back-edge; vcat infers the concrete element type.
    pairs = map(back_edges) do edge
        src_out = _get_outlet(net.elements[edge.src], edge.src_port)
        if isnothing(src_out)
            [net.seed_state.Tt, net.seed_state.Pt]
        else
            [src_out[].Tt, src_out[].Pt]
        end
    end
    isempty(pairs) ? Float64[] : reduce(vcat, pairs)
end
