# ── Helper: gather/scatter independent variables ──────────────────────────────

# Promote a heterogeneous Real vector to its common numeric type (Float64 when
# empty). Element unknowns/residuals may mix Float64 and Dual values.
_promote_vec(v) = isempty(v) ? Float64[] :
    convert(Vector{mapreduce(typeof, promote_type, v)}, v)

function _collect_indeps(net::FlowNetwork)
    x      = Real[]
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
    for l in net.heat_links
        v = indep_vars(l)
        push!(slices, length(x)+1 : length(x)+length(v))
        append!(x, v)
    end
    (_promote_vec(x), slices)
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
    for l in net.heat_links
        set_indep_vars!(l, view(x, slices[k]))
        k += 1
    end
end

function _collect_residuals(net::FlowNetwork)
    F = Real[]
    for el in net.elements
        append!(F, residuals(el))
    end
    for sh in net.shafts
        append!(F, residuals(sh))
    end
    for l in net.heat_links
        append!(F, residuals(l))
    end
    _promote_vec(F)
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
from the outlet states — the Newton solvers and ForwardDiff need a concrete-typed u0.
"""
function _back_edge_z0(net::FlowNetwork, back_edges)
    # Collect pairs [Tt, Pt] per back-edge; vcat infers the concrete element type.
    pairs = map(back_edges) do edge
        src_out = _get_outlet(net.elements[edge.src], edge.src_port)
        if isnothing(src_out)
            g = _guess_state(net, edge)
            [g.Tt, g.Pt]
        else
            [src_out[].Tt, src_out[].Pt]
        end
    end
    isempty(pairs) ? Float64[] : reduce(vcat, pairs)
end
