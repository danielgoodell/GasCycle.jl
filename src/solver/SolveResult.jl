# ── Result type ───────────────────────────────────────────────────────────────

struct SolveResult
    net::FlowNetwork
    status::Symbol
    iterations::Int
    residual_norm::Real  # Real to accept both Float64 and ForwardDiff Duals
end

Base.show(io::IO, r::SolveResult) = begin
    rn = r.residual_norm
    rn_str = rn isa AbstractFloat ? string(round(rn, sigdigits=3)) : string(rn)
    print(io, "SolveResult($(r.status), $(r.iterations) iters, |F|=$rn_str)")
end

function Base.getindex(r::SolveResult, name::String)
    el = findfirst(e -> hasproperty(e, :name) && e.name == name, r.net.elements)
    isnothing(el) && error("Element '$name' not found")
    r.net.elements[el]
end
