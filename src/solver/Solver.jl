"""
Solver for the FlowNetwork.

Two modes of operation:

1. **Fixed-point iteration** (n_indeps == 0):
   Used for cycles where all design parameters are specified but there are
   circular dependencies (e.g., the recuperator hot inlet depends on the
   turbine outlet which depends on the recuperator cold outlet).
   Iterates one_pass! until outlet temperatures converge.

2. **Newton-Raphson** (n_indeps > 0):
   Used for off-design where shaft speed, map operating points, etc. are
   unknowns with associated residual equations.
   Uses a forward-difference Jacobian so element fields can stay Float64.
"""

using LinearAlgebra: norm

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

# ── Newton-Raphson with forward-difference Jacobian ───────────────────────────

function _fd_newton!(eval!::Function, x::Vector{Float64};
                     tol::Float64 = 1e-6, maxiter::Int = 100)
    n  = length(x)
    Fv = zeros(n)
    Fp = zeros(n)
    J  = zeros(n, n)
    ε  = cbrt(eps(Float64))   # ≈ 6e-6, good for FD with Float64

    eval!(Fv, x)
    for iter in 1:maxiter
        rn = norm(Fv)
        rn < tol && return (iter - 1, rn, true)

        for j in 1:n
            h    = ε * max(abs(x[j]), 1.0)
            xj   = x[j]
            x[j] = xj + h
            eval!(Fp, x)
            @. J[:, j] = (Fp - Fv) / h
            x[j] = xj
        end

        x .+= J \ (-Fv)
        eval!(Fv, x)
    end
    (maxiter, norm(Fv), norm(Fv) < tol)
end

# ── Fixed-point iteration (no unknowns, just iterate to convergence) ──────────

"""Collect all element outlet Tt values for convergence tracking."""
function _outlet_Tt_vec(net::FlowNetwork)
    vals = Float64[]
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

# ── Main solve ────────────────────────────────────────────────────────────────

"""
    solve!(net; tol=1e-6, maxiter=100, verbose=false) -> SolveResult

Solve the flow network.  Design-mode networks with circular dependencies
(recuperators) are resolved by fixed-point iteration.  Off-design networks
with explicit unknowns are solved by Newton-Raphson.
"""
function solve!(net::FlowNetwork;
                tol::Float64  = 1e-6,
                maxiter::Int  = 100,
                verbose::Bool = false)

    # Seed element states (needed before collecting indeps or convergence vec)
    one_pass!(net)

    x0, slices = _collect_indeps(net)
    n_x   = length(x0)
    n_res = sum(n_residuals(el) for el in net.elements; init=0) +
            sum(n_residuals(sh) for sh in net.shafts; init=0)

    n_x == n_res || @warn "Solver: $n_x indeps but $n_res residuals — system may be over/under-determined"

    if n_x == 0
        # Fixed-point iteration: keep running one_pass! until outlet Tt converges.
        # This handles recuperators whose hot inlet depends on the turbine outlet
        # which in turn depends on the recuperator cold outlet.
        prev_Tt = _outlet_Tt_vec(net)
        for iter in 1:maxiter
            one_pass!(net)
            new_Tt = _outlet_Tt_vec(net)
            if isempty(new_Tt)
                return SolveResult(net, :success, iter, 0.0)
            end
            δ = maximum(abs, new_Tt .- prev_Tt)
            verbose && println("  iter $iter: max ΔTt = $(round(δ, sigdigits=3)) K")
            δ < tol && return SolveResult(net, :success, iter, δ)
            prev_Tt = new_Tt
        end
        δ = maximum(abs, _outlet_Tt_vec(net) .- prev_Tt)
        @warn "GasCycle solver (fixed-point) did not converge after $maxiter iterations"
        return SolveResult(net, :failed, maxiter, δ)
    end

    # Newton-Raphson path (off-design)
    function eval!(F, x)
        _scatter_indeps!(net, x, slices)
        one_pass!(net)
        r = _collect_residuals(net)
        copyto!(F, r)
    end

    x = copy(x0)
    iters, rn, converged = _fd_newton!(eval!, x; tol=tol, maxiter=maxiter)

    if verbose
        status = converged ? "converged" : "DID NOT CONVERGE"
        println("Solver (Newton): $status in $iters iterations, |F| = $(round(rn, sigdigits=3))")
    end
    converged || @warn "GasCycle solver (Newton) did not converge after $maxiter iterations (|F|=$rn)"

    _scatter_indeps!(net, x, slices)
    one_pass!(net)

    SolveResult(net, converged ? :success : :failed, iters, rn)
end

# ── Result type ───────────────────────────────────────────────────────────────

struct SolveResult
    net::FlowNetwork
    status::Symbol
    iterations::Int
    residual_norm::Float64
end

Base.show(io::IO, r::SolveResult) =
    print(io, "SolveResult($(r.status), $(r.iterations) iters, |F|=$(round(r.residual_norm, sigdigits=3)))")

function Base.getindex(r::SolveResult, name::String)
    el = findfirst(e -> hasproperty(e, :name) && e.name == name, r.net.elements)
    isnothing(el) && error("Element '$name' not found")
    r.net.elements[el]
end

# ── Cycle-level performance quantities ────────────────────────────────────────

"""
    net_power(result) -> Float64

Net shaft power [W] = total turbine work - total compressor work.
This is the power available to the generator in a closed power cycle.
"""
function net_power(r::SolveResult)
    W_turb = 0.0
    W_comp = 0.0
    for el in r.net.elements
        if el isa Turbine && !isnothing(el.inlet)
            W_turb += specific_work(el) * el.inlet[].W
        elseif el isa Compressor && !isnothing(el.inlet)
            W_comp += specific_work(el) * el.inlet[].W
        end
    end
    W_turb - W_comp
end

"""
    cycle_efficiency(result) -> Float64

Thermal efficiency = net power / total heat input from all HeatSource elements.
"""
function cycle_efficiency(r::SolveResult)
    Q_in = sum(el.Q for el in r.net.elements if el isa HeatSource; init=0.0)
    Q_in ≈ 0.0 && return 0.0
    net_power(r) / Q_in
end
