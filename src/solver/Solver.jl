"""
Solver for the FlowNetwork.

Two modes of operation:

1. **Back-edge Newton** (n_indeps == 0):
   Used for cycles where all design parameters are specified but there are
   circular dependencies (e.g., the recuperator hot inlet depends on the
   turbine outlet which depends on the recuperator cold outlet).
   NonlinearSolve.jl Newton treats [Tt, Pt] per back-edge as unknowns;
   the residual is the fixed-point consistency condition.
   Using AutoForwardDiff() for the Newton Jacobian lets outer
   ForwardDiff.jl differentiation thread through correctly for design
   sensitivity — outer Dual numbers propagate as values while the inner
   Newton Jacobian is formed by a separate, tagged ForwardDiff pass.

2. **Newton-Raphson with finite differences** (n_indeps > 0):
   Used for off-design where shaft speed, map operating points, etc. are
   unknowns with associated residual equations.  Any back-edge (Tt, Pt)
   states are appended to the unknown vector and solved simultaneously, so
   closed loops and off-design unknowns share one Newton iteration.
   Uses AutoFiniteDiff() because element struct fields are Float64.
"""

using LinearAlgebra: norm
using NonlinearSolve

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

# ── Convergence tracking ──────────────────────────────────────────────────────

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

# ── Back-edge helpers ─────────────────────────────────────────────────────────

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

# ── Main solve ────────────────────────────────────────────────────────────────

"""
    solve!(net; tol=1e-6, maxiter=100, verbose=false) -> SolveResult

Solve the flow network.

Design-mode networks with circular dependencies (recuperators) are resolved by
a Newton solve on the back-edge [Tt, Pt] states with AutoForwardDiff() Jacobian.
Off-design networks with explicit unknowns use Newton with AutoFiniteDiff().
"""
function solve!(net::FlowNetwork;
                tol::Float64  = 1e-6,
                maxiter::Int  = 100,
                verbose::Bool = false)

    # Initial pass to populate element outlet states (used for z0 below).
    one_pass!(net)

    x0, slices = _collect_indeps(net)
    n_x   = length(x0)
    n_res = sum(n_residuals(el) for el in net.elements; init=0) +
            sum(n_residuals(sh) for sh in net.shafts; init=0)

    n_x == n_res || @warn "Solver: $n_x indeps but $n_res residuals — system may be over/under-determined"

    back_edges = filter(e -> e.back_edge, net.edges)
    n_be = length(back_edges)

    if n_x == 0
        # ── Design mode ───────────────────────────────────────────────────────
        if n_be == 0
            # No circular dependencies — the single pass above is the solution.
            return SolveResult(net, :success, 1, 0.0)
        end

        # Back-edge Newton: z = [Tt₁, Pt₁, Tt₂, Pt₂, …].
        # Residual: each source element's computed outlet must equal the seeded z.
        # Normalise by z so residuals are dimensionless (tol is fractional).
        #
        # OUT-OF-PLACE form (f(z,p) returning F) avoids the nested-Dual
        # pre-allocation conflict when the outer ForwardDiff wraps solve!:
        # in-place f! pre-allocates F at u0's element type, but AutoForwardDiff
        # internally creates Dual{inner, Dual{outer}, N} values that can't be
        # stored back into the outer-Dual-typed buffer.
        function be_residual(z, _p)
            one_pass!(net, z)
            F = similar(z, 2 * length(back_edges))
            for (i, edge) in enumerate(back_edges)
                src_out = _get_outlet(net.elements[edge.src], edge.src_port)
                if isnothing(src_out)
                    F[2i-1] = zero(eltype(z))
                    F[2i]   = zero(eltype(z))
                else
                    F[2i-1] = (src_out[].Tt - z[2i-1]) / z[2i-1]
                    F[2i]   = (src_out[].Pt - z[2i])   / z[2i]
                end
            end
            F
        end

        z0   = _back_edge_z0(net, back_edges)
        prob = NonlinearProblem(be_residual, z0, nothing)
        nls  = NonlinearSolve.solve(prob,
                   NewtonRaphson(autodiff = AutoForwardDiff());
                   abstol   = tol,
                   reltol   = tol,
                   maxiters = maxiter)

        rn        = norm(nls.resid)
        converged = rn < tol
        iters     = nls.stats.nsteps

        # Final pass with converged back-edge states to restore all element fields.
        one_pass!(net, nls.u)

        if verbose
            status = converged ? "converged" : "DID NOT CONVERGE"
            println("Solver (back-edge Newton): $status in $iters iters, |F| = $(round(Float64(rn), sigdigits=3))")
        end
        converged || @warn "GasCycle back-edge Newton did not converge after $maxiter iterations (|F|=$rn)"

        return SolveResult(net, converged ? :success : :failed, iters, rn)
    end

    # ── Off-design mode ───────────────────────────────────────────────────────
    # Unknown vector u = [element/shaft indeps; back-edge (Tt, Pt) pairs].
    # The back-edge states are folded into the same Newton solve so each
    # residual evaluation is a deterministic function of u.  (Seeding back-edges
    # from stored outlets inside one_pass! — the old behaviour — leaves stale
    # state between iterations and corrupts finite-difference Jacobians.)
    #
    # The solve runs in NORMALIZED unknowns v = u ./ uref: raw unknowns mix
    # corrected flows (~0.5), shaft speeds (~10⁴), and back-edge Tt/Pt
    # (~10³/10⁵), which makes the raw Jacobian condition number ~10⁶ and
    # stalls the globalization just above tolerance.  Residuals are already
    # dimensionless, so scaling the columns brings cond(J) to O(10).
    u0   = vcat(x0, _back_edge_z0(net, back_edges))
    uref = [abs(u) > 0 ? abs(u) : 1.0 for u in u0]

    function od_residual!(F, v, _p)
        u = v .* uref
        _scatter_indeps!(net, view(u, 1:n_x), slices)
        z = view(u, n_x+1:length(u))
        one_pass!(net, n_be == 0 ? nothing : z)
        r = _collect_residuals(net)
        copyto!(F, 1, r, 1, length(r))
        for (i, edge) in enumerate(back_edges)
            src_out = _get_outlet(net.elements[edge.src], edge.src_port)
            if isnothing(src_out)
                F[n_res + 2i - 1] = 0.0
                F[n_res + 2i]     = 0.0
            else
                F[n_res + 2i - 1] = (src_out[].Tt - z[2i-1]) / z[2i-1]
                F[n_res + 2i]     = (src_out[].Pt - z[2i])   / z[2i]
            end
        end
        nothing
    end

    # TrustRegion rather than plain Newton for robustness far from the
    # solution (map slope kinks at grid lines make undamped steps unreliable).
    prob = NonlinearProblem(od_residual!, u0 ./ uref, nothing)
    nls  = NonlinearSolve.solve(prob,
               TrustRegion(autodiff = AutoFiniteDiff());
               abstol   = tol,
               reltol   = tol,
               maxiters = maxiter)

    rn        = norm(nls.resid)
    converged = rn < tol
    iters     = nls.stats.nsteps

    u_final = nls.u .* uref
    _scatter_indeps!(net, view(u_final, 1:n_x), slices)
    one_pass!(net, n_be == 0 ? nothing : u_final[n_x+1:end])

    if verbose
        status = converged ? "converged" : "DID NOT CONVERGE"
        println("Solver (off-design Newton): $status in $iters iters, |F| = $(round(Float64(rn), sigdigits=3))")
    end
    converged || @warn "GasCycle off-design Newton did not converge after $maxiter iterations (|F|=$rn)"

    SolveResult(net, converged ? :success : :failed, iters, rn)
end

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

# ── Cycle-level performance quantities ────────────────────────────────────────

"""
    net_power(result) -> Real

Net shaft power [W] = total turbine work - total compressor work.
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
    cycle_efficiency(result) -> Real

Thermal efficiency = net power / total heat input from all HeatSource elements.
"""
function cycle_efficiency(r::SolveResult)
    Q_in = sum(r.net.elements; init=0.0) do el
        if el isa HeatSource
            if !isnothing(el.inlet) && !isnothing(el.outlet)
                (enthalpy(el.outlet[]) - enthalpy(el.inlet[])) * el.inlet[].W
            else
                el.Q
            end
        else
            0.0
        end
    end
    Q_in ≈ 0.0 && return 0.0
    net_power(r) / Q_in
end
