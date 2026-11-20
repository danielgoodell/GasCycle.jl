# ── Built-in nonlinear solvers ───────────────────────────────────────────────
#
# Small dense solvers for the cycle-level Newton problem (a few to a few dozen
# unknowns). The cost of a solve is dominated by residual evaluations (a full
# `one_pass!` through the cycle), so these favour few residual calls and simple,
# predictable behaviour over generality.
#
# Interface: a backend is a subtype of `NonlinearBackend` with a method
#
#     nlsolve(alg, f, u0; tol, maxiters, lower, upper) -> NLResult
#
# where `f(u)` is out-of-place and returns the residual vector, and
# `lower`/`upper` are optional box bounds (`nothing` = unbounded). A package
# extension can add a backend (e.g. sparse or NonlinearSolve.jl-based) by
# defining a new subtype and its `nlsolve` method; `solve!(net; alg=...)`
# passes it through.
#
# AD through the solve: the model may carry ForwardDiff Duals (an outer
# gradient w.r.t. design parameters), so `u` and `f(u)` may be Dual-valued.
# Jacobians are always taken on the PRIMAL values (a plain Float64 matrix):
# the unknowns are often stored in concretely typed element fields, so the
# inner Jacobian must never push its own Duals into the model. Newton steps
# u ← u − J⁻¹ F(u) applied to Dual-valued u then converge the derivative
# parts as a linear fixed-point iteration whose fixed point is the
# implicit-function-theorem result F_u·u̇ + F_p·ṗ = 0 (F's Dual parts are exact;
# J only preconditions). After the primal converges, `_refine_partials` keeps
# iterating until the derivative parts settle — one step for an exact J, a few
# for a finite-difference J. Nested Duals (Hessians) converge the same way.

import ForwardDiff
using LinearAlgebra: norm, dot, lu, issuccess

abstract type NonlinearBackend end

"""
    NewtonBacktracking(; jacobian = :forwarddiff, max_backtracks = 10)

Newton's method with a backtracking safeguard: the full step is taken unless
it fails to reduce `‖F‖` (or produces non-finite residuals), in which case it
is halved up to `max_backtracks` times. Default for design-mode solves.
`jacobian = :forwarddiff` requires that the unknowns are not stored in
concretely typed fields (true for back-edge states); use `:finitediff`
otherwise.
"""
Base.@kwdef struct NewtonBacktracking <: NonlinearBackend
    jacobian::Symbol    = :forwarddiff
    max_backtracks::Int = 10
end

"""
    DoglegTrustRegion(; jacobian = :finitediff, initial_radius = 1.0, max_radius = 1e3)

Powell dogleg trust-region method. More robust than plain Newton far from the
solution (e.g. map-slope kinks at grid lines). `initial_radius` is relative to
`max(‖u0‖, 1)`. Default for off-design solves and solves with free variables /
targets, which run in normalized unknowns.
"""
Base.@kwdef struct DoglegTrustRegion <: NonlinearBackend
    jacobian::Symbol        = :finitediff
    initial_radius::Float64 = 1.0
    max_radius::Float64     = 1e3
end

"""Result of `nlsolve`: solution, residual at the solution, iteration count."""
struct NLResult{V<:AbstractVector, R<:AbstractVector}
    u::V
    resid::R
    nsteps::Int
    converged::Bool
end

# ── Dual helpers ─────────────────────────────────────────────────────────────

# _primal (strip all Dual levels) is defined in network/Targets.jl.
_primal_vec(v::AbstractVector) = Float64[_primal(x) for x in v]

_is_dual(u) = eltype(u) <: ForwardDiff.Dual

# Total magnitude of all derivative parts of x (recursing through nesting).
_absall(x::Real) = abs(Float64(x))
_absall(x::ForwardDiff.Dual) = _absall(ForwardDiff.value(x)) + sum(_absall, ForwardDiff.partials(x); init = 0.0)
_pmag(x::Real) = 0.0
_pmag(x::ForwardDiff.Dual) = _pmag(ForwardDiff.value(x)) + sum(_absall, ForwardDiff.partials(x); init = 0.0)

_finite(F) = all(x -> isfinite(_primal(x)), F)
_nval(x) = norm(_primal_vec(x))                     # primal 2-norm

# Promote u0 to the common numeric type of u0 and F(u0): parameters may carry
# Duals even when the initial guess is Float64.
function _promote_unknowns(u0, F0)
    T = promote_type(eltype(u0), eltype(F0))
    convert(Vector{T}, u0)
end

_project(u, lower, upper) =
    isnothing(lower) ? u : [clamp(u[i], lower[i], upper[i]) for i in eachindex(u)]

# ── Jacobians (always primal Float64 matrices) ──────────────────────────────

function _jacobian(method::Symbol, f, u::AbstractVector, Fu::AbstractVector)
    if method === :forwarddiff
        return _primal.(ForwardDiff.jacobian(f, u))
    elseif method === :finitediff
        # Forward differences on the primal values, one residual call per column.
        F0 = _primal_vec(Fu)
        J  = Matrix{Float64}(undef, length(Fu), length(u))
        up = copy(u)
        for j in eachindex(u)
            h = sqrt(eps(Float64)) * max(abs(_primal(u[j])), 1.0)
            up[j] = u[j] + h
            J[:, j] = (_primal_vec(f(up)) .- F0) ./ h
            up[j] = u[j]
        end
        return J
    else
        throw(ArgumentError("unknown jacobian method :$method (use :forwarddiff or :finitediff)"))
    end
end

_factor(J) = (fac = lu(J; check = false); issuccess(fac) ? fac : nothing)

# Newton step −J⁻¹F for a (possibly Dual-valued) F; `nothing` if J is singular.
function _newton_step(fac, F)
    isnothing(fac) && return nothing
    p = -(fac \ F)
    _finite(p) ? p : nothing
end

# After primal convergence, iterate Newton steps with a fixed primal Jacobian
# factorization until the Dual (derivative) parts of u stop changing.
# A no-op for Float64 problems.
function _refine_partials(f, u, F, fac, lower, upper; maxit = 50, rtol = 1e-12)
    _is_dual(u) || return (u, F)
    isnothing(fac) && return (u, F)
    for _ in 1:maxit
        p = _newton_step(fac, F)
        isnothing(p) && break
        u = _project(u .+ p, lower, upper)
        F = f(u)
        scale = max(1.0, maximum(_pmag, u))
        maximum(_pmag, p) <= rtol * scale && break
    end
    (u, F)
end

# ── Newton with backtracking ─────────────────────────────────────────────────

function nlsolve(alg::NewtonBacktracking, f, u0::AbstractVector;
                 tol::Real, maxiters::Int, lower = nothing, upper = nothing)
    F = f(u0)
    u = _promote_unknowns(u0, F)
    _finite(F) || error("nlsolve: residual is not finite at the initial guess")
    nF = _nval(F)
    fac = nothing
    iters = 0
    while nF >= tol && iters < maxiters
        iters += 1
        fac = _factor(_jacobian(alg.jacobian, f, u, F))
        p = _newton_step(fac, F)
        isnothing(p) && break                       # singular Jacobian
        λ = 1.0
        accepted = false
        for _ in 0:alg.max_backtracks
            u_try = _project(u .+ λ .* p, lower, upper)
            F_try = f(u_try)
            if _finite(F_try) && _nval(F_try) < nF
                u, F, nF = u_try, F_try, _nval(F_try)
                accepted = true
                break
            end
            λ /= 2
        end
        accepted || break                           # no descent along Newton direction
    end
    converged = nF < tol
    if converged && _is_dual(u)
        isnothing(fac) && (fac = _factor(_jacobian(alg.jacobian, f, u, F)))
        u, F = _refine_partials(f, u, F, fac, lower, upper)
    end
    NLResult(u, F, iters, converged)
end

# ── Powell dogleg trust region ───────────────────────────────────────────────

# Dogleg step for the primal problem (J, F primal Float64).
function _dogleg(fac, J, F, Δ)
    pN = _newton_step(fac, F)
    if !isnothing(pN) && norm(pN) <= Δ
        return pN                                   # Newton step fits
    end
    g  = J' * F                                     # gradient of ½‖F‖²
    ng = norm(g)
    ng == 0 && return zero(g)                       # stationary: no progress possible
    Jg = J * g
    pC = -(dot(g, g) / dot(Jg, Jg)) .* g            # Cauchy point
    if isnothing(pN) || norm(pC) >= Δ
        return -(Δ / ng) .* g                       # steepest descent to the boundary
    end
    # Dogleg: pC + τ (pN − pC) with ‖·‖ = Δ, τ ∈ [0, 1]
    d = pN .- pC
    a = dot(d, d); b = 2 * dot(pC, d); c = dot(pC, pC) - Δ^2
    τ = (-b + sqrt(b^2 - 4a * c)) / (2a)
    pC .+ τ .* d
end

function nlsolve(alg::DoglegTrustRegion, f, u0::AbstractVector;
                 tol::Real, maxiters::Int, lower = nothing, upper = nothing)
    F = f(u0)
    u = _promote_unknowns(u0, F)
    _finite(F) || error("nlsolve: residual is not finite at the initial guess")
    nF   = _nval(F)
    Δ    = alg.initial_radius * max(_nval(u), 1.0)
    Δmax = alg.max_radius * max(_nval(u), 1.0)
    J    = nF >= tol ? _jacobian(alg.jacobian, f, u, F) : nothing
    fac  = isnothing(J) ? nothing : _factor(J)
    iters = 0
    while nF >= tol && iters < maxiters
        iters += 1
        # The step is computed on primal values; applied to Dual-valued u it
        # only perturbs derivative parts, which _refine_partials settles.
        Fp = _primal_vec(F)
        p  = _dogleg(fac, J, Fp, Δ)
        u_try = _project(u .+ p, lower, upper)
        p  = _primal_vec(u_try) .- _primal_vec(u)    # actual (projected) step
        np = norm(p)
        np == 0 && break
        pred = nF^2 - norm(Fp .+ J * p)^2           # predicted reduction of ‖F‖²
        F_try = f(u_try)
        ared = _finite(F_try) ? nF^2 - _nval(F_try)^2 : -Inf
        ρ = pred > 0 ? ared / pred : -Inf
        if ρ < 0.25
            Δ = 0.25 * np
        elseif ρ > 0.75 && np >= 0.99 * Δ
            Δ = min(2Δ, Δmax)
        end
        if ρ > 1e-4                                 # accept
            u, F, nF = u_try, F_try, _nval(F_try)
            if nF >= tol
                J   = _jacobian(alg.jacobian, f, u, F)
                fac = _factor(J)
            end
        end
        Δ < 1e-12 * max(_nval(u), 1.0) && break     # radius collapsed: stalled
    end
    converged = nF < tol
    if converged && _is_dual(u)
        isnothing(fac) && (fac = _factor(_jacobian(alg.jacobian, f, u, F)))
        u, F = _refine_partials(f, u, F, fac, lower, upper)
    end
    NLResult(u, F, iters, converged)
end
