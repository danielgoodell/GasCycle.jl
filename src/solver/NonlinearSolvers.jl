# ── Built-in nonlinear solvers ───────────────────────────────────────────────
#
# Small dense solvers for the cycle-level Newton problem (a few to a few dozen
# unknowns). The cost of a solve is dominated by residual evaluations (a full
# `one_pass!` through the cycle), so these favour few residual calls and simple,
# predictable behaviour over generality.
#
# Interface: a backend is a subtype of `NonlinearBackend` with a method
#
#     nlsolve(alg, f, u0; tol, maxiters) -> NLResult
#
# where `f(u)` is out-of-place and returns the residual vector. A package
# extension can add a backend (e.g. sparse or NonlinearSolve.jl-based) by
# defining a new subtype and its `nlsolve` method; `solve!(net; alg=...)`
# passes it through.
#
# AD: `u0` may carry ForwardDiff Duals (an outer gradient through the solve).
# Iterating Newton on Duals converges the derivative part along with the value;
# a final Newton step from the converged point makes the derivative part equal
# to the implicit-function-theorem result. Jacobians from ForwardDiff nest
# correctly inside an outer ForwardDiff pass (different tags).

import ForwardDiff
using LinearAlgebra: norm, dot, lu, issuccess, ldiv!

abstract type NonlinearBackend end

"""
    NewtonBacktracking(; jacobian = :forwarddiff, max_backtracks = 10)

Newton's method with a backtracking safeguard: the full step is taken unless
it fails to reduce `‖F‖` (or produces non-finite residuals), in which case it
is halved up to `max_backtracks` times. Default for design-mode solves.
"""
Base.@kwdef struct NewtonBacktracking <: NonlinearBackend
    jacobian::Symbol    = :forwarddiff
    max_backtracks::Int = 10
end

"""
    DoglegTrustRegion(; jacobian = :finitediff, initial_radius = 1.0, max_radius = 1e3)

Powell dogleg trust-region method. More robust than plain Newton far from the
solution (e.g. map-slope kinks at grid lines). `initial_radius` is relative to
`max(‖u0‖, 1)`. Default for off-design solves, which run in normalized unknowns.
"""
Base.@kwdef struct DoglegTrustRegion <: NonlinearBackend
    jacobian::Symbol       = :finitediff
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

# ── Jacobians ────────────────────────────────────────────────────────────────

function _jacobian(method::Symbol, f, u::AbstractVector, Fu::AbstractVector)
    if method === :forwarddiff
        return ForwardDiff.jacobian(f, u)
    elseif method === :finitediff
        # Forward differences, one residual call per column.
        J = similar(Fu, length(Fu), length(u))
        up = copy(u)
        for j in eachindex(u)
            h = sqrt(eps(Float64)) * max(abs(ForwardDiff.value(u[j])), 1.0)
            up[j] = u[j] + h
            J[:, j] = (f(up) .- Fu) ./ h
            up[j] = u[j]
        end
        return J
    else
        throw(ArgumentError("unknown jacobian method :$method (use :forwarddiff or :finitediff)"))
    end
end

# Solve J p = -F; `nothing` if J is singular.
function _newton_step(J, F)
    fac = lu(J; check = false)
    issuccess(fac) || return nothing
    p = -(fac \ F)
    all(isfinite, p) ? p : nothing
end

_finite(F) = all(isfinite, F)
_nval(x) = ForwardDiff.value(norm(x))     # primal norm, also for Dual vectors
_is_dual(u) = eltype(u) <: ForwardDiff.Dual

# One undamped Newton step from a converged point so the Dual (derivative)
# parts of `u` satisfy the linearized residual exactly — implicit function
# theorem. A no-op for Float64 problems.
function _ift_polish(method, f, u, F)
    _is_dual(u) || return (u, F)
    p = _newton_step(_jacobian(method, f, u, F), F)
    isnothing(p) && return (u, F)
    u2 = u .+ p
    (u2, f(u2))
end

# ── Newton with backtracking ─────────────────────────────────────────────────

function nlsolve(alg::NewtonBacktracking, f, u0::AbstractVector;
                 tol::Real, maxiters::Int)
    u = copy(u0)
    F = f(u)
    _finite(F) || error("nlsolve: residual is not finite at the initial guess")
    nF = _nval(F)
    iters = 0
    while nF >= tol && iters < maxiters
        iters += 1
        p = _newton_step(_jacobian(alg.jacobian, f, u, F), F)
        isnothing(p) && break                       # singular Jacobian
        λ = 1.0
        accepted = false
        for _ in 0:alg.max_backtracks
            u_try = u .+ λ .* p
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
    converged && ((u, F) = _ift_polish(alg.jacobian, f, u, F))
    NLResult(u, F, iters, converged)
end

# ── Powell dogleg trust region ───────────────────────────────────────────────

function _dogleg(J, F, Δ)
    pN = _newton_step(J, F)
    if !isnothing(pN) && _nval(pN) <= Δ
        return pN                                   # Newton step fits
    end
    g  = J' * F                                     # gradient of ½‖F‖²
    ng = _nval(g)
    ng == 0 && return zero(g)                       # stationary: no progress possible
    Jg = J * g
    pC = -(dot(g, g) / dot(Jg, Jg)) .* g            # Cauchy point
    if isnothing(pN) || _nval(pC) >= Δ
        return -(Δ / ng) .* g                       # steepest descent to the boundary
    end
    # Dogleg: pC + τ (pN − pC) with ‖·‖ = Δ, τ ∈ [0, 1]
    d = pN .- pC
    a = dot(d, d); b = 2 * dot(pC, d); c = dot(pC, pC) - Δ^2
    τ = (-b + sqrt(b^2 - 4a * c)) / (2a)
    pC .+ τ .* d
end

function nlsolve(alg::DoglegTrustRegion, f, u0::AbstractVector;
                 tol::Real, maxiters::Int)
    u = copy(u0)
    F = f(u)
    _finite(F) || error("nlsolve: residual is not finite at the initial guess")
    nF   = _nval(F)
    Δ    = alg.initial_radius * max(_nval(u), 1.0)
    Δmax = alg.max_radius * max(_nval(u), 1.0)
    J    = nF >= tol ? _jacobian(alg.jacobian, f, u, F) : nothing
    iters = 0
    while nF >= tol && iters < maxiters
        iters += 1
        p  = _dogleg(J, F, Δ)
        np = _nval(p)
        np == 0 && break
        pred = nF^2 - _nval(F .+ J * p)^2           # predicted reduction of ‖F‖²
        u_try = u .+ p
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
            nF >= tol && (J = _jacobian(alg.jacobian, f, u, F))
        end
        Δ < 1e-12 * max(_nval(u), 1.0) && break     # radius collapsed: stalled
    end
    converged = nF < tol
    converged && ((u, F) = _ift_polish(alg.jacobian, f, u, F))
    NLResult(u, F, iters, converged)
end
