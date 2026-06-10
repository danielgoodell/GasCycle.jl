# ── Design-mode solve ────────────────────────────────────────────────────────

function _solve_design!(net::FlowNetwork, back_edges;
                        tol::Float64, maxiter::Int, verbose::Bool)
    if isempty(back_edges)
        # No circular dependencies — the initial one_pass! is the solution.
        return SolveResult(net, :success, 1, 0.0)
    end

    # Back-edge Newton: z = [Tt₁, Pt₁, Tt₂, Pt₂, …].
    # Residual: each source element's computed outlet must equal the seeded z.
    # Normalise by z so residuals are dimensionless (tol is fractional).
    #
    # OUT-OF-PLACE form avoids the nested-Dual pre-allocation conflict when an
    # outer ForwardDiff pass wraps solve!.
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

    SolveResult(net, converged ? :success : :failed, iters, rn)
end
