# ── Off-design solve ─────────────────────────────────────────────────────────

function _solve_offdesign!(net::FlowNetwork, x0::AbstractVector, slices,
                           n_x::Int, n_res::Int, back_edges, n_be::Int;
                           tol::Float64, maxiter::Int, verbose::Bool)
    # Unknown vector u = [element/shaft indeps; back-edge (Tt, Pt) pairs].
    # The solve runs in normalized unknowns v = u ./ uref so map coordinates,
    # shaft speeds, and thermodynamic states have similarly scaled columns.
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
