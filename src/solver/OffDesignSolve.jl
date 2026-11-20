# ── General solve: off-design unknowns + user free variables + back-edges ────
#
# Unknown vector  u = [element/shaft indeps; free variables; back-edge (Tt, Pt)]
# Residual vector F = [element/shaft residuals; targets;     back-edge closure]
#
# Used for off-design (map coordinates, shaft speed) and for any solve with
# add_free!/add_target! pairs, in design or off-design mode. The solve runs in
# normalized unknowns v = u ./ uref so map coordinates, shaft speeds, flows and
# thermodynamic states have similarly scaled columns.

_ref_scale(x) = (a = abs(Float64(_primal(x)))) > 0 ? a : 1.0

function _solve_general!(net::FlowNetwork, x0::AbstractVector, slices, back_edges;
                         tol::Float64, maxiter::Int, verbose::Bool,
                         alg::NonlinearBackend = DoglegTrustRegion())
    fvs  = net.free_vars
    n_x  = length(x0)
    n_f  = length(fvs)
    n_be = length(back_edges)

    f0 = _promote_vec(Real[fv.get() for fv in fvs])
    z0 = _back_edge_z0(net, back_edges)
    u0 = _promote_vec(Real[x0; f0; z0])

    uref = Float64[_ref_scale(u) for u in u0]
    for (k, fv) in enumerate(fvs)
        isnan(fv.ref) || (uref[n_x + k] = fv.ref)
    end

    # Box bounds (normalized); only free variables are bounded.
    bounded = any(fv -> isfinite(fv.lower) || isfinite(fv.upper), fvs)
    lower = upper = nothing
    if bounded
        lower = fill(-Inf, length(u0)); upper = fill(Inf, length(u0))
        for (k, fv) in enumerate(fvs)
            lower[n_x + k] = fv.lower / uref[n_x + k]
            upper[n_x + k] = fv.upper / uref[n_x + k]
        end
    end

    ix = 1:n_x
    iff = n_x+1 : n_x+n_f
    iz = n_x+n_f+1 : n_x+n_f+2n_be

    function apply!(u)
        _scatter_indeps!(net, view(u, ix), slices)
        for (k, fv) in enumerate(fvs)
            fv.set!(u[n_x + k])
        end
        one_pass!(net, n_be == 0 ? nothing : u[iz])
    end

    function residual(v)
        u = v .* uref
        apply!(u)
        z = u[iz]
        r_el = _collect_residuals(net)
        r_t  = Real[(t.f(net) - t.value) / t.scale for t in net.targets]
        r_be = Real[]
        for (i, edge) in enumerate(back_edges)
            src_out = _get_outlet(net.elements[edge.src], edge.src_port)
            if isnothing(src_out)
                push!(r_be, 0.0, 0.0)
            else
                push!(r_be, (src_out[].Tt - z[2i-1]) / z[2i-1],
                            (src_out[].Pt - z[2i])   / z[2i])
            end
        end
        _promote_vec(Real[r_el; r_t; r_be])
    end

    nls = nlsolve(alg, residual, u0 ./ uref; tol = tol, maxiters = maxiter,
                  lower = lower, upper = upper)

    rn        = norm(nls.resid)
    converged = nls.converged
    iters     = nls.nsteps

    # Final pass with the converged unknowns to restore all element fields.
    apply!(nls.u .* uref)

    if verbose
        status = converged ? "converged" : "DID NOT CONVERGE"
        println("Solver (Newton, $(length(u0)) unknowns): $status in $iters iters, " *
                "|F| = $(round(Float64(_primal(rn)), sigdigits=3))")
        for fv in fvs
            println("  free   $(rpad(fv.name, 24)) = $(round(Float64(_primal(fv.get())), sigdigits=8))")
        end
        for t in net.targets
            println("  target $(rpad(t.name, 24)) = $(round(Float64(_primal(t.f(net))), sigdigits=8))" *
                    "  (want $(round(Float64(_primal(t.value)), sigdigits=8)))")
        end
    end
    converged || @warn "GasCycle Newton did not converge after $maxiter iterations (|F|=$(_primal(rn)))"

    SolveResult(net, converged ? :success : :failed, iters, rn)
end
