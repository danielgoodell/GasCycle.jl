"""
Solver for the FlowNetwork.

Two modes of operation:

1. **Back-edge Newton** (n_indeps == 0):
   Used for cycles where all design parameters are specified but there are
   circular dependencies (e.g., the recuperator hot inlet depends on the
   turbine outlet which depends on the recuperator cold outlet).

2. **Newton-Raphson with finite differences** (n_indeps > 0):
   Used for off-design where shaft speed, map operating points, etc. are
   unknowns with associated residual equations.
"""

using LinearAlgebra: norm
using NonlinearSolve

include("SolveResult.jl")
include("ResidualAssembly.jl")
include("DesignSolve.jl")
include("OffDesignSolve.jl")
include("CycleMetrics.jl")

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

    # Initial pass to populate element outlet states and seed back-edge guesses.
    one_pass!(net)

    x0, slices = _collect_indeps(net)
    n_x   = length(x0)
    n_res = sum(n_residuals(el) for el in net.elements; init=0) +
            sum(n_residuals(sh) for sh in net.shafts; init=0)

    n_x == n_res || @warn "Solver: $n_x indeps but $n_res residuals — system may be over/under-determined"

    back_edges = filter(e -> e.back_edge, net.edges)
    n_be = length(back_edges)

    if n_x == 0
        return _solve_design!(net, back_edges; tol=tol, maxiter=maxiter, verbose=verbose)
    end

    _solve_offdesign!(net, x0, slices, n_x, n_res, back_edges, n_be;
                      tol=tol, maxiter=maxiter, verbose=verbose)
end
