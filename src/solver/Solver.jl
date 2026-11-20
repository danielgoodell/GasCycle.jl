"""
Solver for the FlowNetwork.

Two modes of operation:

1. **Back-edge Newton** (n_indeps == 0):
   Used for cycles where all design parameters are specified but there are
   circular dependencies (e.g., the recuperator hot inlet depends on the
   turbine outlet which depends on the recuperator cold outlet).

2. **Trust-region Newton with finite differences** (n_indeps > 0):
   Used for off-design where shaft speed, map operating points, etc. are
   unknowns with associated residual equations.

Both use the small built-in solvers in `NonlinearSolvers.jl`.
"""

using LinearAlgebra: norm

include("NonlinearSolvers.jl")
include("SolveResult.jl")
include("ResidualAssembly.jl")
include("DesignSolve.jl")
include("OffDesignSolve.jl")
include("CycleMetrics.jl")

"""
    solve!(net; tol=1e-6, maxiter=100, verbose=false, alg=nothing) -> SolveResult

Solve the flow network.

Design-mode networks whose only unknowns are back-edge [Tt, Pt] states
(recuperators, loop closures) are resolved by Newton (`NewtonBacktracking`,
ForwardDiff Jacobian). Networks with element unknowns (off-design maps, shaft
speed) or user free variables/targets (`add_free!`/`add_target!`) are solved
as one normalized system by a dogleg trust region (`DoglegTrustRegion`,
finite-difference Jacobian). Pass `alg` to override the backend.

The system must be square: element unknowns + free variables must equal
element residuals + targets.
"""
function solve!(net::FlowNetwork;
                tol::Float64  = 1e-6,
                maxiter::Int  = 100,
                verbose::Bool = false,
                alg::Union{Nothing,NonlinearBackend} = nothing)

    # Initial pass to populate element outlet states and seed back-edge guesses,
    # then start heat links / pumps at their model values and pass again.
    one_pass!(net)
    if !isempty(net.heat_links) || any(el -> el isa Pump, net.elements)
        _init_couplings!(net)
        one_pass!(net)
    end

    x0, slices = _collect_indeps(net)
    n_x   = length(x0)
    n_res = sum(n_residuals(el) for el in net.elements; init=0) +
            sum(n_residuals(sh) for sh in net.shafts; init=0) +
            length(net.heat_links)
    n_f = length(net.free_vars)
    n_t = length(net.targets)

    if n_x + n_f != n_res + n_t
        error("Solver: system is not square — $n_x element unknowns + $n_f free variables " *
              "vs $n_res element residuals + $n_t targets. Each target needs a free " *
              "variable (add_free!) and vice versa.")
    end

    back_edges = filter(e -> e.back_edge, net.edges)

    if n_x == 0 && n_f == 0
        return _solve_design!(net, back_edges; tol=tol, maxiter=maxiter, verbose=verbose,
                              alg = something(alg, NewtonBacktracking()))
    end

    _solve_general!(net, x0, slices, back_edges;
                    tol=tol, maxiter=maxiter, verbose=verbose,
                    alg = something(alg, DoglegTrustRegion()))
end
