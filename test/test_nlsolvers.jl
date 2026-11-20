using Test
using GasCycle
using ForwardDiff

const _nlsolve = GasCycle.nlsolve

@testset "Built-in nonlinear solvers" begin
    # Classic test problems with known roots.
    rosenbrock(u) = [1 - u[1], 10 * (u[2] - u[1]^2)]
    powell_badly_scaled(u) = [1e4 * u[1] * u[2] - 1, exp(-u[1]) + exp(-u[2]) - 1.0001]
    trig(u) = [cos(u[1]) - u[2], u[1] - sin(u[2]) - 0.5]

    @testset "$(nameof(typeof(alg))) / $jac" for alg in (NewtonBacktracking, DoglegTrustRegion),
                                                jac in (:forwarddiff, :finitediff)
        a = alg(jacobian = jac)

        r = _nlsolve(a, rosenbrock, [-1.2, 1.0]; tol = 1e-10, maxiters = 100)
        @test r.converged
        @test r.u ≈ [1.0, 1.0] atol = 1e-8

        r = _nlsolve(a, trig, [0.0, 0.0]; tol = 1e-12, maxiters = 100)
        @test r.converged
        @test maximum(abs, trig(r.u)) < 1e-12
    end

    @testset "trust region on a badly scaled problem" begin
        r = _nlsolve(DoglegTrustRegion(jacobian = :forwarddiff), powell_badly_scaled,
                     [0.0, 1.0]; tol = 1e-10, maxiters = 500)
        @test r.converged
        @test r.u ≈ [1.098159e-5, 9.106146] rtol = 1e-5
    end

    @testset "already converged initial guess" begin
        r = _nlsolve(NewtonBacktracking(), u -> u .- 2.0, [2.0]; tol = 1e-8, maxiters = 10)
        @test r.converged
        @test r.nsteps == 0
    end

    @testset "singular / unsolvable reports failure" begin
        r = _nlsolve(NewtonBacktracking(), u -> [u[1]^2 + 1.0], [0.5]; tol = 1e-10, maxiters = 50)
        @test !r.converged
        r = _nlsolve(DoglegTrustRegion(jacobian = :forwarddiff), u -> [u[1]^2 + 1.0], [0.5];
                     tol = 1e-10, maxiters = 50)
        @test !r.converged
    end

    # Derivatives through the solve: root of g(u, p) = 0 differentiated w.r.t. p
    # must match the implicit function theorem du/dp = -g_u⁻¹ g_p.
    @testset "AD through the solve (IFT), $(nameof(typeof(alg))) / $(alg.jacobian)" for alg in
            (NewtonBacktracking(), NewtonBacktracking(jacobian = :finitediff),
             DoglegTrustRegion(jacobian = :forwarddiff), DoglegTrustRegion())
        function root(p)
            g(u) = [u[1]^3 + u[1] - p[1], u[2] - p[2] * u[1]^2]
            _nlsolve(alg, g, [one(eltype(p)), one(eltype(p))]; tol = 1e-12, maxiters = 50).u
        end
        p = [3.0, 0.5]
        J = ForwardDiff.jacobian(root, p)
        u = root(p)
        gu = [3u[1]^2+1  0.0; -2p[2]*u[1]  1.0]
        gp = [-1.0 0.0; 0.0 -u[1]^2]
        @test J ≈ -(gu \ gp) rtol = 1e-10
    end
end
