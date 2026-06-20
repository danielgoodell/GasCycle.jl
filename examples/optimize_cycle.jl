"""
Gradient-based DESIGN OPTIMIZATION of a recuperated He-Xe Brayton cycle.

The "killer feature" demo. GasCycle propagates ForwardDiff through the entire
cycle solve -- including the Newton loop-closure -- so the derivatives of any
performance metric with respect to the design variables are available *exactly*
and cheaply. Not just the gradient: the full Hessian computes too (ForwardDiff
over ForwardDiff, nested two levels deep, straight through the nonlinear solve).
We hand those exact derivatives to a standard optimizer (Optim.jl) and let it
size the cycle. This is the advantage pyCycle gets from hand-coded analytic
derivatives -- here it comes for free from AD, and it extends to second order,
which a finite-difference tool like NPSS cannot do in practice.

Problem (a textbook closed-Brayton trade):

    maximize   net shaft power  W_net(PR_comp, TIT)
    over       PR_comp in [1.2, 5.0]
               TIT     in [900, 1150] K     (1150 K = turbine material limit)
    fixed      recuperator eps = 0.90, inlet 400 K / 500 kPa, 10 kg/s He-Xe

Net power increases monotonically with TIT (so the optimum pegs TIT at its
material limit) but has a genuine *interior* maximum in PR_comp: too little PR
and the turbine does little work; too much and the compressor eats it back.
The optimizer must discover that interior PR while driving TIT to its constraint.

We solve it two ways -- first-order (L-BFGS, exact gradient) and second-order
(interior-point Newton, exact gradient + Hessian) -- then VALIDATE both against
a brute-force grid sweep.

NOTE: Optim.jl is intentionally NOT a dependency of GasCycle (the core package
only needs ForwardDiff). This example assumes Optim is available in the active
or stacked environment (`] add Optim`).
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using GasCycle
using ForwardDiff
using Optim
using Printf

# ── Fixed operating conditions ────────────────────────────────────────────────
const FLUID = HeXeIdealGas(0.47)           # 47 mol% He CBC mixture
const T0, P0, W = 400.0, 500e3, 10.0        # compressor inlet T/P, mass flow
const EPS_RECUP = 0.90

# Count every cycle solve so we can compare optimizer vs grid fairly.
const NSOLVE = Ref(0)

# ── Build + solve the recuperated cycle for design vars x = [PR_comp, TIT] ─────
# Parametric in eltype(x) so ForwardDiff Duals flow straight through solve!.
function brayton_solve(x::AbstractVector)
    NSOLVE[] += 1
    PR_comp, TIT = x[1], x[2]
    comp   = Compressor("Comp";  PR=PR_comp, η_poly=0.88)
    recup  = HeatExchanger("Recup"; ε=EPS_RECUP, dPqP_hot=0.01, dPqP_cold=0.01)
    heater = HeatSource("Heat";  TtExit=TIT, dPqP=0.02)
    turb   = Turbine("Turb"; mode=:pressure_closure,
                              P_exit=P0 / ((1 - 0.01) * (1 - 0.01)), η_poly=0.90)
    net = FlowNetwork()
    add!(net, comp, recup, heater, turb)
    connect!(net, comp => recup => heater => turb => comp)
    add_hx_pair!(net, recup; hot=turb)
    set_state!(net, comp; Pt=P0, Tt=T0, W=W, fluid=FLUID)
    solve!(net; maxiter=200)
end

net_power_kW(x) = net_power(brayton_solve(x)) / 1000

# ── Design space (physical) and a normalized [0,1]^2 working space ────────────
# Optimizing in normalized coordinates puts PR (~O(1)) and TIT (~O(1e3)) on the
# same scale, which conditions the optimizer well.
const LB = [1.2,  900.0]
const UB = [5.0, 1150.0]
phys(u)  = LB .+ u .* (UB .- LB)            # u in [0,1]^2  ->  physical x
norml(x) = (x .- LB) ./ (UB .- LB)

# Objective handed to Optim: MINIMIZE negative power, in normalized coords.
# The gradient and Hessian closures are pure ForwardDiff -- exact, no finite diff.
fobj(u)     = -net_power_kW(phys(u))
gobj!(G, u) = ForwardDiff.gradient!(G, fobj, u)
hobj!(H, u) = ForwardDiff.hessian!(H, fobj, u)

# ── Sanity: exact AD gradient vs finite differences at the start point ────────
x_start = [2.5, 1000.0]
println("="^72)
println("Exact AD gradient vs finite difference   (at PR = 2.5, TIT = 1000 K)")
∇phys = ForwardDiff.gradient(net_power_kW, x_start)
hh = [1e-4, 1e-2]
fd = [(net_power_kW(x_start .+ [i == 1 ? hh[1] : 0.0, i == 2 ? hh[2] : 0.0]) -
       net_power_kW(x_start .- [i == 1 ? hh[1] : 0.0, i == 2 ? hh[2] : 0.0])) / (2hh[i])
      for i in 1:2]
@printf("  dW/dPR  : AD = %+9.4f   FD = %+9.4f   kW per unit PR\n", ∇phys[1], fd[1])
@printf("  dW/dTIT : AD = %+9.5f   FD = %+9.5f   kW per K\n",       ∇phys[2], fd[2])

# ── Optimizer A: first order (L-BFGS + box constraints, exact gradient) ───────
u0 = norml(x_start)
NSOLVE[] = 0
resA = optimize(fobj, gobj!, zeros(2), ones(2), u0,
                Fminbox(LBFGS()), Optim.Options(g_tol = 1e-4))
xA, WA, nA = phys(Optim.minimizer(resA)), -Optim.minimum(resA), NSOLVE[]

# ── Optimizer B: second order (interior-point Newton, exact grad + Hessian) ───
df  = TwiceDifferentiable(fobj, gobj!, hobj!, u0)
dfc = TwiceDifferentiableConstraints(zeros(2), ones(2))   # box constraints
NSOLVE[] = 0
resB = optimize(df, dfc, u0, IPNewton())
xB, WB, nB = phys(Optim.minimizer(resB)), -Optim.minimum(resB), NSOLVE[]

println("="^72)
println("OPTIMIZER RESULTS   (TIT* driven to the 1150 K material limit; PR* interior)")
@printf("  A  L-BFGS  (exact grad)        : PR* = %.4f  TIT* = %.1f  W* = %.3f kW  (%d solves)\n",
        xA[1], xA[2], WA, nA)
@printf("  B  IPNewton (exact grad+Hess)  : PR* = %.4f  TIT* = %.1f  W* = %.3f kW  (%d solves)\n",
        xB[1], xB[2], WB, nB)
∇opt = ForwardDiff.gradient(net_power_kW, xB)
@printf("  grad at optimum: dW/dPR = %+.4f (~0, interior),  dW/dTIT = %+.4f (>0, pinned)\n",
        ∇opt[1], ∇opt[2])

# ── Validate against a brute-force grid sweep ─────────────────────────────────
println("="^72)
println("VALIDATION: brute-force grid sweep")
nPR, nT = 60, 26
best = (-Inf, 0.0, 0.0)
NSOLVE[] = 0
for PR in range(LB[1], UB[1], nPR), TIT in range(LB[2], UB[2], nT)
    global best
    sol = brayton_solve([PR, TIT])
    sol.status == :success || continue
    Wk = net_power(sol) / 1000
    Wk > best[1] && (best = (Wk, PR, TIT))
end
nGrid = NSOLVE[]
@printf("  grid best : W_net = %.3f kW  at PR = %.3f, TIT = %.1f K   (%d solves)\n",
        best[1], best[2], best[3], nGrid)
@printf("  optimizer : W_net = %.3f kW  at PR = %.4f (exact)            (%d solves)\n",
        WB, xB[1], nB)
@printf("  agreement : dW = %.4f kW (%.4f%%) -- and the optimizer's PR is exact,\n",
        abs(WB - best[1]), 100 * abs(WB - best[1]) / best[1])
println("              not snapped to the grid's discrete PR spacing.")

# ── Why this matters: scaling with the number of design variables ─────────────
println("="^72)
println("WHY IT MATTERS -- a 2-variable grid is cheap; real cycle MDO is not.")
println("Cost of a grid sweep at 40 points/axis vs gradient-based optimization:")
println("  design vars   grid solves (40^k)     gradient-based (forward-mode AD)")
for k in (2, 4, 6, 8)
    @printf("     %d          %1.1e              ~%d solves\n", k, Float64(40)^k, 60 * k)
end
println("""
The grid (or any DOE / gradient-free sweep) explodes as resolution^k; this is
the wall pyCycle was built to get past. Gradient-based optimization with exact
AD derivatives stays tractable. With forward-mode AD the gradient costs ~k
sweeps; switch to reverse-mode/adjoint and that cost becomes independent of k
entirely -- the natural next step once the design space gets large.""")
println("="^72)
