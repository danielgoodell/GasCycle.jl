"""
ForwardDiff sensitivity analysis on a recuperated He-Xe Brayton cycle.

Demonstrates exact first-order derivatives of cycle performance metrics
(net power, cycle efficiency) with respect to design parameters:
  - Compressor pressure ratio (PR_comp)
  - Recuperator effectiveness (ε_recup)
  - Turbine inlet temperature (TIT)

Uses IdealGasFluid (constant-Cp monatomic ideal gas) so all inversions
are closed-form and ForwardDiff Dual numbers propagate exactly through
the entire cycle computation.
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using GasCycle
using ForwardDiff

# ── Fluid ──────────────────────────────────────────────────────────────────────
fluid = HeXeIdealGas(0.47)   # 47 mol% He, ~39 g/mol (medium-weight CBC mix)

# ── Helper: build and solve a recuperated Brayton cycle ────────────────────────
# params = [PR_comp, ε_recup, TIT_K]
function brayton_solve(params::AbstractVector)
    PR_comp, ε_recup, TIT = params[1], params[2], params[3]

    T0  = 400.0    # compressor inlet temperature [K]
    P0  = 500e3    # compressor inlet pressure    [Pa]
    W   = 10.0     # mass flow rate               [kg/s]

    net = FlowNetwork()

    comp   = Compressor("Comp";  PR=PR_comp, η_poly=0.88)
    recup  = HeatExchanger("Recup"; ε=ε_recup, dPqP_hot=0.01, dPqP_cold=0.01)
    heater = HeatSource("Heat";  TtExit=TIT, dPqP=0.02)
    turb   = Turbine("Turb";     mode=:pressure_closure,
                                  P_exit = P0 / ((1-0.01)*(1-0.01)),
                                  η_poly = 0.90)

    add!(net, comp, recup, heater, turb)
    connect!(net, comp => recup => heater => turb => comp)
    add_hx_pair!(net, recup; hot=turb)
    set_state!(net, comp; Pt=P0, Tt=T0, W=W, fluid=fluid)

    solve!(net; maxiter=200)
end

function net_power_kW(params::AbstractVector)
    sol = brayton_solve(params)
    net_power(sol) / 1000
end

function cycle_eff_pct(params::AbstractVector)
    sol = brayton_solve(params)
    cycle_efficiency(sol) * 100
end

# ── Design point ───────────────────────────────────────────────────────────────
x0 = [2.5,   # PR_comp
      0.90,   # ε_recup
      1100.0] # TIT [K]

sol0 = brayton_solve(x0)
W_net  = net_power(sol0) / 1000
η_cyc  = cycle_efficiency(sol0) * 100
println("Design-point results:")
println("  Net shaft power   : $(round(W_net,  digits=2)) kW")
println("  Cycle efficiency  : $(round(η_cyc,  digits=2)) %")
println()

# ── ForwardDiff gradients ─────────────────────────────────────────────────────
println("ForwardDiff.gradient — net power w.r.t. [PR_comp, ε_recup, TIT]:")
∇W = ForwardDiff.gradient(net_power_kW, x0)
println("  ∂W_net/∂PR_comp   = $(round(∇W[1], digits=3)) kW per unit PR")
println("  ∂W_net/∂ε_recup   = $(round(∇W[2], digits=3)) kW per unit ε")
println("  ∂W_net/∂TIT       = $(round(∇W[3], sigdigits=3)) kW/K")
println()

println("ForwardDiff.gradient — cycle efficiency w.r.t. [PR_comp, ε_recup, TIT]:")
∇η = ForwardDiff.gradient(cycle_eff_pct, x0)
println("  ∂η/∂PR_comp       = $(round(∇η[1], digits=3)) % per unit PR")
println("  ∂η/∂ε_recup       = $(round(∇η[2], digits=3)) % per unit ε")
println("  ∂η/∂TIT           = $(round(∇η[3], sigdigits=3)) %/K")
println()

# ── Verify via finite differences ─────────────────────────────────────────────
println("Finite-difference verification (h = 1e-4 * x0):")
h = 1e-4
fd = [(net_power_kW(x0 .+ [i==j ? h*x0[j] : 0.0 for j in 1:3]) -
       net_power_kW(x0 .+ [i==j ? -h*x0[j] : 0.0 for j in 1:3])) / (2*h*x0[i])
      for i in 1:3]
labels = ["PR_comp", "ε_recup ", "TIT     "]
for i in 1:3
    ad_val = ∇W[i]
    fd_val = fd[i]
    rel_err = abs(ad_val - fd_val) / (abs(fd_val) + 1e-12) * 100
    println("  $(labels[i]): AD = $(round(ad_val, sigdigits=5)),  FD = $(round(fd_val, sigdigits=5)),  err = $(round(rel_err, sigdigits=2))%")
end
println()

# ── Jacobian: all performance metrics × all design params ─────────────────────
function perf_vector(params::AbstractVector)
    sol = brayton_solve(params)
    [net_power(sol)/1000, cycle_efficiency(sol)*100]
end

J = ForwardDiff.jacobian(perf_vector, x0)
println("Jacobian  [W_net_kW, η_cyc_pct]  ×  [PR_comp, ε_recup, TIT]:")
println("  $(round.(J; sigdigits=4))")
println()
println("(rows = outputs, cols = inputs)")
println()
println("Interpretation:")
println("  Increasing PR by 1 changes net power by $(round(J[1,1],digits=2)) kW")
println("  Increasing ε by 0.01 changes net power by $(round(J[1,2]*0.01, digits=3)) kW")
println("  Increasing TIT by 10 K changes efficiency by $(round(J[2,3]*10, digits=3)) %")
