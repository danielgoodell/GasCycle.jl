"""
Simple (unrecuperated) He-Xe Brayton cycle — design point.

Validates that cycle thermal efficiency matches the ideal Brayton formula:
  η_ideal = 1 - T_cold_in / T_hot_in  (only holds for ideal gas, η=1)

For a real cycle with polytropic efficiencies:
  η_real ≈ 1 - (T_2 / T_1) * (1 + 1/η_turb) / (1 + η_comp)
  (approximate; exact value comes from the solver)
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using GasCycle

# ── Fluid ──────────────────────────────────────────────────────────────────────
# Use ideal gas fallback (no FPT file needed for this example)
# HeXe84 ≈ M=84 g/mol → Cp ≈ 247 J/(kg·K)
fluid = HeXeIdealGas(0.47)   # ~47% He by mole → M ≈ 71.5 g/mol

# For the FPT-based version, use:
# fluid = FPTFluid(joinpath(@__DIR__, "..", "data", "HeXe84.fpt"))

# ── Cycle parameters ────────────────────────────────────────────────────────────
T_in   = 400.0    # K  compressor inlet temperature
P_in   = 500e3    # Pa compressor inlet pressure
W_flow = 10.0     # kg/s mass flow

PR_comp  = 2.5    # compressor pressure ratio
TIT      = 1100.0 # K  turbine inlet temperature (= reactor outlet)
η_comp   = 0.87   # polytropic efficiency
η_turb   = 0.90

# ── Build network ───────────────────────────────────────────────────────────────
net = FlowNetwork()

comp    = Compressor("Comp";    PR=PR_comp, η_poly=η_comp)
reactor = HeatSource("Reactor"; TtExit=TIT, dPqP=0.02)
turb    = Turbine("Turb";       mode=:pressure_closure, P_exit=P_in, η_poly=η_turb)
shaft   = Shaft("Main";         N=15000.0)

add!(net, comp, reactor, turb)
connect!(net, comp => reactor => turb => comp)    # closed loop
add_shaft!(net, shaft; drives=comp, driven_by=turb)
set_state!(net, comp; Pt=P_in, Tt=T_in, W=W_flow, fluid=fluid)

# ── Solve ───────────────────────────────────────────────────────────────────────
sol = solve!(net; verbose=true)

println("\n=== Simple Brayton Cycle Results ===")
println("Compressor PR:         $(round(comp.PR, digits=3))")
println("Turbine PR:            $(round(pressure_ratio(turb), digits=3))")
println("Compressor Tt_out:     $(round(comp.outlet[].Tt, digits=1)) K")
println("Turbine Tt_out:        $(round(turb.outlet[].Tt, digits=1)) K")
println("Net power:             $(round(net_power(sol)/1000, digits=2)) kW")
println("Heat input:            $(round(reactor.Q/1000, digits=2)) kW")
println("Cycle efficiency:      $(round(cycle_efficiency(sol)*100, digits=2)) %")

# Ideal Brayton efficiency (ideal gas, η_poly=1)
η_ideal = 1.0 - T_in / TIT
println("Ideal Brayton η:       $(round(η_ideal*100, digits=2)) %")
println("(Real < ideal due to compressor/turbine losses)")
