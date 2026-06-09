"""
Recuperated He-Xe Brayton cycle — design point using actual FPT fluid data.

Demonstrates:
  - FPTFluid loading from HeXe84.fpt
  - HeatExchanger (recuperator) with ε-NTU method
  - cycle_efficiency() improvement over unrecuperated cycle

Topology (closed loop):
  Comp → [Recup cold side] → Reactor → Turb → [Recup hot side] → Comp
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using GasCycle

# ── Fluid ──────────────────────────────────────────────────────────────────────
fpt_path = joinpath(@__DIR__, "..", "HeXe84.fpt")
fluid = isfile(fpt_path) ? FPTFluid(fpt_path) : HeXeIdealGas(0.47)
println("Using fluid: ", fluid isa FPTFluid ? fluid.name : "HeXe ideal gas")

# ── Cycle parameters ────────────────────────────────────────────────────────────
T_in     = 400.0    # K   compressor inlet (= recuperator cold outlet)
P_in     = 500e3    # Pa  compressor inlet pressure
W_flow   = 10.0     # kg/s
PR_comp  = 2.5
TIT      = 1100.0   # K   reactor outlet / turbine inlet
ε_recup  = 0.92     # recuperator effectiveness
η_comp   = 0.87
η_turb   = 0.90

# ── Build network ───────────────────────────────────────────────────────────────
net = FlowNetwork()

comp    = Compressor("Comp";    PR=PR_comp,   η_poly=η_comp)
recup   = HeatExchanger("Recup"; ε=ε_recup)
reactor = HeatSource("Reactor";  TtExit=TIT,  dPqP=0.02)
turb    = Turbine("Turb";        mode=:pressure_closure, P_exit=P_in, η_poly=η_turb)
shaft   = Shaft("Main";          N=15000.0)

add!(net, comp, recup, reactor, turb)

# Flow order: cold side of recup sits between comp and reactor;
# hot side connects turb outlet back to comp inlet.
# The FlowNetwork handles the two-port HeatExchanger via add_hx_pair!.
# For the sequential flow_order, the recup appears where the cold side is.
connect!(net, comp => recup => reactor => turb => comp)
add_shaft!(net, shaft; drives=comp, driven_by=turb)

# Wire the heat exchanger: hot stream = turbine outlet, cold stream = comp outlet
add_hx_pair!(net, recup; hot=turb, cold=comp)

set_state!(net, comp; Pt=P_in, Tt=T_in, W=W_flow, fluid=fluid)

# ── Solve ───────────────────────────────────────────────────────────────────────
sol = solve!(net; verbose=true)

println("\n=== Recuperated Brayton Cycle Results ===")
println("Compressor PR:            $(round(comp.PR,           digits=3))")
println("Turbine PR:               $(round(pressure_ratio(turb), digits=3))")
println("Comp outlet Tt:           $(round(comp.outlet[].Tt,  digits=1)) K")
println("Recup cold outlet Tt:     $(round(recup.cold_outlet[].Tt, digits=1)) K")
println("Turbine inlet Tt:         $(round(reactor.outlet[].Tt,    digits=1)) K")
println("Turbine outlet Tt:        $(round(turb.outlet[].Tt,  digits=1)) K")
println("Recup hot outlet Tt:      $(round(recup.hot_outlet[].Tt,  digits=1)) K")
println()
println("Recuperator Q:            $(round(Q_transferred(recup)/1000, digits=2)) kW")
println("Reactor heat input:       $(round(reactor.Q/1000, digits=2)) kW")
println("Net shaft power:          $(round(net_power(sol)/1000, digits=2)) kW")
println("Cycle efficiency:         $(round(cycle_efficiency(sol)*100, digits=2)) %")
println()
println("Ideal Brayton η (no recup): $(round((1 - T_in/TIT)*100, digits=2)) %")
println("(Recuperation should raise η above ideal simple-cycle)")
