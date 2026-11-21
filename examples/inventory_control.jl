"""
Inventory (charge-pressure) control of a closed Brayton loop.

A closed loop's power is set by how much gas is in it: at fixed speed and
turbine inlet temperature the turbomachines run at the same corrected
conditions, so flow, pressures and power scale with the charge. This example
sizes a recuperated loop with a water cooler at its design point, gives every
component a gas volume, then — off-design on scaled maps at fixed (alternator-
locked) speed — specifies the loop's inventory and lets the solver find the
flow. Sweeping the charge down shows the near-proportional power control. The
departure from exact proportionality comes from the water cooler: sized
(fixed UA) for full flow, it cools the smaller flow further, so the compressor
inlet runs colder at low charge — power falls a little more slowly than the
charge and efficiency rises slightly. (With the recuperator in `:scaled_UA`
mode its effectiveness would also rise at low flow.)
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using GasCycle
using Printf
include(joinpath(@__DIR__, "..", "test", "synthetic_maps.jl"))

const gas   = IdealGasFluid(M_molar = 83.8)
const water = ConstantPropertyLiquid(joinpath(@__DIR__, "..", "data", "H2O.fpt"))
const N = 36000.0
P1, T1g, W_des = 163.4e3, 300.0, 0.6
dh, dq = 0.022, 0.01

function loop(; od = nothing, W = W_des)
    comp = isnothing(od) ? Compressor("Comp"; PR = 1.9, η_poly = 0.80) :
                           Compressor("Comp"; η_poly = 0.80, map = od.cmap, mode = :off_design)
    turb = isnothing(od) ? Turbine("Turb"; η_poly = 0.87, mode = :pressure_closure,
                                   P_exit = P1 / ((1 - dh) * (1 - dq))) :
                           Turbine("Turb"; PR = od.PR_t, η_poly = 0.87, map = od.tmap, mode = :off_design)
    recup = HeatExchanger("Recup"; ε = 0.95, dPqP_hot = dh, dPqP_cold = 0.011)
    heat  = HeatSource("Heater"; TtExit = 1144.0, dPqP = 0.027)
    cool  = HeatExchanger("Cooler"; UA = 400.0, dPqP_hot = dq, dPqP_cold = 0.005)
    shaft = Shaft("Shaft"; N = N)
    net = FlowNetwork()
    add!(net, comp, recup, heat, turb, cool)
    connect!(net, comp => recup => heat => turb)
    add_hx_pair!(net, recup; hot = turb)
    connect_port!(net, recup, :hot_outlet, cool, :hot_inlet)
    connect_port!(net, cool, :hot_outlet, comp, :inlet; back_edge = true)
    set_boundary!(net, cool, :cold_inlet; Pt = 600e3, Tt = 290.0, W = 1.0, fluid = water)
    add_shaft!(net, shaft; drives = comp, driven_by = turb)
    set_state!(net, comp; Pt = P1, Tt = T1g, W = W, fluid = gas)
    # Gas volumes [m³]: ducting lumped onto the components.
    add_volume!(net, comp; V = 0.002)
    add_volume!(net, recup; V = 0.015, side = :cold)
    add_volume!(net, recup; V = 0.015, side = :hot)
    add_volume!(net, heat; V = 0.008)
    add_volume!(net, turb; V = 0.003)
    add_volume!(net, cool; V = 0.006, side = :hot)
    (; net, comp, turb, heat)
end

# Design point, then maps through it.
d = loop()
solve!(d.net)
s1, s3 = d.comp.inlet[], d.turb.inlet[]
od = (cmap = synthetic_comp_map(Nc_des = corrected_speed(N, s1.Tt), Wc_des = corrected_flow(s1.W, s1.Tt, s1.Pt),
                                PR_des = 1.9, eta_des = 0.80),
      tmap = synthetic_turb_map(Np_des = corrected_speed(N, s3.Tt), Wp_des = corrected_flow(s3.W, s3.Tt, s3.Pt),
                                PR_des = pressure_ratio(d.turb), eta_des = 0.87),
      PR_t = pressure_ratio(d.turb))
M_des = loop_mass(d.net, d.comp)
P_des = net_power(d.net)
println("Design: inventory $(round(M_des * 1e3, digits = 2)) g, net power $(round(P_des / 1e3, digits = 3)) kW\n")
inventory_report(d.net)
println()

@printf("%9s %10s %11s %12s %11s %10s\n", "charge", "W [kg/s]", "P1 [kPa]", "T1 [K]", "P_net [kW]", "η [%]")
for frac in (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)
    m = loop(; od)
    set_inventory!(m.net, m.comp; mass = frac * M_des)
    sol = solve!(m.net)
    s = m.comp.inlet[]
    @printf("%8.0f%% %10.4f %11.2f %12.2f %11.3f %10.2f\n", 100frac, s.W, s.Pt / 1e3, s.Tt,
            net_power(sol) / 1e3, 100 * cycle_efficiency(sol))
end
