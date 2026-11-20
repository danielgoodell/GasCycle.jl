"""
Assemblies: building a two-unit power system from reusable blocks.

Builder functions return `Assembly` objects — a turbine-alternator-compressor
(TAC), and a power conversion unit (PCU) that nests a TAC with its
recuperator, heater, and water cooler, plus its own rectifier and load. Two
PCUs from the same builder go into one network; names are prefixed
("PCU1.TAC.Comp"), each unit balances its own shaft through solver pairs
declared inside its assembly, and the loss ledger can be read per unit.

Copy and adapt the builders: they are templates, not fixed API.
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using GasCycle

const DATA = joinpath(@__DIR__, "..", "data")

"""Turbine–alternator–compressor on one shaft; gas enters at the alternator."""
function tac_assembly(; N = 36000.0, PR = 1.9, η_c = 0.80, η_t = 0.87,
                        P_exit = 163.4e3 / ((1 - 0.022) * (1 - 0.01)))
    tac   = Assembly("TAC")
    alt   = Alternator("Alt"; V_rated = 400.0, N_rated = N, poles = 4, η = 0.95)
    comp  = Compressor("Comp"; PR = PR, η_poly = η_c)
    turb  = Turbine("Turb"; η_poly = η_t, mode = :pressure_closure, P_exit = P_exit)
    shaft = Shaft("Shaft"; N = N)
    add!(tac, alt, comp, turb)
    connect!(tac, alt => comp)                     # alternator cooled by compressor-inlet gas
    add_shaft!(tac, shaft; drives = [comp, alt], driven_by = turb)
    expose!(tac, :inlet, alt)
    expose!(tac, :comp_out, comp, :outlet)
    expose!(tac, :turb_in, turb, :inlet)
    expose!(tac, :outlet, turb)
    expose!(tac, :elec, alt)
    tac
end

"""Power conversion unit: TAC + recuperator + heater + water cooler + rectifier and load."""
function pcu_assembly(name; TIT = 1144.0, water_T = 290.0)
    gas   = IdealGasFluid(M_molar = 83.8)
    water = ConstantPropertyLiquid(joinpath(DATA, "H2O.fpt"))
    pcu   = Assembly(name)
    tac   = tac_assembly()
    recup = HeatExchanger("Recup"; ε = 0.95, dPqP_hot = 0.022, dPqP_cold = 0.011)
    heat  = HeatSource("Heater"; TtExit = TIT, dPqP = 0.027)
    cool  = HeatExchanger("Cooler"; UA = 400.0, dPqP_hot = 0.01, dPqP_cold = 0.005)
    add!(pcu, tac, recup, heat, cool)
    connect_port!(pcu, tac, :comp_out, recup, :inlet)
    connect!(pcu, recup => heat)
    connect_port!(pcu, heat, :outlet, tac, :turb_in)
    add_hx_pair!(pcu, recup; hot = tac)
    connect_port!(pcu, recup, :hot_outlet, cool, :hot_inlet)
    connect_port!(pcu, cool, :hot_outlet, tac, :inlet; back_edge = true)
    set_boundary!(pcu, cool, :cold_inlet; Pt = 600e3, Tt = water_T, W = 1.0, fluid = water)
    # Initial guess for this unit's loop (the loop closure overrides it).
    set_boundary!(pcu, tac["Alt"], :inlet; Pt = 163.4e3, Tt = 300.0, W = 0.6, fluid = gas)

    rect = Rectifier("Rect"; η = 0.97)
    user = ElectricalLoad("User"; P = 5e3)
    connect_electrical!(pcu, tac => rect => user)
    shaft = tac["Shaft"]
    add_free!(pcu, user, :P)                                   # each unit delivers what it can
    add_target!(pcu, n -> power_balance(shaft), 0.0; scale = 1e3, name = "shaft balance")
    pcu
end

pcu1 = pcu_assembly("PCU1")
pcu2 = pcu_assembly("PCU2"; TIT = 1100.0)                      # a slightly derated second unit
net = FlowNetwork()
add!(net, pcu1)
add!(net, pcu2)
set_state!(net, pcu1["TAC"]; Pt = 163.4e3, Tt = 300.0, W = 0.6, fluid = IdealGasFluid(M_molar = 83.8))

sol = solve!(net)
summary(sol)
println()
electrical_report(sol)
println()
for p in (pcu1, pcu2)
    println(rpad(p.name, 6), "delivers ", round(p["User"].P / 1e3, digits = 3), " kW   ",
            "compressor inlet ", round(p["TAC.Comp"].inlet[].Tt, digits = 2), " K")
end
println()
L = loss_ledger(sol)
println("Exergy losses by unit:")
for (unit, x) in losses_by_assembly(L, net)
    println("  ", rpad(unit, 10), round(x / 1e3, digits = 3), " kW")
end
println("  … and one level deeper:")
for (unit, x) in losses_by_assembly(L, net; depth = 2)
    println("  ", rpad(unit, 10), round(x / 1e3, digits = 3), " kW")
end
