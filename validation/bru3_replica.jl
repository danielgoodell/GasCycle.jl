"""
Rung-2 of the NPSS cross-validation ladder (see PLAN.md): full BRU3.mdl
replica solving the *NPSS-equivalent* constraint set, not the textbook one.

From solver.bad + the .mdl, the (successful-run) NPSS system is:

  unknowns:  TIT, turbine PR, tear (Tt, Pt), oil inlet Tt, pump Pout
  residuals: comp-inlet closure (Tt, Pt), tear closure (Tt, Pt),
             oil-loop Pt closure, shaft power balance vs HPX at 36 krpm

i.e. **TIT is an output** (whatever delivers the HPX load), turbine PR is
set by pressure closure, and the oil temperature floats so the gas returns
to the compressor at exactly 540 °R.  The system is triangular, so this
script solves it as: PR_t analytic (pressure closure) → 1-D bisection on
TIT (shaft balance) → oil Tt direct (cold-end closure).

Fluid: CEAT-equivalent.  CEAT.fpt is a live CEA passthrough with He weight
fraction 0.0181 → M = 83.328 g/mol; CEA noble-gas thermo is exactly ideal
monatomic, so IdealGasFluid replicates it to machine precision.

HPX = 10.9/.92 + parasitics = 13.4198 in the .mdl, with ambiguous units
("horsepower" comment, kW-looking magnitude).  Both interpretations are
solved; the heater-Q and station anchors decide which one the successful
NPSS run used.
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using GasCycle
using Printf

# ── CEAT-equivalent fluid ─────────────────────────────────────────────────────
M_He, M_Xe = 4.002602, 131.293
w_He  = 0.0181                                  # CEAT.fpt: Wreac1
M_mix = 1 / (w_He / M_He + (1 - w_He) / M_Xe)   # 83.328 g/mol
fluid = IdealGasFluid(M_molar = M_mix)
oil   = ConstantPropertyLiquid(joinpath(@__DIR__, "..", "Oil.fpt"))  # Dow 200

# ── BRU3.mdl parameters ───────────────────────────────────────────────────────
T0, P0 = R_to_K(540.0), psia_to_Pa(23.7)
W      = lbps_to_kgps(1.32)
W_oil  = lbps_to_kgps(0.14)
P_oil  = psia_to_Pa(70 + 14.2 + 6)              # pump outlet

PR_c, η_c, η_t = 1.9, 0.80, 0.87
ε_recup        = 0.95
ε_sink         = 0.946
dP_rc, dP_rh   = 0.011, 0.022                   # recup cold/hot
dP_heat        = 0.027
dP_sink_gas    = (14.45 - 14.2) / 14.45         # ≈ 0.0173
dP_sink_oil    = 0.005
frac_bld       = 0.02
T_bld          = R_to_K(459.0 + 100.0)          # bleed forced to 559 °R
N_des          = 36_000.0

# Turbine exit pressure for exact gas-loop closure (NPSS Start.dep_Pt):
P_exit = P0 / ((1 - dP_rh) * (1 - dP_sink_gas))

# ── Network (BRU3 topology incl. bleed and cold end) ──────────────────────────
comp   = Compressor("Comp"; PR = PR_c, η_poly = η_c, η_type = :isentropic)
bsplit = Splitter("BldSplit"; fracs = [1 - frac_bld, frac_bld])
recup  = HeatExchanger("Recup"; ε = ε_recup, dPqP_hot = dP_rh, dPqP_cold = dP_rc)
bcool  = HeatSource("BldCool"; TtExit = T_bld, mode = :fixed_TtExit, dPqP = 0.0)
heater = HeatSource("Heater"; TtExit = R_to_K(2060.0), dPqP = dP_heat)
bmix   = Mixer("BldMix")
turb   = Turbine("Turb"; mode = :pressure_closure, P_exit = P_exit,
                 η_poly = η_t, η_type = :isentropic)
sink   = HeatExchanger("HeatSinkHx"; ε = ε_sink,
                       dPqP_hot = dP_sink_gas, dPqP_cold = dP_sink_oil)
shaft  = Shaft("Shaft"; N = N_des)

net = FlowNetwork()
add!(net, comp, bsplit, recup, bcool, heater, bmix, turb, sink)
connect!(net, comp => bsplit => recup => heater => bmix => turb)
connect_port!(net, bsplit, :bleed_outlet, bcool, :inlet)
connect_port!(net, bcool, :outlet, bmix, :bleed_inlet)
add_shaft!(net, shaft; drives = comp, driven_by = turb)
add_hx_pair!(net, recup; hot = turb)
connect_port!(net, recup, :hot_outlet, sink, :hot_inlet)
set_state!(net, comp; Pt = P0, Tt = T0, W = W, fluid = fluid)
set_boundary!(net, sink, :cold_inlet;
              Pt = P_oil, Tt = R_to_K(527.0), W = W_oil, fluid = oil)

solve_at!(TIT_K) = (heater.TtExit = TIT_K; solve!(net))

# ── Shaft balance: bisect TIT so net power = HPX ──────────────────────────────
function tit_for_load(P_load_W)
    lo, hi = R_to_K(1500.0), R_to_K(2600.0)
    for _ in 1:60
        mid = 0.5 * (lo + hi)
        sol = solve_at!(mid)
        sol.status == :success || error("inner solve failed at TIT=$(K_to_R(mid)) °R")
        net_power(sol) < P_load_W ? (lo = mid) : (hi = mid)
    end
    0.5 * (lo + hi)
end

HPX = 10.9 / 0.92 + 0.65 + 0.236 + 0.272 + 0.064 + 0.1 + 0.05 + 0.2   # 13.4198

function report(label, P_load_W)
    TIT = tit_for_load(P_load_W)
    sol = solve_at!(TIT)

    # Oil inlet temperature that closes the cold end at exactly 540 °R
    # (NPSS Dow200Start.ind_Tt ↔ Start.dep_Tt).  Gas is C_min in the sink HX.
    Th_in  = sink.hot_inlet[].Tt
    T_oil  = Th_in - (Th_in - T0) / ε_sink
    set_boundary!(net, sink, :cold_inlet;
                  Pt = P_oil, Tt = T_oil, W = W_oil, fluid = oil)
    sol = solve_at!(TIT)

    println("\n=== $label:  shaft load = $(round(P_load_W/1e3, digits=3)) kW ===")
    @printf("%-38s %10s %12s\n", "quantity", "GasCycle", "BRU3.mdl")
    @printf("%-38s %10.1f %12s\n", "TIT [°R]  (NPSS independent)",
            K_to_R(TIT), "2060 (init)")
    @printf("%-38s %10.2f %12s\n", "comp outlet Tt [°R]",
            K_to_R(comp.outlet[].Tt), "737")
    @printf("%-38s %10.2f %12s\n", "comp outlet Pt [psia]",
            Pa_to_psia(comp.outlet[].Pt), "45.03")
    @printf("%-38s %10.2f %12s\n", "turb inlet Pt [psia]",
            Pa_to_psia(turb.inlet[].Pt), "43.2")
    @printf("%-38s %10.4f %12s\n", "turbine PR (pressure closure)",
            pressure_ratio(turb), "1.75 (init)")
    @printf("%-38s %10.2f %12s\n", "turb outlet Tt [°R] (tear)",
            K_to_R(turb.outlet[].Tt), "1701 (init)")
    @printf("%-38s %10.2f %12s\n", "turb outlet Pt [psia] (tear)",
            Pa_to_psia(turb.outlet[].Pt), "24.69")
    @printf("%-38s %10.2f %12s\n", "recup hot outlet Tt [°R]",
            K_to_R(recup.hot_outlet[].Tt), "786")
    @printf("%-38s %10.2f %12s\n", "sink gas outlet Tt [°R] (closure)",
            K_to_R(sink.hot_outlet[].Tt), "540")
    @printf("%-38s %10.2f %12s\n", "oil inlet Tt [°R] (NPSS independent)",
            K_to_R(T_oil), "527 (init)")
    @printf("%-38s %10.2f %12s\n", "heater Q [kW]",
            heater.Q / 1e3, "~33.1")
    @printf("%-38s %10.2f %12s\n", "net shaft power [kW]",
            net_power(sol) / 1e3, "= HPX")
    nothing
end

# Four cases: {2% bleed, no bleed} × {HPX as kW, HPX as hp}.  The no-bleed
# rows test the hypothesis that the successful NPSS run's interstage bleed
# moved no flow (the .mdl header's open issue).
for (blabel, fracs) in (("2% bleed", [1 - frac_bld, frac_bld]),
                        ("no bleed", [1.0, 0.0]))
    bsplit.fracs = fracs
    report("$blabel, HPX as kW", HPX * 1e3)          # 13.420 kW
    report("$blabel, HPX as hp (NPSS native)", hp_to_W(HPX))  # 10.008 kW
end

println("""

Verdict (2026-06-10, see RESULTS.md): "no bleed, HPX as kW" matches every
.mdl anchor simultaneously — comp 737.58/45.03 (≈737/45.03), tear Pt 24.66
(24.69), recup hot out 785.0 (786), oil 526.0 (527), heater 33.75 kW
(~33.1), TIT 2045.8 °R (2060 design).  The successful NPSS run had no
effective bleed flow and a 13.42 kW shaft extraction.""")
