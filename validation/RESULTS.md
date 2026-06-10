# NPSS Cross-Validation Results

Running log. See PLAN.md for the ladder definition. English units (°R, psia)
for direct comparison with NPSS prints.

## 2026-06-10 — Rung 0: entropy interpolation artifact found and fixed

Probing FPT property lookups at the BRU station states (`HeXe84.fpt`):
cp, γ, ρ all within ~1 % of ideal monatomic M = 83.8, but the s-table's
pressure derivative was wrong:

    ds/dlnP = −113.23 J/(kg·K)   measured at 300 K, 165 kPa (mid-cell)
            = −99.22             required (−R; real-gas correction ≪ 1 % here)

Root cause (confirmed to 4 significant figures): bilinear interpolation
**linear in P** of s ∝ −R·ln P across the coarse Pt grid (nodes at 101.6
and 198.2 kPa bracket the BRU compressor inlet at 163.4 kPa — worst case).
Predicted chord slope −R·P·ln(P₂/P₁)/(P₂−P₁) = −113.23 — identical to the
measured value. This single artifact inflated the implied isentropic
exponent to γ_eff ≈ 1.74 (> 5/3, unphysical) and produced **+14 °R at the
compressor outlet** — most of the historical "fluid model difference".

**Fix:** `FPTFluid(...; s_interp=:log_pressure)` (now the default)
interpolates the detrended σ = s + R·ln(Pt/P_ref) and restores the log term
analytically; R is fitted from the table itself (median −Δs/ΔlnP over all
node pairs: 99.27 vs 99.22 from M = 83.8). Exact at table nodes; exact in P
for ideal gas. `s_interp=:linear` keeps the legacy behavior for
apples-to-apples runs against interpolators that are linear in P (NPSS's
scheme: to be confirmed at rung 0 with the run listing).

After the fix the isentropic step at the compressor inlet implies
γ_eff = 1.666 vs table γ = 1.667. Regression tests in `test_thermo.jl`.

## 2026-06-10 — Rung 1: isolation results (GasCycle side; NPSS couts pending)

`julia validation/bru3_isolation.jl`, HeXe84.fpt, effDes (isentropic):

| Element (BRU3.mdl isolation inputs)        | GasCycle | NPSS      | anchor |
|--------------------------------------------|----------|-----------|--------|
| Compressor out (540 °R/23.7 psia, PR 1.9)  | 737.16 °R| PENDING   | ≈737 °R (in-loop) |
| Turbine out (2060 °R/43.2 psia, PR 1.7497) | 1700.83 °R| PENDING  | 1701 °R (tear value) |
| Recup cold out (737 + 1701 °R, ε 0.95)     | 1652.74 °R| PENDING  | ~1651 °R (.mdl note) |
| Recup hot out                              | 793.51 °R| PENDING   | 786 °R (in-loop) |
| HeatSinkHx gas out (786 °R + 527 °R oil)   | 540.92 °R| PENDING   | ≈540 °R (loop closure) |

Sub-°R agreement with every isolation-consistent anchor. Factor sizes at
the compressor outlet (°R): s-interp artifact +14.3 (fixed), semantics
isen−poly −6.3, HeXe84-vs-ideal −0.4.

Key inference: the BRU3.mdl 1701 °R HotStart tear value **is** the no-bleed
isolated turbine output (we get 1700.83). The .mdl pins station 6 there
even though its own 2 % bleed (cooled to 559 °R, reinjected at the turbine
inlet) lowers the effective turbine inlet to ~2030 °R and hence the outlet
to ~1673 °R — exactly the .mdl header's open issue ("look at Interstage
Bleed Port").

## 2026-06-10 — Rung 2 preview: full loop (examples/bru_10kw.jl, isentropic)

| Quantity        | GasCycle | NPSS/anchor | status |
|-----------------|----------|-------------|--------|
| Comp outlet     | 737.0 °R | ≈737 °R     | matches at print precision |
| Turb outlet     | 1673 °R  | 1701 °R     | Δ = bleed premix vs no-bleed tear (see above) |
| Recup cold out  | 1627 °R  | ~1651 °R    | cascades from turb outlet Δ |
| Heater input    | 34.95 kW | ~33.1 kW    | cascades from recup cold out Δ |
| Net shaft power | 13.12 kW | "13.4" (HPX) | HPX hp-vs-kW ambiguity pending |
| Est. electrical | 10.6 kW  | 10.5 kW design | — |

Everything left on the list traces to (a) the bleed/tear treatment — needs
the NPSS run listing to see what NPSS actually converges to, and ideally a
no-bleed run of both models, and (b) the CEAT-vs-HeXe84 table at hot-end
conditions — bounded small by the turbine isolation row but only directly
checkable with the artifacts in PLAN.md.

## 2026-06-10 — NPSS artifacts received: CEAT.fpt and solver.bad

**CEAT.fpt is not a table.** It is a live passthrough: every property
function sets an NPSS FlowStation (He weight fraction Wreac1 = 0.0181) and
calls the CEA equilibrium package. Consequences:

- CEA noble-gas thermo is exactly ideal monatomic (Cp = 2.5R, γ = 5/3, no
  real-gas terms, no excitation below ~8000 K), so **CEAT is replicable to
  machine precision** by `IdealGasFluid(M_molar = 83.328)`. No table file
  needed after all.
- w_He = 0.0181 ⇒ **M = 83.328 g/mol, not 83.8**: the NPSS model's fluid is
  0.57 % off the BRU spec (and off HeXe84.fpt) in R and cp. Temperature
  ratios are unaffected (γ exactly 5/3); powers and flows scale by it.

**solver.bad is the diagnostic of a run that cannot converge.** Its active
set is 5×5 — {HeatSourceHX.Tout, HotStart Tt, HotStart Pt, Dow200Start Tt,
Pump Pout} vs {Start Tt+Pt closure, HotStart Tt+Pt closure, Dow200Start Pt
closure} — with **no turbine-PR independent and no shaft-balance
dependent**. With every PR fixed, the gas loop returns 23.798 psia against
the required 23.700 (0.41 % error, 41× the 1e-4 tolerance) with no knob
that can move it: structurally unsatisfiable. This matches the .mdl's
mid-debug state. The successful run that produced the station comments must
have had `Turb.ind_PRdes` active (declared in the .mdl): exact closure then
gives PR_t = 1.7572 and turbine exit 24.66 psia ≈ the 24.69 comment.

Confirmed picture of the successful NPSS solve (6×6): unknowns {TIT, PR_t,
tear Tt, tear Pt, oil inlet Tt, pump Pout}; residuals {comp-inlet Tt and Pt
closure, tear Tt and Pt closure, oil Pt closure, shaft balance vs HPX at
36 krpm}. TIT, turbine PR, and the oil temperature are all *outputs*.

## 2026-06-10 — Rung 2: full NPSS-equivalent replica (validation/bru3_replica.jl)

CEAT-equivalent fluid, NPSS constraint set (triangular: PR_t analytic from
pressure closure → TIT from shaft balance → oil Tt from cold-end closure).
Four cases: {2 % bleed, no bleed} × {HPX = 13.42 kW, HPX = 13.42 hp}.

**"No bleed + HPX as kW" matches every anchor simultaneously:**

| quantity | replica | .mdl anchor |
|---|---|---|
| comp outlet | 737.58 °R / 45.03 psia | ≈737 / 45.03 |
| turb inlet Pt | 43.33 psia | 43.2 |
| turbine PR | 1.7572 | 1.75 (initial) |
| tear (turb out) | 1686.5 °R / 24.66 psia | 1701 (initial) / 24.69 |
| recup hot out | 785.0 °R | 786 |
| oil inlet Tt | 526.0 °R | 527 (initial) |
| heater Q | 33.75 kW | ~33.1 |
| TIT (output!) | 2045.8 °R | 2060 design spec |

The alternatives fail decisively: HPX-as-hp drives TIT to 1812–1837 °R and
heater Q to ≤31 kW; with-bleed unbalances the recuperator (hot out 803 °R
vs 786) and shifts heater Q to 35.5 kW.

**Conclusions:**
1. The successful NPSS run's interstage bleed moved no flow — precisely the
   .mdl header's open issue ("look at Interstage Bleed Port"). The 1701 °R
   tear and 786 °R station-8 comments are no-bleed values.
2. HPX was effectively a kW quantity in that run (13.42 kW extraction).
3. Residual gaps are now ~1 °R on temperatures and ~2 % on heater Q
   (33.75 vs "33.1", a rounded comment), with TIT 2045.8 vs the 2060 spec —
   i.e. the as-built NPSS model delivers its HPX load at 14 °R below the
   BRU design TIT. Closing these last digits needs the actual NPSS output
   listing (still the top artifact request) and CEA's exact constants.

**⚠ Correction (later the same day, after the spreadsheet arrived):** the
.mdl comment values (737 / 1701 / 786 / 45.03 / 24.69) are the **Test
Data** column of BRU-ModelOutput-ParameterBookKeeping.xlsx — i.e. NASA
TN D-5815 design/test values, *not* NPSS outputs. The table above therefore
shows the replica matching the test/design data (which is self-consistent
with isentropic 0.80/0.87, no bleed, PR ≈ 1.75 — hence the clean sweep),
and conclusions 1–2 describe that data's provenance rather than an NPSS
run. The actual NPSS output is analyzed in the next section.

## 2026-06-10 — Excel received: first true NPSS output column (rung 2 proper)

`reference/BRU-ModelOutput-ParameterBookKeeping.xlsx` holds a full
station-by-station **Test Data vs NPSS** comparison (with the collaborator's
diff formulas flagging stations 4, 6/7, 12) plus a parasitic-loss table
(236 + 650 + 179 + 95 + 665 = 1825 W — kW-scale quantities, consistent
with HPX-as-kW) and an oil-Cp bookkeeping note.

Forensics on the NPSS column (`validation/bru3_excel_run.jl`):

- TIT pinned at 2060 °R (no shaft balance in play), bleed flowing
  (0.03 lb/s = 2.27 % at 540 °R — not heated to 559 °R), **heater ΔP = 0**
  (st 4→5 pressure unchanged ⇒ turbine PR = 1.8285, not 1.75), sink
  gas-side dPqP = 0.005 (not the .mdl's 0.0173), and the HotStart tear
  destroys 0.04 lb/s (st6 W = 1.32 → st7 W = 1.28).

Element-by-element replication at the exact NPSS port states:

| element | GasCycle (best match) | NPSS | Δ |
|---|---|---|---|
| Turbine | **polytropic η = 0.87**: 1642.04 °R | 1642.07 | **−0.03 °R** |
| HeatSinkHx gas out | ε-NTU, C_min = gas: 540.08 °R | 540.00 | +0.08 °R |
| Compressor | isen 0.784 → 741.61; poly 0.809 → 741.33 | 741.70 | η ≠ 0.80 under either semantics |
| Recup cold out | ε = 0.95: 1590.40 °R | 1580.02 | implied ε ≈ 0.938–0.941 |

So: the NPSS turbine is **polytropic** (matched at solver precision); the
test/design column is what's isentropic-consistent. The earlier README
inference ("NPSS effDes is isentropic") came from matching the 737 anchor —
which is test data. Open items now have sharp, single-number questions:

1. Compressor: what η/semantics produced 741.70 °R? (0.784 isen / 0.809
   poly; the .mdl says effDes = 0.80 — possibly 0.81 polytropic in that run)
2. Recuperator: what is NPSS HeatExchanger's `effect` definition? The run
   behaves as C_min-based ε ≈ 0.941 despite `effect = 0.95`.
3. Oil cp actually used: implied 0.674 BTU/(lbm·R) = 2.82 kJ/(kg·K); the
   sheet lists 1.8 ("current"), 0.884 ("correct" Dow 200), and Oil.fpt has
   0.8 — none match. (Gas loop unaffected: oil inlet Tt floated to close.)
4. HPX units: still open — this run pinned TIT, so the shaft balance never
   acted. The loss table supports kW.

The spreadsheet's flagged Test-vs-NPSS discrepancies fully decompose:
station 6/7 Δ = 58.9 °R = 24.9 (PR 1.8285 from missing heater ΔP)
+ 28.1 (bleed dilution) + 5.5 (residual semantics/η); station 4's 73 °R is
the recuperator amplifying station 7 (0.95 × 59 ≈ 56) plus comp-outlet and
flow-bookkeeping effects; station 12's 14.4 °R is the oil-cp issue.
**In other words: the collaborator's flagged discrepancies are not
mysteries — they are the .mdl's own known issues (bleed port, missing
heater ΔP, oil cp), now quantified individually.**

## 2026-06-10 — HeXe.out: SOLVER-PRECISION MATCH ACHIEVED (campaign goal met)

Daniel got NPSS 3.3 running and produced `reference/HeXe.out`: BRU3.mdl
converged (4 iterations) with **the same HeXe84.fpt table GasCycle reads**.
Replication: `validation/bru3_hexe_run.jl`.

**Rung 0 — exact.** ht and s match every printed digit at every gas
station; γ matches all five printed digits (1.66743 / 1.66880 / 1.66778 /
1.66746). Table parsing, unit constants, and reference states are
identical.

**NPSS's FPT interpolation is linear-in-P — confirmed to 0.002 °R.** NPSS
station 1 prints 751.47 °R: exactly GasCycle's `s_interp=:linear` value
(751.468), i.e. NPSS carries the same mid-cell entropy artifact we found
and fixed at rung 0 (+14 °R at this state vs the physical 737.16).
GasCycle's default `:log_pressure` mode is the *more correct* answer;
`:linear` exists precisely to reproduce NPSS.

**Run configuration recovered from the listing:** machines isentropic
(eff = 0.8000/0.8700 printed exactly; efPoly 0.8176/0.8452 derived);
turbine PR floated to 1.757 = our predicted pressure-closure value; TIT
floated to 2024.04 °R via shaft balance; heater ΔP present; sink gas
dPqP = 0.0173; tear mass-consistent (W = 1.32 both sides — fixed since the
Excel run). **HPX resolved:** displays 18.00 hp = 13.42 kW — the .mdl's
expression is kW-tagged and NPSS stores hp internally. **Oil cp resolved:**
the listing's oil-side energy balance implies 0.7999 BTU/(lbm·R) and
ht = 0.8·T exactly — Oil.fpt's constants. **Recup `effect` resolved:**
C_min-based ε (cold-side prediction exact to 0.08 °R at ε = 0.95).

**Bleed bookkeeping decoded (dhb/dh = 0, dPb/dP = 1):** NPSS's interstage
bleed does no work in either machine — comp pwr 17.1 = main-flow-only Δh
(full-flow is 17.40), and the turbine expands the main flow alone with the
559 °R bleed diluting at the *exit* (main-only 1646.4 °R → mixed 1624.6 ≈
1624.41). Physically the bleed is granted free recompression to 45 psia;
GasCycle replicates it as a split-before-compressor bypass with a free
repressurization, mixed after the turbine.

**Final comparison (both codes solving their own shaft balance, TIT an
output of each):**

| quantity | GasCycle (:linear) | NPSS | Δ |
|---|---|---|---|
| TIT (solved) | 2024.015 °R | 2024.04 | −0.025 |
| comp outlet | 751.468 °R | 751.47 | −0.002 |
| turbine outlet (mixed) | 1624.395 °R | 1624.41 | −0.015 |
| recup cold out | 1580.691 °R | 1580.77 | −0.079 |
| recup hot out | 810.789 °R | 810.74 | +0.049 |
| turbine PR / exit Pt | 1.757 / 24.660 | 1.757 / 24.660 | exact |
| comp / turb power | 17.051 / 30.471 kW | "17.1" / "30.5" | within print rounding |
| net = HPX | 13.420 kW | 13.42 | exact |

**That is the campaign's success criterion: agreement within solver/print
precision, with every contributing modeling choice named.** The one
remaining open item: the sink HX behaves as ε = 0.9405 vs its 0.946
setting (the recuperator honors its setting exactly) — it does not affect
the gas loop (the oil inlet temperature is NPSS's closure knob and floats
to absorb it), but it is worth one collaborator question for completeness.

**Physical-mode footnote:** with GasCycle's default `:log_pressure`
entropy interpolation (the artifact-free answer on the same table), the
same model solves to TIT = 2083 °R, comp outlet 737.2 °R, net 13.42 kW —
i.e. NPSS's converged TIT is ~59 °R below the physically-consistent value
because of its s-table interpolation. For mission work (SR-1 off-design,
inventory, transients) GasCycle's default mode is the one to trust;
`:linear` is for NPSS cross-checks.
