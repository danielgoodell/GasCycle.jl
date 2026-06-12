# data/ — fluid property table files

NPSS-format FPT files consumed by the tests, examples, benchmarks, and
validation scripts (`FPTFluid` reads the table files; `ConstantPropertyLiquid`
reads the function-style ones). They were moved here from the repository root
on 2026-06-11; all code references point at `data/`.

## HeXe84.fpt

He-Xe binary mixture property table, NPSS ENGLISH units (°R, psia, BTU/lbm).
Generated in 2026 by **Joel Krakower (NASA)** with a script implementing the
same property methods as this package's `NobleGasMixture` backend:

- Thermodynamics (virial EOS): J.-M. Tournier, M. S. El-Genk and B. M. Gallo,
  *"Best Estimates of Binary Gas Mixtures Properties for Closed Brayton Cycle
  Space Applications,"* AIAA 2006-4154.
- Transport (μ, k): P. K. Johnson, *"A Method for Calculating Viscosity and
  Thermal Conductivity of a Helium-Xenon Gas Mixture,"* NASA/CR-2006-214394.

Both papers are in `reference/`.

Grid: 32 Pt nodes (0.725–435.1 psia ≈ 5 kPa–3.0 MPa) × 32 Tt nodes
(450–2520 °R = 250–1400 K). Thirteen tables: forward `h_T`, `s_T`, `Cp`, `Cv`,
`R`, `gam`, `rho`, `mu`, `k`, `Pr` and inverse `T_h`, `T_s`, `h_s`. The
`FPTFluid` reader uses the thermo and inverse tables; the transport tables
(`mu`, `k`, `Pr`) are present in the file but not currently read — transport
comes from `NobleGasMixture`.

### Audit vs `NobleGasMixture` (2026-06-11)

- **Composition: the table is M ≈ 84.07 kg/kmol, not the BRU-spec 83.8.**
  The `R`, `Cp` and `Cv` tables at the ideal-gas corner all imply
  M = 84.06–84.08 self-consistently. Comparisons in `test/test_noblegas.jl`
  use `HeXe(83.8)` with ~1% tolerances, which absorbs the 0.3% offset.
- **`h_T`, `s_T`, `rho` match the El-Genk method.** Against `HeXe(84.07)`
  over the BRU cycle envelope (285–1340 K, 20–310 kPa): Δh within 0.05%,
  ρ within 0.35%, Δs within ~1% (interpolation-dominated).
- **`mu` and `k` match the Johnson/El-Genk method** to ≤0.06% and ≤0.9%
  respectively across the full grid — confirming the same-methods provenance.
- **The `Cp` table (and the derived `gam`, `Pr`) carries a wrong-signed
  real-gas departure.** At every probed state, `Cp_table ≈ Cp_ideal − 0.75·δ`
  where the virial method (and the file's own `h_T` table, via dh/dT) gives
  `Cp_ideal + δ`. E.g. at 300 K / 1 MPa: d(h_T)/dT = 258.2 J/(kg·K) (virial:
  259.1) but the Cp table reads 238.3. The error is ≤0.6% at BRU cycle
  conditions (≥300 K, ≤300 kPa) but reaches −42% at the 250 K / 3 MPa corner
  (near Xe condensation). GasCycle cycle solves are h/s-based and essentially
  unaffected (cp enters only ε-NTU capacity rates and reporting); flag this
  to the table's author before using `Cp`/`gam`/`Pr` columns at high
  pressure and low temperature.

## H2O.fpt

Function-style constant-property file: Cp = 1.0 BTU/(lbm·°R),
ρ = 62.37 lbm/ft³ — liquid water at ~60 °F. Note the `description` string
inside the file ("Methane tables from Justin Gray") is a leftover from the
template it was copied from; the values are water.

## Oil.fpt

Function-style constant-property file used as the heat-sink coolant by the
collaborator's `BRU3.mdl`/`HeatSinkHX.mdl` (NPSS composition name `"Oil"`):
Cp = 0.8 BTU/(lbm·°R) = 3.35 kJ/(kg·K), ρ = 62.424 lbm/ft³ = 999.9 kg/m³.

**These values do not match Dow Corning 200**, the silicone coolant used in
the NASA BRU heat-rejection loop, despite the file's intent: DC-200
(polydimethylsiloxane, any viscosity grade) has cp ≈ 0.37–0.50 BTU/(lbm·°R)
(1.5–2.1 kJ/(kg·K)) and specific gravity 0.87–0.97 — Oil.fpt's cp is ~2×
too high and its density is exactly water's. The collaborator's bookkeeping
spreadsheet listed yet another value (0.884 BTU/(lbm·°R), labeled
"correct" Dow 200) that matches neither; see `validation/RESULTS.md`
(2026-06-10 entries). The NPSS `HeXe.out` run verifiably used 0.8
(its oil-side energy balance gives ht = 0.8·T exactly), so this file is kept
byte-identical — and under its NPSS name — for replication. Do not "fix" it;
a physically-correct DC-200 coolant should be a new file or a
`ConstantPropertyLiquid(cp=..., rho=...)` constructed directly.
