# Transport implementation notes (roadmap item 4, second half)

Extracted 2026-06-10 from El-Genk AIAA 2006-4154 (implementation source)
and Johnson NASA/CR-2006-214394 (validation oracle). The thermo half is
done and committed (`src/thermo/NobleGasMixture.jl`, 6a13300). This file
preserves everything needed to implement μ, k, Pr without re-reading the
PDFs.

## Equations (El-Genk)

Pure dilute:
- Eq 3:  μ°(T) = Aμ·(T − Tμ)^n        — Aμ, Tμ, n per gas already in the
  `NobleGas` struct (Table 1), SI (Pa·s).
- Eq 7:  λ°(T) = (15/4)·(Rg/M)·μ°(T)

Pure dense (ρr = ρ/ρcr):
- Ψμ(ρr) = 0.221ρr + 1.062ρr² − 0.509ρr³ + 0.225ρr⁴   (Fig. 2 fit)
- Eq 4: Ψλ(ρr) = 0.645ρr + 0.33ρr² + 0.0368ρr³ − 0.0128ρr⁴
- Eq 2:  μ = μ°(T) + Δμcr·Ψμ(ρr)   — Δμcr per gas in struct (He: none;
  helium viscosity is pressure-independent). [(1−1/2.3)·μcr ≡ Δμcr]
- Eq 6:  λ = λ°(T) + (1 − 1/2.94)·λcr·Ψλ(ρr) — λcr per gas from Table 1:
  He 34.32e−3, Ne 35.47e−3, Ar 28.223e−3, Kr 19.828e−3, Xe 15.966e−3 W/m·K
  (NOT yet in the struct — add a λcr field, or resolve Eq 5's units).

Dilute mixture:
- Eq 30/31 (Sutherland-Wassiljewa) for μ°mix with A*₁₂ = 1.10 (constant):
  μ° = μ°₁/(1 + φ₁₂x₂/x₁) + μ°₂/(1 + φ₂₁x₁/x₂)
  φ₁₂ = (μ°₁/μ₁₂)·[2m₁m₂/(m₁+m₂)²]·[5/(3A*₁₂) + m₂/m₁]
  φ₂₁ = (μ°₂/μ₁₂)·[2m₁m₂/(m₁+m₂)²]·[5/(3A*₁₂) + m₁/m₂]
- Eq 32: μ₁₂(T) = Aμ(T − Tμ)^n with the PAIR coefficients (Tables 2-3 below).
- Eq 34/35 full first-order λ°mix with A*₁₂ = B*₁₂ = 1.10:
  λ° = [x₁²/L₁₁ + 2x₁x₂·(−L₁₂... ) ...]  — standard determinant form:
  λ° = (x₁²/L₁₁ + 2x₁x₂L₁₂/(L₁₁L₂₂) + x₂²/L₂₂)... use the paper's exact form:
  λ° = [x₁²/L₁₁ − 2x₁x₂L₁₂/(L₁₁L₂₂) + x₂²/L₂₂] × [1 − L₁₂²/(L₁₁L₂₂)]⁻¹
  L₁₁ = x₁²/λ°₁ + (x₁x₂/2λ₁₂)·[(15/2)m₁² + (25/4)m₂² − 3m₂²B*₁₂ + 4m₁m₂A*₁₂]/[(m₁+m₂)²A*₁₂]
  L₂₂ = same with 1↔2
  L₁₂ = −(x₁x₂/2λ₁₂)·[m₁m₂/((m₁+m₂)²A*₁₂)]·(55/4 − 3B*₁₂ − 4A*₁₂)
  (mind the sign convention; verify against the dilute Chapman-Enskog
  values in Johnson's tables)
- Eq 36: λ₁₂ = (15k/4m₁₂)·μ₁₂·f₁₂,  m₁₂ = 2m₁m₂/(m₁+m₂); f₁₂ in Table 3.

Dense mixture (argument of Ψ is 0.291·V*·ρ̂ — reduced density via the
pseudo-critical volume; ρ̂ = molar density from the thermo half):
- Eq 23a: μ = μ°mix(T) + 0.565·μ*·Ψμ(0.291V*ρ̂)
- Eq 23b: μ* = 0.204e−7·√(M·T*)/(0.291·V*)^(2/3)
  UNITS VERIFIED: M kg/mol, V* m³/mol → Pa·s (reproduces Table 1 μcr,
  e.g. Xe 52.25e−6 to −0.08%).
- viscosity mixing rules: V* = x₁V*₁₁ + x₂V*₂₂ (Eq 24a→25);
  V*T* = Σᵢ xᵢV*ᵢᵢT*ᵢᵢ (Eq 24b).
- He-mixture special case Eq 28 (He defect viscosity is nil, Eq 27):
  μ = μ°mix(T) + 0.565·x_heavy·μ*_heavy,pure·Ψμ(0.291V*ρ̂)
- Eq 33a: λ = λ°mix(T) + (1 − 1/2.94)·λ*·Ψλ(0.291V*ρ̂)
- Eq 33b: λ* = 0.304e−4·(T*)^0.277/(M^0.465·(0.291V*)^0.415)
  ⚠ UNITS UNRESOLVED: plugging Xe pure values does NOT reproduce Table 1
  λcr = 15.966e−3 W/m·K in any obvious unit combination tried (kg vs g,
  m³ vs cm³). Resolve by calibrating against the Table 1 λcr column and
  its "λ*cr deviation %" row (He +6.5%, Ne +0.47%, Ar −0.50%, Kr +0.99%,
  Xe −0.28%) — those deviations define the intended formula output.
- conductivity mixing rules: V* per Eq 24a(=25); V*T* = ΣᵢΣⱼ xᵢxⱼV*ᵢⱼT*ᵢⱼ
  (Eq 26, van der Waals 2nd rule) with V*ᵢⱼ, T*ᵢⱼ from Eqs 20-21.

## Table 2 — pair coefficients for μ₁₂ (Eq 32)

Aμ×10⁷ (Pa·s/K^n):
  He-Ne 8.81837, He-Ar 6.57562, He-Kr 11.2472, He-Xe 3.40998,
  Ne-Ar 10.4450, Ne-Kr 15.2944, Ne-Xe 18.2681,
  Ar-Kr 18.5326, Ar-Xe 16.9684, Kr-Xe 11.2303
Tμ (K):
  He-Ne 63.63, He-Ar 51.87, He-Kr 121.27, He-Xe 45.89,
  Ne-Ar 47.66, Ne-Kr 94.86, Ne-Xe 125.0,
  Ar-Kr 148.93, Ar-Xe 151.58, Kr-Xe 116.71

## Table 3 — exponent n (Eq 32) and correction factor f₁₂ (Eq 36)

n:
  He-Ne 0.614098, He-Ar 0.602712, He-Kr 0.508158, He-Xe 0.658754,
  Ne-Ar 0.620956, Ne-Kr 0.565158, Ne-Xe 0.522034,
  Ar-Kr 0.541593, Ar-Xe 0.542416, Kr-Xe 0.631571
f₁₂:
  He-Ne 0.895, He-Ar 0.940, He-Kr 1.020, He-Xe 1.060,
  Ne-Ar 0.845, Ne-Kr 0.900, Ne-Xe 0.940,
  Ar-Kr 0.850, Ar-Xe 0.935, Kr-Xe 0.871

## Validation oracles

- Johnson NASA/CR-2006-214394 Tables 4-6: third-order μ and k for He-Xe at
  M = 20.183, 39.94, 83.8 g/mol, T = 400-1200 K — extract from the PDF
  (text dump worked: /tmp/johnson.txt; tables start ~line 378 for Ω
  tables; data tables 4-6 later). Dilute-limit oracle (0.1 MPa).
  Johnson Table 7: Prandtl vs Taylor-1988 experiment.
- El-Genk §V spot values: at 300 K / 2 MPa — Xe: λ +12.3%, μ +4.3% vs
  dilute; He-Xe M=40: λ ~+0.8%, μ ~+0.5%. At 400 K / 2 MPa the M=40
  mixture errors are < 1%.
- Sanity: Pr of monatomic dilute gas = cp·μ/λ = (5/2 R)·μ/((15/4)Rμ) = 2/3
  exactly for a PURE gas (Eq 7); mixtures deviate (that's their benefit —
  Pr ~ 0.2-0.3 for He-Xe).

## Implementation plan

1. Add `viscosity(fp, T, P)`, `conductivity(fp, T, P)`, `prandtl(fp, T, P)`
   to the FluidProperties interface (error fallbacks in base, like cp).
2. Add λcr to `NobleGas` (Table 1 values above) — sidesteps the Eq 33b
   unit question for pure gases; resolve 33b only for mixtures (calibrate
   so pure-gas endpoints reproduce Table 1 within the printed deviations).
3. Pair-coefficient lookup keyed on (gas1, gas2) for Tables 2-3.
4. Tests: Johnson Tables 4-6 at M = 83.8 (and 40 if present), Pr → 2/3
   pure-gas dilute limit, He-Xe Pr ≈ 0.2-0.3 band, El-Genk §V spot values,
   AD-through-T smoothness, speed (transport is HX-sizing frequency, so
   the gate is loose: < 5 μs is fine).
5. Transport for other backends: optional (FPT tables lack μ/k; H2O/Oil
   function files return NaN/constants — leave unimplemented).
"""
