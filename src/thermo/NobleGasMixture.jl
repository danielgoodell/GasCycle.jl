"""
Direct property backend for noble gases and their binary mixtures
(roadmap item 4), replacing FPT tables for He-Xe (and Ar/Kr/Ne stand-ins).

Thermodynamics from Tournier, El-Genk & Gallo, AIAA 2006-4154 ("Best
Estimates of Binary Gas Mixtures Properties for Closed Brayton Cycle Space
Applications"): virial EOS

    P = Rg·T·(ρ̂ + B·ρ̂² + C·ρ̂³)                                  (Eq. 8)

with corresponding-states correlations for the reduced virial coefficients
(θ = T/Tcr; equation images verified against Boyle-point physics and the
paper's figures — the published text extraction garbles the exponents):

    B/V*  = −102.6 + (102.732 − 0.001θ − 0.44/θ^1.22)·tanh(4.5√θ)   (Eq. 9)
    B_He  = 8.4 − 0.0018T + 115/√T − 835/T   [cm³/mol]             (Eq. 10)
    C/V*² = 0.0757 + (−0.0862 − 3.6e−5·θ + 0.0237/θ^0.059)·tanh(0.84θ) (Eq. 11)

Enthalpy/Cp/Cv from B(T), C(T) derivatives (Eqs. 12–15); mixture B, C from
Eqs. 16–17 with combining rules Eqs. 20–21 for heavy pairs.  The paper fits
He-pair interaction coefficients B₁₂ to data without publishing them, so
He pairs use the classical Lennard-Jones second-virial series with the LJ
parameters from Johnson (NASA/CR-2006-214394) and standard combining rules
(σ₁₂ = (σ₁+σ₂)/2, ε₁₂ = √(ε₁ε₂)) — this reproduces the paper's Figure 6
He-Xe data.  Third virial cross-terms use Eq. 17b with C_He = 0 (paper:
"for He, 3rd virial coefficient can be neglected").

Entropy (not given explicitly in the paper) from the Helmholtz departure
of the same EOS:

    ŝ(T,P) = ŝ°(T,P) + Rg·ln Z − Rg·[(B + T·B′)ρ̂ + (C + T·C′)ρ̂²/2]

All functions are closed-form/AD-generic; ForwardDiff Duals propagate
through every method, including through the mole fraction x₁ (the
constructor is generic), enabling d(cycle)/d(mixture ratio) optimization.

Transport (μ, k, Pr) from the same paper's Eqs. 2–7 and 23–36 with
Tables 1–3 constants; Johnson's tables serve as validation oracles.
"""

using ForwardDiff: ForwardDiff, derivative

# ── Per-gas constants (El-Genk Table 1; LJ parameters per Johnson) ───────────
struct NobleGas
    name::String
    M::Float64     # kg/mol
    Tcr::Float64   # K
    Pcr::Float64   # Pa
    ρcr::Float64   # kg/m³
    σ::Float64     # LJ collision diameter [m]
    εk::Float64    # LJ well depth ε/k [K]
    Aμ::Float64    # dilute viscosity μ° = Aμ(T−Tμ)^n  [Pa·s] (Eq. 3)
    Tμ::Float64    # [K]
    nμ::Float64
    Δμcr::Float64  # excess viscosity at the critical point [Pa·s] (NaN: none)
    λcr::Float64   # critical-point thermal conductivity [W/(m·K)]
end

const _Rg = 8.31441   # J/(mol·K), the paper's value

const HELIUM  = NobleGas("He", 0.004003,   5.2, 0.2275e6,  69.64, 2.576e-10,  10.22,
                         3.0629e-7, -21.33, 0.7243,   NaN,     34.32e-3)
const NEON    = NobleGas("Ne", 0.020179,  44.5, 2.678e6,  481.9,  2.789e-10,  35.7,
                         8.4528e-7,  16.47, 0.642584, 8.9e-6,  35.47e-3)
const ARGON   = NobleGas("Ar", 0.039948, 150.7, 4.863e6,  535.6,  3.418e-10, 124.0,
                         6.9891e-7,  65.70, 0.63977,  16.0e-6, 28.223e-3)
const KRYPTON = NobleGas("Kr", 0.0838,   209.5, 5.51e6,   908.4,  3.610e-10, 190.0,
                         6.9629e-7,  71.07, 0.667,    23.3e-6, 19.828e-3)
const XENON   = NobleGas("Xe", 0.13129,  289.7, 5.84e6,  1110.0,  4.055e-10, 229.0,
                         7.5683e-7, 112.31, 0.655473, 29.7e-6, 15.966e-3)

# ── Pair interaction coefficients for transport (Tables 2-3) ─────────────────
# μ₁₂(T) = Aμ·(T − Tμ)^n  (Eq. 32) and the λ₁₂ correction factor f₁₂ (Eq. 36).
struct PairTransport
    Aμ::Float64    # [Pa·s/K^n]
    Tμ::Float64    # [K]
    n::Float64
    f::Float64     # f₁₂ of Eq. 36
end

const _PAIR_TRANSPORT = Dict(
    ("He", "Ne") => PairTransport(8.81837e-7,  63.63, 0.614098, 0.895),
    ("He", "Ar") => PairTransport(6.57562e-7,  51.87, 0.602712, 0.940),
    ("He", "Kr") => PairTransport(11.2472e-7, 121.27, 0.508158, 1.020),
    ("He", "Xe") => PairTransport(3.40998e-7,  45.89, 0.658754, 1.060),
    ("Ne", "Ar") => PairTransport(10.4450e-7,  47.66, 0.620956, 0.845),
    ("Ne", "Kr") => PairTransport(15.2944e-7,  94.86, 0.565158, 0.900),
    ("Ne", "Xe") => PairTransport(18.2681e-7, 125.0,  0.522034, 0.940),
    ("Ar", "Kr") => PairTransport(18.5326e-7, 148.93, 0.541593, 0.850),
    ("Ar", "Xe") => PairTransport(16.9684e-7, 151.58, 0.542416, 0.935),
    ("Kr", "Xe") => PairTransport(11.2303e-7, 116.71, 0.631571, 0.871),
)

function _pair_transport(g1::NobleGas, g2::NobleGas)
    g1 === g2 && return PairTransport(g1.Aμ, g1.Tμ, g1.nμ, 1.0)
    pt = get(_PAIR_TRANSPORT, (g1.name, g2.name), nothing)
    isnothing(pt) && (pt = get(_PAIR_TRANSPORT, (g2.name, g1.name), nothing))
    isnothing(pt) && error("no transport pair coefficients for $(g1.name)-$(g2.name)")
    pt
end

_Vstar(g::NobleGas) = _Rg * g.Tcr / g.Pcr   # characteristic volume [m³/mol]

# ── Pure-gas virial coefficients [m³/mol], [m⁶/mol²] ─────────────────────────
# Each returns (value, d/dT, d²/dT²).  Derivatives are hand-derived so a
# property call costs a single pass (the speed gate vs FPT lookups); all
# primitives are smooth, so outer ForwardDiff still flows through.

"""ΨB(θ) of Eq. 9 with first and second θ-derivatives."""
function _ΨB2(θ)
    s  = sqrt(θ)
    τ  = tanh(4.5 * s)
    σ2 = 1 - τ^2                       # sech²
    p  = 0.44 / θ^2.22                 # 0.44·θ^−2.22
    g  = 102.732 - 0.001 * θ - p * θ   # 102.732 − 0.001θ − 0.44/θ^1.22
    g′ = -0.001 + 1.22 * p
    g″ = -2.7084 * p / θ               # −1.22·2.22·0.44·θ^−3.22
    τ′ = σ2 * 2.25 / s
    τ″ = σ2 * (-2 * τ * 5.0625 / θ - 1.125 / (s * θ))
    (-102.6 + g * τ, g′ * τ + g * τ′, g″ * τ + 2 * g′ * τ′ + g * τ″)
end

"""ΨC(θ) of Eq. 11 with first and second θ-derivatives."""
function _ΨC2(θ)
    w  = tanh(0.84 * θ)
    σ2 = 1 - w^2
    p  = 0.0237 / θ^1.059              # 0.0237·θ^−1.059
    q  = -0.0862 - 3.6e-5 * θ + p * θ
    q′ = -3.6e-5 - 0.059 * p
    q″ = 0.062481 * p / θ              # 0.059·1.059·0.0237·θ^−2.059
    w′ = 0.84 * σ2
    w″ = -2 * 0.84^2 * w * σ2
    (0.0757 + q * w, q′ * w + q * w′, q″ * w + 2 * q′ * w′ + q * w″)
end

function _B_pure2(g::NobleGas, T)
    if g === HELIUM                                              # Eq. 10
        is = 1 / sqrt(T)
        B  = (8.4 - 0.0018 * T + 115.0 * is - 835.0 / T) * 1e-6
        B′ = (-0.0018 - 57.5 * is / T + 835.0 / T^2) * 1e-6
        B″ = (86.25 * is / T^2 - 1670.0 / T^3) * 1e-6
        (B, B′, B″)
    else                                                         # Eqs. 9, 18
        V = _Vstar(g)
        Ψ, Ψ′, Ψ″ = _ΨB2(T / g.Tcr)
        (V * Ψ, V * Ψ′ / g.Tcr, V * Ψ″ / g.Tcr^2)
    end
end

function _C_pure2(g::NobleGas, T)
    g === HELIUM && return (zero(T), zero(T), zero(T))           # paper: neglect
    V2 = _Vstar(g)^2
    Ψ, Ψ′, Ψ″ = _ΨC2(T / g.Tcr)
    (V2 * Ψ, V2 * Ψ′ / g.Tcr, V2 * Ψ″ / g.Tcr^2)
end

_B_pure(g::NobleGas, T) = _B_pure2(g, T)[1]
_C_pure(g::NobleGas, T) = _C_pure2(g, T)[1]

# ── He-pair B₁₂: Lennard-Jones second-virial series ──────────────────────────
# B*(T*) = Σ b_j·T*^(−(2j+1)/4); b_j = −2^(j+1/2)/(4·j!)·Γ((2j−1)/4).
# Verified against Hirschfelder tables: B*(1) = −2.538, B*(10) = +0.4609.
const _LJ_B_COEFF = (
    +1.733000920185e+00, -2.563693352041e+00, -8.665004600924e-01,
    -4.272822253401e-01, -2.166251150231e-01, -1.068205563350e-01,
    -5.054586017206e-02, -2.289011921465e-02, -9.928651105225e-03,
    -4.132938191534e-03, -1.654775184204e-03, -6.387268114189e-04,
    -2.381873371203e-04, -8.598245538331e-05, -3.010059754817e-05,
    -1.023600659325e-05, -3.386317224169e-06, -1.091338938251e-06)

"""B*(T*) with T*-derivatives, using a running power (2 pows total):
the exponents −(2j+1)/4 step by −1/2 each term."""
function _Bstar_LJ2(Ts)
    t    = Ts^(-1 / 4)         # current term power
    step = 1 / sqrt(Ts)        # Ts^(−1/2)
    s = s′ = s″ = zero(Ts)
    e = -1 / 4
    for b in _LJ_B_COEFF
        s  += b * t
        s′ += b * e * t / Ts
        s″ += b * e * (e - 1) * t / Ts^2
        t  *= step
        e  -= 1 / 2
    end
    (s, s′, s″)
end

const _NA = 6.022045e23   # paper's Avogadro number

function _B12_2(g1::NobleGas, g2::NobleGas, T)
    if g1 === HELIUM || g2 === HELIUM
        σ12 = 0.5 * (g1.σ + g2.σ)
        ε12 = sqrt(g1.εk * g2.εk)
        b0  = (2π / 3) * _NA * σ12^3            # m³/mol
        s, s′, s″ = _Bstar_LJ2(T / ε12)
        (b0 * s, b0 * s′ / ε12, b0 * s″ / ε12^2)
    else
        # Prausnitz combining rules (Eqs. 19-21), validated for heavy pairs
        V1, V2 = _Vstar(g1), _Vstar(g2)
        V12 = 0.5 * (V1 + V2)                                    # Eq. 20
        β   = V1 / V2
        T12 = 4 * sqrt(β) / (1 + β)^2 * sqrt(g1.Tcr * g2.Tcr)    # Eq. 21
        Ψ, Ψ′, Ψ″ = _ΨB2(T / T12)                                # Eq. 19
        (V12 * Ψ, V12 * Ψ′ / T12, V12 * Ψ″ / T12^2)
    end
end

_B12(g1::NobleGas, g2::NobleGas, T) = _B12_2(g1, g2, T)[1]

# ── Mixture backend ───────────────────────────────────────────────────────────
"""
    NobleGasMixture(gas1, gas2, x1; name) <: FluidProperties

Binary noble-gas mixture with mole fraction `x1` of `gas1`.  `x1` may be a
ForwardDiff Dual for mixture-ratio design optimization.  Pure gases:
`NobleGasMixture(XENON, HELIUM, 1.0)` or the `NobleGasFluid(gas)` helper.

    HeXe(M_molar)  — He-Xe mixture specified by molecular weight [kg/kmol]
"""
struct NobleGasMixture{X<:Real} <: FluidProperties
    gas1::NobleGas
    gas2::NobleGas
    x1::X
    M::X            # mixture molar mass [kg/mol]
    name::String
    pair::PairTransport   # resolved once: Tables 2-3 lookup for transport
end

function NobleGasMixture(gas1::NobleGas, gas2::NobleGas, x1::Real;
                         name::String = "$(gas1.name)$(gas2.name)")
    0 <= ForwardDiff.value(x1) <= 1 ||
        error("NobleGasMixture: x1 must be in [0,1], got $x1")
    M = x1 * gas1.M + (1 - x1) * gas2.M                          # Eq. 22
    NobleGasMixture(gas1, gas2, promote(x1, M)..., name, _pair_transport(gas1, gas2))
end

NobleGasFluid(gas::NobleGas) = NobleGasMixture(gas, gas, 1.0; name = gas.name)

"""He-Xe mixture by molecular weight in kg/kmol (e.g. HeXe(83.8))."""
function HeXe(M_molar::Real)
    M = M_molar * 1e-3
    x_He = (XENON.M - M) / (XENON.M - HELIUM.M)
    # name from the primal value only, so M_molar may be a ForwardDiff Dual
    name = "HeXe$(round(ForwardDiff.value(M_molar), digits = 1))"
    NobleGasMixture(HELIUM, XENON, x_He; name)
end

"""
    _virial(fp, T) -> (B, B′, B″, C, C′, C″)

Mixture virial coefficients and T-derivatives in one pass (Eqs. 16-17).
The Eq. 17b geometric-mean C cross-terms reduce to scaled copies of the
heavy-gas ΨC when one component is He (C_He ≡ 0 would otherwise put a
cbrt-at-zero singularity in the chain rule).
"""
function _virial(fp::NobleGasMixture, T)
    if fp.gas1 === fp.gas2
        B, B′, B″ = _B_pure2(fp.gas1, T)
        C, C′, C″ = _C_pure2(fp.gas1, T)
        return (B, B′, B″, C, C′, C″)
    end
    x1 = fp.x1
    x2 = 1 - x1
    B1, B1′, B1″ = _B_pure2(fp.gas1, T)
    B2, B2′, B2″ = _B_pure2(fp.gas2, T)
    Bx, Bx′, Bx″ = _B12_2(fp.gas1, fp.gas2, T)
    w11, w12, w22 = x1^2, 2 * x1 * x2, x2^2
    B  = w11 * B1  + w12 * Bx  + w22 * B2
    B′ = w11 * B1′ + w12 * Bx′ + w22 * B2′
    B″ = w11 * B1″ + w12 * Bx″ + w22 * B2″

    # C cross-terms: Cijk = (Ci·Cj·Ck)^(1/3); each pure C = V*²·ΨC(θ), so
    # the geometric means scale the same ΨC evaluated at each gas's θ.
    # With one He component (C ≡ 0) the cross-terms vanish identically.
    C1, C1′, C1″ = _C_pure2(fp.gas1, T)
    C2, C2′, C2″ = _C_pure2(fp.gas2, T)
    if fp.gas1 === HELIUM
        w = x2^3
        C, C′, C″ = w * C2, w * C2′, w * C2″
    elseif fp.gas2 === HELIUM
        w = x1^3
        C, C′, C″ = w * C1, w * C1′, w * C1″
    else
        # both gases use ΨC at their own θ; cube-root mean per Eq. 17b
        r1, r2 = cbrt(C1), cbrt(C2)
        C112 = r1^2 * r2
        C122 = r1 * r2^2
        # d(u^⅓) = u′/(3r²);  d²(u^⅓) = u″/(3r²) − 2u′²/(9r⁵)
        d1, d2 = C1′ / (3 * r1^2), C2′ / (3 * r2^2)
        e1 = C1″ / (3 * r1^2) - 2 * C1′^2 / (9 * r1^5)
        e2 = C2″ / (3 * r2^2) - 2 * C2′^2 / (9 * r2^5)
        C112′ = 2 * r1 * d1 * r2 + r1^2 * d2
        C122′ = d1 * r2^2 + 2 * r1 * r2 * d2
        C112″ = 2 * (d1^2 + r1 * e1) * r2 + 4 * r1 * d1 * d2 + r1^2 * e2
        C122″ = e1 * r2^2 + 4 * d1 * r2 * d2 + 2 * r1 * (d2^2 + r2 * e2)
        w111, w112, w122, w222 = x1^3, 3 * x1^2 * x2, 3 * x1 * x2^2, x2^3
        C  = w111 * C1  + w112 * C112  + w122 * C122  + w222 * C2
        C′ = w111 * C1′ + w112 * C112′ + w122 * C122′ + w222 * C2′
        C″ = w111 * C1″ + w112 * C112″ + w122 * C122″ + w222 * C2″
    end
    (B, B′, B″, C, C′, C″)
end

_B_mix(fp::NobleGasMixture, T) = _virial(fp, T)[1]
_C_mix(fp::NobleGasMixture, T) = _virial(fp, T)[4]

# Molar density from Newton on the virial EOS (ideal-gas start; the cubic
# is mildly nonlinear at CBC conditions — 4 iterations reach machine eps)
function _rhom_BC(B, C, T, P)
    ρ̂ = P / (_Rg * T)
    for _ in 1:4
        f  = _Rg * T * (ρ̂ + B * ρ̂^2 + C * ρ̂^3) - P
        f′ = _Rg * T * (1 + 2B * ρ̂ + 3C * ρ̂^2)
        ρ̂ -= f / f′
    end
    ρ̂
end

function _rhom(fp::NobleGasMixture, T, P)
    v = _virial(fp, T)
    _rhom_BC(v[1], v[4], T, P)
end

# ── FluidProperties interface (all analytic, Eqs. 12-15) ─────────────────────
density(fp::NobleGasMixture, T, P) = _rhom(fp, T, P) * fp.M

function enthalpy(fp::NobleGasMixture, T, P)
    B, B′, _, C, C′, _ = _virial(fp, T)
    ρ̂ = _rhom_BC(B, C, T, P)
    (2.5 * _Rg * T +
     ρ̂ * _Rg * T * ((B - T * B′) + ρ̂ * (C - T * C′ / 2))) / fp.M  # Eq. 12
end

function cp(fp::NobleGasMixture, T, P)
    B, B′, B″, C, C′, C″ = _virial(fp, T)
    ρ̂ = _rhom_BC(B, C, T, P)
    # Eq. 14: density change along the isobar
    dρ̂dT = -((ρ̂ + B * ρ̂^2 + C * ρ̂^3) / T + B′ * ρ̂^2 + C′ * ρ̂^3) /
            (1 + 2B * ρ̂ + 3C * ρ̂^2)
    # Eq. 13
    ĉp = 2.5 * _Rg +
         ρ̂ * _Rg * ((B - T * B′ - T^2 * B″) + ρ̂ * (C - T^2 * C″ / 2)) +
         _Rg * T * ((B - T * B′) + ρ̂ * (2C - T * C′)) * dρ̂dT
    ĉp / fp.M
end

function entropy(fp::NobleGasMixture, T, P)
    B, B′, _, C, C′, _ = _virial(fp, T)
    ρ̂ = _rhom_BC(B, C, T, P)
    Z = P / (ρ̂ * _Rg * T)
    (2.5 * _Rg * log(T / 298.15) - _Rg * log(P / 101325.0) + _Rg * log(Z) -
     _Rg * ((B + T * B′) * ρ̂ + (C + T * C′) * ρ̂^2 / 2)) / fp.M
end

function gamma(fp::NobleGasMixture, T, P)
    B, B′, B″, C, C′, C″ = _virial(fp, T)
    ρ̂ = _rhom_BC(B, C, T, P)
    dρ̂dT = -((ρ̂ + B * ρ̂^2 + C * ρ̂^3) / T + B′ * ρ̂^2 + C′ * ρ̂^3) /
            (1 + 2B * ρ̂ + 3C * ρ̂^2)
    ĉp = 2.5 * _Rg +
         ρ̂ * _Rg * ((B - T * B′ - T^2 * B″) + ρ̂ * (C - T^2 * C″ / 2)) +
         _Rg * T * ((B - T * B′) + ρ̂ * (2C - T * C′)) * dρ̂dT
    # Eq. 15
    ĉv = 1.5 * _Rg -
         ρ̂ * _Rg * T * ((2 * B′ + T * B″) + ρ̂ * (C′ + T * C″ / 2))
    ĉp / ĉv
end

# Newton inversions with exact derivatives (2-4 iterations; h is nearly
# linear in T for monatomic gases)
function T_from_h(fp::NobleGasMixture, h_target, P; T_guess = nothing)
    T = isnothing(T_guess) ? h_target * fp.M / (2.5 * _Rg) : float(T_guess)
    for _ in 1:6
        r = enthalpy(fp, T, P) - h_target
        T -= r / cp(fp, T, P)
        abs(r) < 1e-9 * abs(h_target) + 1e-12 && break
    end
    T
end

function T_from_s(fp::NobleGasMixture, s_target, P; T_guess = 500.0)
    cp0 = 2.5 * _Rg / fp.M
    # ideal-gas closed-form start, then Newton with ds/dT = cp/T
    T = 298.15 * exp((s_target + (_Rg / fp.M) * log(P / 101325.0)) / cp0)
    for _ in 1:6
        r = entropy(fp, T, P) - s_target
        T -= r * T / cp(fp, T, P)
        abs(r) < 1e-11 * cp0 && break
    end
    T
end

# ── Transport: μ, k, Pr (Eqs. 2-7, 23-36) ────────────────────────────────────
# Dilute baselines (Eqs. 3, 7) and the excess-property shape functions of
# reduced density (Fig. 2 fit; Eq. 4):
_μ0(g::NobleGas, T) = g.Aμ * (T - g.Tμ)^g.nμ
_λ0(g::NobleGas, T) = 3.75 * (_Rg / g.M) * _μ0(g, T)
_Ψμ(ρr) = ρr * (0.221 + ρr * (1.062 + ρr * (-0.509 + 0.225 * ρr)))
_Ψλ(ρr) = ρr * (0.645 + ρr * (0.33 + ρr * (0.0368 - 0.0128 * ρr)))

_μ12(pt::PairTransport, T) = pt.Aμ * (T - pt.Tμ)^pt.n            # Eq. 32

const _Astar = 1.10   # A*₁₂ and B*₁₂: near-constant over CBC temperatures
const _Bstar = 1.10   # (paper takes both as 1.10 throughout)

"""Dilute mixture viscosity, Sutherland-Wassiljewa (Eqs. 30-31)."""
function _μ0_mix(fp::NobleGasMixture, T)
    g1, g2 = fp.gas1, fp.gas2
    x1 = fp.x1
    x2 = 1 - x1
    μ1, μ2 = _μ0(g1, T), _μ0(g2, T)
    μ12 = _μ12(fp.pair, T)
    m1, m2 = g1.M, g2.M
    mm  = 2 * m1 * m2 / (m1 + m2)^2
    φ12 = (μ1 / μ12) * mm * (5 / (3 * _Astar) + m2 / m1)
    φ21 = (μ2 / μ12) * mm * (5 / (3 * _Astar) + m1 / m2)
    μ1 / (1 + φ12 * x2 / x1) + μ2 / (1 + φ21 * x1 / x2)
end

"""Dilute mixture conductivity, first-order Hirschfelder (Eqs. 34-36)."""
function _λ0_mix(fp::NobleGasMixture, T)
    g1, g2 = fp.gas1, fp.gas2
    x1 = fp.x1
    x2 = 1 - x1
    λ1, λ2 = _λ0(g1, T), _λ0(g2, T)
    m1, m2 = g1.M, g2.M
    m12 = 2 * m1 * m2 / (m1 + m2)                                # harmonic-ish mean
    λ12 = 3.75 * (_Rg / m12) * _μ12(fp.pair, T) * fp.pair.f      # Eq. 36
    s   = (m1 + m2)^2 * _Astar
    L11 = x1^2 / λ1 + x1 * x2 / (2λ12) *
          (7.5 * m1^2 + 6.25 * m2^2 - 3 * m2^2 * _Bstar + 4 * m1 * m2 * _Astar) / s
    L22 = x2^2 / λ2 + x1 * x2 / (2λ12) *
          (7.5 * m2^2 + 6.25 * m1^2 - 3 * m1^2 * _Bstar + 4 * m1 * m2 * _Astar) / s
    L12 = -x1 * x2 / (2λ12) * (m1 * m2 / s) * (55 / 4 - 3 * _Bstar - 4 * _Astar)
    (x1^2 * L22 - 2 * x1 * x2 * L12 + x2^2 * L11) / (L11 * L22 - L12^2)
end

# Pseudo-critical mixing for the dense corrections.  ρr = 0.291·V*·ρ̂: with
# Zc = 0.291, 0.291·V* is the (pseudo-)critical molar volume, so the Ψ
# argument is ρ/ρcr.  Viscosity T* uses the linear rule (Eq. 24b); the
# conductivity T* uses the van der Waals double sum (Eq. 26) with the
# Eq. 20-21 cross terms.
function _Vstar_mix(fp::NobleGasMixture)
    fp.x1 * _Vstar(fp.gas1) + (1 - fp.x1) * _Vstar(fp.gas2)      # Eq. 24a/25
end

function _mustar(M, Tstar, Vstar)                                # Eq. 23b
    0.204e-7 * sqrt(M * Tstar) / (0.291 * Vstar)^(2 / 3)
end

function _lamstar(M, Tstar, Vstar)                               # Eq. 33b
    # SI throughout (M kg/mol, V* m³/mol): reproduces the Table 1 λcr
    # column to its printed "λ*cr deviation" row for all five gases.
    0.304e-4 * Tstar^0.277 / (M^0.465 * (0.291 * Vstar)^0.415)
end

"""Pure-gas viscosity (Eq. 2).  He's excess viscosity is nil (Table 1)."""
function _viscosity_pure(fp::NobleGasMixture, g::NobleGas, T, P)
    μ = _μ0(g, T)
    isnan(g.Δμcr) && return μ + zero(_rhom(fp, T, P))   # keep type uniform
    μ + g.Δμcr * _Ψμ(_rhom(fp, T, P) * g.M / g.ρcr)
end

"""Pure-gas conductivity (Eq. 6)."""
function _conductivity_pure(fp::NobleGasMixture, g::NobleGas, T, P)
    _λ0(g, T) + (1 - 1 / 2.94) * g.λcr * _Ψλ(_rhom(fp, T, P) * g.M / g.ρcr)
end

function viscosity(fp::NobleGasMixture, T, P)
    g1, g2 = fp.gas1, fp.gas2
    xv = ForwardDiff.value(fp.x1)
    (g1 === g2 || xv == 1) && return _viscosity_pure(fp, g1, T, P)
    xv == 0 && return _viscosity_pure(fp, g2, T, P)
    μ0 = _μ0_mix(fp, T)
    V  = _Vstar_mix(fp)
    ρr = 0.291 * V * _rhom(fp, T, P)
    if g1 === HELIUM || g2 === HELIUM
        # Eq. 28: He contributes no excess viscosity (Eq. 27), so the dense
        # correction is the heavy gas's μ*, weighted by its mole fraction.
        heavy, xh = g1 === HELIUM ? (g2, 1 - fp.x1) : (g1, fp.x1)
        μstar = _mustar(heavy.M, heavy.Tcr, _Vstar(heavy))
        return μ0 + 0.565 * xh * μstar * _Ψμ(ρr)
    end
    # Eq. 23 with the Eq. 24 mixing rules
    x1 = fp.x1
    VT = x1 * _Vstar(g1) * g1.Tcr + (1 - x1) * _Vstar(g2) * g2.Tcr   # Eq. 24b
    μ0 + 0.565 * _mustar(fp.M, VT / V, V) * _Ψμ(ρr)
end

function conductivity(fp::NobleGasMixture, T, P)
    g1, g2 = fp.gas1, fp.gas2
    xv = ForwardDiff.value(fp.x1)
    (g1 === g2 || xv == 1) && return _conductivity_pure(fp, g1, T, P)
    xv == 0 && return _conductivity_pure(fp, g2, T, P)
    λ0 = _λ0_mix(fp, T)
    x1 = fp.x1
    x2 = 1 - x1
    V  = _Vstar_mix(fp)
    # Eq. 26 (van der Waals second rule) with Eq. 20-21 cross terms
    V1, V2 = _Vstar(g1), _Vstar(g2)
    V12 = 0.5 * (V1 + V2)                                        # Eq. 20
    β   = V1 / V2
    T12 = 4 * sqrt(β) / (1 + β)^2 * sqrt(g1.Tcr * g2.Tcr)        # Eq. 21
    VT  = x1^2 * V1 * g1.Tcr + 2 * x1 * x2 * V12 * T12 + x2^2 * V2 * g2.Tcr
    ρr  = 0.291 * V * _rhom(fp, T, P)
    λ0 + (1 - 1 / 2.94) * _lamstar(fp.M, VT / V, V) * _Ψλ(ρr)    # Eq. 33
end

prandtl(fp::NobleGasMixture, T, P) =
    cp(fp, T, P) * viscosity(fp, T, P) / conductivity(fp, T, P)
