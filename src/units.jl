"""
Unit-conversion helpers.

NPSS models, FPT property tables, and the NASA reference reports all use
English engineering units (°R, psia, lbm/s, BTU); GasCycle works in SI.
These helpers centralize the conversions so examples and scripts don't each
redefine them.

All factors are exact by definition (International Table BTU, international
pound, international foot) except where noted.
"""

# ── Temperature ───────────────────────────────────────────────────────────────
R_to_K(T_R) = T_R * (5 / 9)
K_to_R(T_K) = T_K * (9 / 5)

# ── Pressure ──────────────────────────────────────────────────────────────────
const _PA_PER_PSI = 6.894757293168361e3   # 4.4482216152605 N / 0.0254² m²

psia_to_Pa(P_psia) = P_psia * _PA_PER_PSI
Pa_to_psia(P_Pa)   = P_Pa / _PA_PER_PSI

# ── Mass and mass flow ────────────────────────────────────────────────────────
const _KG_PER_LBM = 0.45359237

lbm_to_kg(m_lbm)   = m_lbm * _KG_PER_LBM
kg_to_lbm(m_kg)    = m_kg / _KG_PER_LBM
lbps_to_kgps(W_lb) = W_lb * _KG_PER_LBM
kgps_to_lbps(W_kg) = W_kg / _KG_PER_LBM

# ── Specific energy and specific heat (FPT table units) ───────────────────────
const _JKG_PER_BTULBM   = 2326.0    # 1 BTU/lbm  = 2.326 kJ/kg
const _JKGK_PER_BTULBMR = 4186.8    # 1 BTU/(lbm·°R) = 4.1868 kJ/(kg·K)

btulbm_to_Jkg(h)   = h * _JKG_PER_BTULBM
Jkg_to_btulbm(h)   = h / _JKG_PER_BTULBM
btulbmR_to_JkgK(c) = c * _JKGK_PER_BTULBMR
JkgK_to_btulbmR(c) = c / _JKGK_PER_BTULBMR

# ── Density ───────────────────────────────────────────────────────────────────
const _KGM3_PER_LBMFT3 = _KG_PER_LBM / 0.3048^3   # ≈ 16.01846

lbmft3_to_kgm3(ρ) = ρ * _KGM3_PER_LBMFT3
kgm3_to_lbmft3(ρ) = ρ / _KGM3_PER_LBMFT3

# ── Rotational speed ──────────────────────────────────────────────────────────
rpm_to_radps(N) = N * (π / 30)
radps_to_rpm(ω) = ω * (30 / π)

# ── Power ─────────────────────────────────────────────────────────────────────
const _W_PER_HP = 745.69987158227022   # mechanical horsepower (550 ft·lbf/s)

hp_to_W(P) = P * _W_PER_HP
W_to_hp(P) = P / _W_PER_HP
