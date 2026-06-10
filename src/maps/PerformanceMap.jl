"""
2-D turbomachinery performance map.

Stores corrected mass flow (Wc) vs. corrected speed (Nc) grid with
corresponding pressure ratio (PR) and isentropic efficiency (η) values.

Follows NPSS corrected variable convention:
  Nc = N / sqrt(Tt / T_std)          [rpm / sqrt(K/K)]
  Wc = W * sqrt(Tt / T_std) / Pt * P_std   [kg/s]

where T_std = 288.15 K, P_std = 101325 Pa.
"""

using Interpolations

const T_STD = 288.15   # K
const P_STD = 101325.0 # Pa

struct PerformanceMap
    Nc_axis::Vector{Float64}    # corrected speed grid  (nN points)
    Wc_axis::Vector{Float64}    # corrected flow grid   (nW points)
    PR_grid::Matrix{Float64}    # pressure ratio        [nN × nW]
    eta_grid::Matrix{Float64}   # isentropic efficiency [nN × nW]
    bounds::Symbol              # :error, :warn, or :clamp

    itp_PR::Any
    itp_eta::Any

    function PerformanceMap(Nc_axis, Wc_axis, PR_grid, eta_grid; bounds::Symbol = :error)
        size(PR_grid)  == (length(Nc_axis), length(Wc_axis)) ||
            error("PR_grid size must be (nNc, nWc)")
        size(eta_grid) == (length(Nc_axis), length(Wc_axis)) ||
            error("eta_grid size must be (nNc, nWc)")
        bounds in (:error, :warn, :clamp) ||
            error("PerformanceMap bounds must be :error, :warn, or :clamp, got :$bounds")

        itp_PR  = interpolate((Nc_axis, Wc_axis), PR_grid,  Gridded(Linear()))
        itp_eta = interpolate((Nc_axis, Wc_axis), eta_grid, Gridded(Linear()))
        new(collect(Nc_axis), collect(Wc_axis),
            PR_grid, eta_grid, bounds, itp_PR, itp_eta)
    end
end

function _clamp_map_point(m::PerformanceMap, Nc, Wc)
    Nc_lo, Nc_hi = m.Nc_axis[1], m.Nc_axis[end]
    Wc_lo, Wc_hi = m.Wc_axis[1], m.Wc_axis[end]
    in_bounds = Nc_lo <= Nc <= Nc_hi && Wc_lo <= Wc <= Wc_hi
    in_bounds && return Nc, Wc

    msg = "PerformanceMap query outside map bounds: requested " *
          "Nc=$Nc, Wc=$Wc; valid ranges are " *
          "Nc=[$Nc_lo, $Nc_hi], Wc=[$Wc_lo, $Wc_hi]. " *
          "Construct the map with bounds=:warn or bounds=:clamp to clamp intentionally."

    if m.bounds == :error
        throw(DomainError((Nc=Nc, Wc=Wc), msg))
    elseif m.bounds == :warn
        @warn msg
    end

    clamp(Nc, Nc_lo, Nc_hi), clamp(Wc, Wc_lo, Wc_hi)
end

"""
    query(map, Nc, Wc) -> (PR, η)

Interpolate the map at corrected speed Nc and corrected flow Wc.
Bounds behavior is controlled by `PerformanceMap(...; bounds=...)`:
`:error` throws on off-map queries, `:warn` warns then clamps, and `:clamp`
silently preserves the old clamping behavior.
"""
function query(m::PerformanceMap, Nc::Float64, Wc::Float64)
    Nc_c, Wc_c = _clamp_map_point(m, Nc, Wc)
    (m.itp_PR(Nc_c, Wc_c), m.itp_eta(Nc_c, Wc_c))
end

"""
    corrected_speed(N, Tt) -> Nc

Convert physical shaft speed N [rpm] and inlet total temperature Tt [K]
to corrected speed.
"""
corrected_speed(N::Float64, Tt::Float64) = N / sqrt(Tt / T_STD)

"""
    corrected_flow(W, Tt, Pt) -> Wc

Convert mass flow W [kg/s], total temperature Tt [K], and total pressure
Pt [Pa] to corrected mass flow.
"""
corrected_flow(W::Float64, Tt::Float64, Pt::Float64) =
    W * sqrt(Tt / T_STD) / Pt * P_STD
