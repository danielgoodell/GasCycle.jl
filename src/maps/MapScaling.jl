"""
Scale a generic performance map to a specific design point.

Given a baseline map and the desired design-point operating conditions
(Nc_des, Wc_des, PR_des, eta_des), returns a new PerformanceMap whose
design point matches exactly while the overall map shape is preserved.

Scale factors:
  s_Nc  = Nc_des  / Nc_ref     (corrected speed at design point on base map)
  s_Wc  = Wc_des  / Wc_ref
  s_PR  = (PR_des - 1) / (PR_ref - 1)   (preserve PR=1 at zero flow)
  s_eta = eta_des / eta_ref
"""
function scale_map(base::PerformanceMap;
                   Nc_des::Float64,
                   Wc_des::Float64,
                   PR_des::Float64,
                   eta_des::Float64,
                   Nc_ref::Float64 = base.Nc_axis[end ÷ 2 + 1],
                   Wc_ref::Float64 = base.Wc_axis[end ÷ 2 + 1])

    PR_ref, eta_ref = query(base, Nc_ref, Wc_ref)

    s_Nc  = Nc_des / Nc_ref
    s_Wc  = Wc_des / Wc_ref
    s_PR  = (PR_des - 1.0) / (PR_ref - 1.0)
    s_eta = eta_des / eta_ref

    Nc_new  = base.Nc_axis .* s_Nc
    Wc_new  = base.Wc_axis .* s_Wc
    PR_new  = 1.0 .+ (base.PR_grid  .- 1.0) .* s_PR
    eta_new = base.eta_grid .* s_eta

    PerformanceMap(Nc_new, Wc_new, PR_new, eta_new)
end
