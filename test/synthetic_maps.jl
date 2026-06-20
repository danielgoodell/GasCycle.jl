# Synthetic native-coordinate turbomachine maps for the off-design testsets,
# built as FunctionMap "scripts" that pass exactly through a given design point.
#
#   Compressor (R-line native):  (Nc, Rline) -> (Wc, PR, eff)
#     Wc rises with R and speed; PR falls with R (physical compressor slope).
#   Turbine (PR native):         (Np, PR)    -> (Wp, eff)
#     Wp rises monotonically with PR so flow continuity has a unique root.
#
# Each closure equals the design outputs at (design speed, R_des / PR_des), so an
# off-design solve at design boundary conditions reproduces the design point.

function synthetic_comp_map(; Nc_des, Wc_des, PR_des, eta_des, R_des = 2.0)
    f = function (Nc, R, rc)
        sN  = Nc / Nc_des
        Wc  = Wc_des * sN * (1 + 0.30 * (R - R_des))
        PR  = 1 + (PR_des - 1) * sN^2 * (1 - 0.10 * (R - R_des))
        eff = eta_des - 0.05 * (R - R_des)^2 - 0.10 * (sN - 1)^2
        (; Wc, PR, eff)
    end
    FunctionMap(f; line_des = R_des)
end

function synthetic_turb_map(; Np_des, Wp_des, PR_des, eta_des)
    f = function (Np, PR, rc)
        sN  = Np / Np_des
        # Linear, monotonically increasing in PR (well-conditioned flow root);
        # efficiency tracks speed/loading, decoupled from the PR unknown.
        Wp  = Wp_des * (1 + 0.8 * (PR - PR_des) / (PR_des - 1)) * (1 + 0.05 * (sN - 1))
        eff = eta_des - 0.10 * (sN - 1)^2
        (; Wp, eff)
    end
    FunctionMap(f)
end
