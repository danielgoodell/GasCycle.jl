# ── Cycle-level performance quantities ────────────────────────────────────────

"""
    net_power(result) -> Real

Net shaft power [W] = total turbine work - total compressor work.
"""
function net_power(r::SolveResult)
    W_turb = 0.0
    W_comp = 0.0
    for el in r.net.elements
        if el isa Turbine && !isnothing(el.inlet)
            W_turb += specific_work(el) * el.inlet[].W
        elseif el isa Compressor && !isnothing(el.inlet)
            W_comp += specific_work(el) * el.inlet[].W
        end
    end
    W_turb - W_comp
end

"""
    cycle_efficiency(result) -> Real

Thermal efficiency = net power / total heat input from all HeatSource elements.
"""
function cycle_efficiency(r::SolveResult)
    Q_in = sum(r.net.elements; init=0.0) do el
        if el isa HeatSource
            if !isnothing(el.inlet) && !isnothing(el.outlet)
                (enthalpy(el.outlet[]) - enthalpy(el.inlet[])) * el.inlet[].W
            else
                el.Q
            end
        else
            0.0
        end
    end
    Q_in ≈ 0.0 && return 0.0
    net_power(r) / Q_in
end
