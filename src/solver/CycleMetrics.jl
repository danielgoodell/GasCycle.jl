# ── Cycle-level performance quantities ────────────────────────────────────────

"""
    net_power(result) -> Real
    net_power(net)    -> Real

Net shaft power [W] = total turbine work - total compressor work.
The `net` form evaluates the current state, so it can be used as a target:
`add_target!(net, net_power, 10e3)`.
"""
net_power(r::SolveResult) = net_power(r.net)

function net_power(net::FlowNetwork)
    W_turb = 0.0
    W_comp = 0.0
    for el in net.elements
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
    cycle_efficiency(net)    -> Real

Thermal efficiency = net power / total heat input from all HeatSource elements.
"""
cycle_efficiency(r::SolveResult) = cycle_efficiency(r.net)

function cycle_efficiency(net::FlowNetwork)
    Q_in = sum(net.elements; init=0.0) do el
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
    net_power(net) / Q_in
end
