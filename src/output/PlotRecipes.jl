"""
Plots.jl recipes (via RecipesBase — no Plots dependency).

  tsdiagram(sol)        T-s diagram of the solved cycle, stations labeled.
  mapplot(map)          Performance map: PR vs Wc speed lines.
  mapplot(el)           Same, with the turbomachine's operating point marked.

Requires `using Plots` (or any RecipesBase-compatible backend) on the caller's
side.
"""

using RecipesBase

# ── T-s diagram ───────────────────────────────────────────────────────────────

@userplot TsDiagram

@recipe function f(h::TsDiagram)
    length(h.args) == 1 && h.args[1] isa SolveResult ||
        error("tsdiagram expects a single SolveResult argument")
    r   = h.args[1]
    net = r.net

    sts, closed_explicit = _station_paths(net; branches = false)
    labels = first.(sts)
    s_vals = [_scalar(entropy(st)) / 1e3 for (_, st) in sts]
    T_vals = [_scalar(st.Tt)             for (_, st) in sts]

    # A closed Brayton model either walks back to the seed explicitly or ends
    # at a dangling outlet whose return to the seed state is the (unmodeled)
    # heat-rejection leg; back-edges only exist in closed-loop models.
    close_loop = closed_explicit || any(e -> e.back_edge, net.edges)

    xguide --> "s  [kJ/(kg·K)]"
    yguide --> "T  [K]"
    title  --> "T-s diagram"
    legend --> false

    @series begin
        seriestype := :path
        marker     := :circle
        markersize --> 4
        s_vals, T_vals
    end

    if close_loop && length(sts) > 2
        @series begin
            seriestype := :path
            linestyle  := :dash
            [s_vals[end], s_vals[1]], [T_vals[end], T_vals[1]]
        end
    end

    @series begin
        seriestype         := :scatter
        markersize         := 0
        markeralpha        := 0
        series_annotations := [(lab, 8, :left, :bottom) for lab in labels]
        s_vals, T_vals
    end
end

# ── Performance map ───────────────────────────────────────────────────────────

@userplot MapPlot

@recipe function f(h::MapPlot)
    el = nothing
    m  = if h.args[1] isa PerformanceMap
        h.args[1]
    elseif h.args[1] isa Union{Compressor,Turbine} && !isnothing(h.args[1].map)
        el = h.args[1]
        el.map
    else
        error("mapplot expects a PerformanceMap, or a Compressor/Turbine with a map attached")
    end

    xguide --> "Wc  [kg/s, corrected]"
    yguide --> "PR"
    title  --> (isnothing(el) ? "Performance map" : "$(el.name) map")
    legend --> :topleft

    for (i, Nc) in enumerate(m.Nc_axis)
        @series begin
            seriestype := :path
            label      := "Nc = $(round(Nc, sigdigits=4))"
            m.Wc_axis, m.PR_grid[i, :]
        end
    end

    # Operating point from the last solve
    if !isnothing(el) && !isnothing(el.inlet)
        s  = el.inlet[]
        Wc = el.mode == :off_design && el.Wc_map > 0 ? el.Wc_map :
             corrected_flow(_scalar(s.W), _scalar(s.Tt), _scalar(s.Pt))
        @series begin
            seriestype  := :scatter
            marker      := :star5
            markersize  := 10
            markercolor := :red
            label       := "operating point"
            [Wc], [_scalar(pressure_ratio(el))]
        end
    end
end
