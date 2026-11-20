"""
Plots.jl recipes (via RecipesBase — no Plots dependency).

  tsdiagram(sol)        T-s diagram of the solved cycle, stations labeled.
  mapplot(map)          Performance map: PR vs Wc speed lines.
  mapplot(el)           Same, with the turbomachine's operating point marked.
  losswaterfall(L)      Exergy waterfall of a loss ledger: input → losses → net power.

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
    m  = if h.args[1] isa TurbomachineMap
        h.args[1]
    elseif h.args[1] isa Union{Compressor,Turbine} && !isnothing(h.args[1].map)
        el = h.args[1]
        el.map
    else
        error("mapplot expects a TurbomachineMap, or a Compressor/Turbine with a map attached")
    end
    m isa FunctionMap &&
        error("mapplot cannot draw a FunctionMap (no tabulated speed lines)")

    xguide --> "corrected flow  [kg/s]"
    yguide --> "PR"
    title  --> (isnothing(el) ? "Performance map" : "$(el.name) map")
    legend --> :topleft

    for (spd, x, PR) in speed_lines(m)
        @series begin
            seriestype := :path
            label      := "Nc = $(round(spd, sigdigits=4))"
            x, PR
        end
    end

    # Operating point from the last solve
    if !isnothing(el) && !isnothing(el.inlet)
        s  = el.inlet[]
        Wc = corrected_flow(_scalar(s.W), _scalar(s.Tt), _scalar(s.Pt))
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

# ── Loss waterfall ────────────────────────────────────────────────────────────

@userplot LossWaterfall

# Rectangles [x0,x1]×[y0,y1] as one :shape series (NaN-separated).
function _rects(xs, lo, hi; w = 0.8)
    X, Y = Float64[], Float64[]
    for (x, a, b) in zip(xs, lo, hi)
        append!(X, [x - w/2, x + w/2, x + w/2, x - w/2, x - w/2, NaN])
        append!(Y, [a, a, b, b, a, NaN])
    end
    X, Y
end

@recipe function f(h::LossWaterfall)
    length(h.args) == 1 && h.args[1] isa Union{LossLedger,SolveResult} ||
        error("losswaterfall expects a LossLedger (or a SolveResult)")
    L = h.args[1] isa SolveResult ? loss_ledger(h.args[1]) : h.args[1]

    Ex   = _primal(sum(last, L.exergy_in; init = 0.0)) / 1e3
    rows = ledger_rows(L)
    labels = ["exergy in"; ["$(r[1]) ($(r[2]))" for r in rows];
              L.electrical ? "net electric" : "net power"]
    n = length(labels)

    # Running level: each loss bar hangs from the level left by the previous one.
    level = Ex
    lo, hi = Float64[], Float64[]
    for r in rows
        v = r[3] / 1e3
        push!(lo, level - v); push!(hi, level)
        level -= v
    end

    xticks  --> (1:n, labels)
    xrotation --> 45
    yguide  --> "exergy  [kW]"
    title   --> "Loss waterfall (T0 = $(round(_primal(L.T0), digits = 1)) K)"
    legend  --> false

    @series begin
        seriestype := :shape
        fillcolor  --> :steelblue
        _rects([1], [0.0], [Ex])
    end
    @series begin
        seriestype := :shape
        fillcolor  --> :indianred
        _rects(2:n-1, lo, hi)
    end
    @series begin
        seriestype := :shape
        fillcolor  --> :seagreen
        _rects([n], [0.0], [_primal(L.P_out) / 1e3])
    end
end
