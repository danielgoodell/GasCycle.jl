"""
Pure-Julia BRU TIT sweep plot.

`Plots.jl` is intentionally optional: GasCycle depends only on RecipesBase
for plotting recipes.  Install Plots in your user or working environment,
then run:

    julia --project=. examples/plot_bru_tit_sweep.jl
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

try
    @eval using Plots
catch err
    err isa ArgumentError || rethrow()
    error("""
    Plots.jl is optional and is not a GasCycle dependency.

    Install it in an environment on Julia's LOAD_PATH, for example:

        julia -e 'using Pkg; Pkg.add("Plots")'

    Then rerun:

        julia --project=. examples/plot_bru_tit_sweep.jl
    """)
end

include(joinpath(@__DIR__, "bru_tit_sweep.jl"))

rows = run_bru_tit_sweep(print_rows = false)
isempty(rows) && error("No converged TIT sweep points to plot")

TIT_R    = [r.TIT_R for r in rows]
W_shaft  = [r.W_shaft_kW for r in rows]
W_elec   = [r.W_elec_kW for r in rows]
eta_pct  = [r.eta_cycle_pct for r in rows]

function crossing_x(x, y, y_target)
    for i in 1:length(x)-1
        y1, y2 = y[i], y[i+1]
        (y1 - y_target) * (y2 - y_target) <= 0 || continue
        y1 == y2 && return x[i]
        return x[i] + (y_target - y1) * (x[i+1] - x[i]) / (y2 - y1)
    end
    nothing
end

colors = (shaft = "#1a4e8a",
          elec  = "#c0392b",
          goal  = "#2c7a2c",
          eta   = "#7d3c98")

default(; framestyle = :box,
          grid = true,
          gridstyle = :dash,
          minorgrid = true,
          minorgridstyle = :dot,
          linewidth = 1.8,
          markersize = 4,
          legendfontsize = 8,
          guidefontsize = 11,
          tickfontsize = 9,
          titlefontsize = 10)

p = plot(TIT_R, W_shaft;
         label = "Net shaft power",
         color = colors.shaft,
         marker = :circle,
         xlabel = "Turbine inlet temperature  (degR)",
         ylabel = "Net power  (kW)",
         xlim = (1420, 2100),
         ylim = (0, 17),
         xticks = 1400:100:2100,
         yticks = 0:2:18,
         legend = :topleft,
         size = (900, 650))

plot!(p, TIT_R, W_elec;
      label = "Est. net electrical power (eta_alt=0.92, parasitic=1.57 kW)",
      color = colors.elec,
      marker = :square)
hline!(p, [10.5]; label = "Design goal 10.5 kW", color = colors.goal,
       linestyle = :dash, linewidth = 1.2)
vline!(p, [2060]; label = false, color = :gray, linestyle = :dot, linewidth = 1.0)
annotate!(p, 2060, 2.0, text("Design TIT\n2060 degR", 8, :gray, :center))

p_eta = twinx(p)
plot!(p_eta, TIT_R, eta_pct;
      label = "Cycle thermal efficiency",
      color = colors.eta,
      linestyle = :dash,
      linewidth = 1.4,
      ylabel = "Cycle thermal efficiency  (%)",
      ylim = (0, 55),
      yticks = 0:10:50,
      legend = :bottomright)

cross_R = crossing_x(TIT_R, W_elec, 10.5)
if !isnothing(cross_R)
    annotate!(p, cross_R - 90, 12.6,
              text("$(round(Int, cross_R)) degR\n($(round(Int, R_to_K(cross_R))) K)",
                   8, colors.elec, :center))
    plot!(p, [cross_R - 50, cross_R], [12.0, 10.5];
          label = false, color = colors.elec, arrow = true, linewidth = 1.0)
end

title!(p,
       "BRU 10 kW Brayton Cycle - Net Power vs Turbine Inlet Temperature\n" *
       "GasCycle.jl (HeXe84.fpt), PR_c=1.9, mdot=0.599 kg/s, epsilon_rec=0.95")

outpath = joinpath(@__DIR__, "plot_bru_tit_sweep_julia.png")
savefig(p, outpath)
println("Saved -> $outpath")
