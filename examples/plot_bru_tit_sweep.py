"""
Plot BRU TIT sweep results, reproducing the style of Figure 6 from
Wong et al., NASA TN D-5815 (May 1970).

Run after: julia examples/bru_tit_sweep.jl > /tmp/bru_sweep.csv
Or pipe directly:  ~/.juliaup/bin/julia examples/bru_tit_sweep.jl | python3 examples/plot_bru_tit_sweep.py
"""

import sys
import io
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Run Julia sweep if not piped ───────────────────────────────────────────────
def load_data():
    if not sys.stdin.isatty():
        raw = sys.stdin.read()
    else:
        print("Running Julia sweep...", flush=True)
        import os
        julia = os.path.expanduser("~/.juliaup/bin/julia")
        script = os.path.join(os.path.dirname(__file__), "bru_tit_sweep.jl")
        result = subprocess.run([julia, script], capture_output=True, text=True)
        raw = result.stdout

    # Find the CSV header line
    lines = [l for l in raw.splitlines() if l.startswith("TIT_R") or
             (len(l) > 0 and l[0].isdigit())]
    text = "\n".join(lines)
    data = np.genfromtxt(io.StringIO(text), delimiter=",", names=True)
    return data

data = load_data()

TIT_R        = data["TIT_R"]
W_shaft      = data["W_shaft_kW"]
W_elec       = data["W_elec_kW"]
Q_heater     = data["Q_heater_kW"]
eta_pct      = data["eta_cycle_pct"]

# ── Figure setup ──────────────────────────────────────────────────────────────
fig, ax1 = plt.subplots(figsize=(7.5, 5.5))

# Gray style reminiscent of NASA report line drawings
colors = {"shaft": "#1a4e8a", "elec": "#c0392b", "goal": "#2c7a2c"}

ax1.plot(TIT_R, W_shaft, "-o", color=colors["shaft"], linewidth=1.8,
         markersize=4, label="Net shaft power")
ax1.plot(TIT_R, W_elec,  "-s", color=colors["elec"],  linewidth=1.8,
         markersize=4, label="Est. net electrical power\n(η$_{alt}$=0.92, parasitic=1.57 kW)")

# Design goal line
ax1.axhline(10.5, color=colors["goal"], linewidth=1.2, linestyle="--",
            label="Design goal  10.5 kW")

# Design-point marker at 2060 °R
ax1.axvline(2060, color="gray", linewidth=0.8, linestyle=":")
ax1.annotate("Design TIT\n2060 °R", xy=(2060, 2.0),
             ha="center", va="bottom", fontsize=8.5, color="gray")

ax1.set_xlabel("Turbine inlet temperature  (°R)", fontsize=11)
ax1.set_ylabel("Net power  (kW)", fontsize=11)
ax1.set_xlim(1420, 2100)
ax1.set_ylim(0, 17)
ax1.yaxis.set_major_locator(ticker.MultipleLocator(2))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(1))
ax1.xaxis.set_major_locator(ticker.MultipleLocator(100))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax1.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.6)
ax1.grid(True, which="minor", linestyle=":",  linewidth=0.3, alpha=0.4)

# Secondary y-axis: cycle thermal efficiency
ax2 = ax1.twinx()
ax2.plot(TIT_R, eta_pct, "--", color="#7d3c98", linewidth=1.4, alpha=0.8,
         label="Cycle thermal efficiency")
ax2.set_ylabel("Cycle thermal efficiency  (%)", color="#7d3c98", fontsize=11)
ax2.tick_params(axis="y", labelcolor="#7d3c98")
ax2.set_ylim(0, 55)
ax2.yaxis.set_major_locator(ticker.MultipleLocator(10))

# ── Second x-axis in Kelvin ────────────────────────────────────────────────────
ax3 = ax1.twiny()
ax3.set_xlim(ax1.get_xlim())
K_ticks_R = np.arange(1460, 2100, 100)
ax3.set_xticks(K_ticks_R)
ax3.set_xticklabels([f"{int(T*(5/9))}" for T in K_ticks_R], fontsize=8)
ax3.set_xlabel("Turbine inlet temperature  (K)", fontsize=10)

# ── Legend ─────────────────────────────────────────────────────────────────────
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2,
           loc="upper left", fontsize=8.5, framealpha=0.9)

# ── Annotations ────────────────────────────────────────────────────────────────
# Mark where electrical output crosses design goal
from scipy.interpolate import interp1d
f_elec = interp1d(TIT_R, W_elec)
# find crossings numerically
cross_R = None
for i in range(len(TIT_R)-1):
    if (W_elec[i] - 10.5) * (W_elec[i+1] - 10.5) < 0:
        cross_R = float(np.interp(10.5, [W_elec[i], W_elec[i+1]], [TIT_R[i], TIT_R[i+1]]))
        break

if cross_R is not None:
    ax1.annotate(f"{cross_R:.0f} °R\n({cross_R*5/9:.0f} K)",
                 xy=(cross_R, 10.5), xytext=(cross_R - 100, 12.5),
                 arrowprops=dict(arrowstyle="->", color=colors["elec"], lw=1.0),
                 color=colors["elec"], fontsize=8, ha="center")

fig.suptitle(
    "BRU 10 kW Brayton Cycle — Net Power vs Turbine Inlet Temperature\n"
    "GasCycle.jl (HeXe84.fpt)  ·  PR$_c$=1.9, $\\dot{m}$=0.599 kg/s, "
    "$\\varepsilon_{rec}$=0.95  ·  ref: NASA TN D-5815",
    fontsize=9.5
)

plt.tight_layout(rect=[0, 0, 1, 0.95])

outpath = __file__.replace(".py", ".png")
plt.savefig(outpath, dpi=150, bbox_inches="tight")
print(f"Saved → {outpath}")
plt.show()
