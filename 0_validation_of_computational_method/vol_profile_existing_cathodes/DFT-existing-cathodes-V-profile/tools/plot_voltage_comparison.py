"""Plot original (VP/voltage.csv) vs recomputed (variant B) voltage profiles, one panel per system.

Usage: python tools/plot_voltage_comparison.py <clean_root>
Writes <clean_root>/recomputed/voltage_comparison.{png,pdf}
"""
import csv, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(sys.argv[1]).resolve()
rows = list(csv.DictReader((root / "recomputed" / "voltage_comparison.csv").open()))

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
SERIES = {  # fixed categorical order: slot 1 blue, slot 2 orange
    "original": dict(color="#2a78d6", ls="-", label="Original (VP/voltage.csv)"),
    "B_qe_energies+qe_mu": dict(color="#eb6834", ls=(0, (4, 2)), label="Recomputed (variant B)"),
}
ORDER = ["LiCoO2", "NaCoO2", "LiFePO4", "Li2MnO3", "MgV2O4", "KVPO4F"]
ION = {"LiCoO2": "Li", "NaCoO2": "Na", "LiFePO4": "Li", "Li2MnO3": "Li", "MgV2O4": "Mg", "KVPO4F": "K"}

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
fig, axes = plt.subplots(2, 3, figsize=(11, 6.2), facecolor=SURFACE)
for ax, s in zip(axes.flat, ORDER):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    avgs = {}
    for var, st in SERIES.items():
        rr = [r for r in rows if r["system"] == s and r["variant"] == var]
        xs, ys = [], []
        for r in rr:
            xs += [float(r["x1"]), float(r["x2"])]
            ys += [float(r["V"])] * 2
        ax.plot(xs, ys, color=st["color"], ls=st["ls"], lw=2, solid_capstyle="round",
                label=st["label"], zorder=3 if var == "original" else 4)
        avgs[var] = float(rr[0]["avg_V"])
    d = round(avgs["B_qe_energies+qe_mu"] - avgs["original"], 3) + 0.0
    ax.set_title(s, loc="left", fontsize=10, color=INK, fontweight="bold")
    ax.text(0.98, 0.95, f"avg {avgs['original']:.3f} → {avgs['B_qe_energies+qe_mu']:.3f} V\n(Δ {d:+.3f} V)",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color=INK2)
    ax.set_xlim(0, 1)
    ax.set_xlabel(f"x = n({ION[s]}) / n_max")
    ax.set_ylabel("Voltage (V)")
handles, labels = axes.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.0))
fig.suptitle("0 K voltage profiles: original vs recomputed (one Ry→eV constant, QE-based μ)",
             y=0.96, fontsize=11, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.93))
for ext in ("png", "pdf"):
    fig.savefig(root / "recomputed" / f"voltage_comparison.{ext}", dpi=200, facecolor=SURFACE)
print("written", root / "recomputed" / "voltage_comparison.png")
