from __future__ import unicode_literals

from functools import partial
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend import Legend
from matplotlib.patches import Patch

import pandas as pd
from scipy.io import mmread
from scipy.optimize import curve_fit
from uncertainties import ufloat

from plot_scaling import UncertaintyScalePlot
from plot_sparse import LogSparsePlot

def amdahl(p, sigma, baseline=1):
    return baseline * p / (1 + sigma * (p - 1))

def USL(p, sigma, kappa, baseline=1):
    return baseline * p / (1 + sigma * (p - 1) + kappa * p * (p - 1))

if __name__ == "__main__":
    all = pd.read_json("all_73655cb.json")

    all.loc[all["suite"] == "no-pysparse", "suite"] = "trilinos"

    by = ["suite", "fipy_rev", "nx"]
    data_set="elapsed / s"
    groups = all.groupby(by + ["tasks"])
    groups = groups.agg(data_mean=(data_set, "mean")).reset_index()
    groups = groups.groupby(by)
    for nx in all["nx"].unique():
        baseline = groups.get_group(("pysparse", "371d28468", nx)).iloc[0]["data_mean"]
        all.loc[all["nx"] == nx, "speedup"] = baseline / all[data_set]
        
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1)
    scaling = UncertaintyScalePlot(ax=ax,
                                   line_style={
                                       'a5f233aa7': ("", "full"),
                                       '371d28468': ("", "none"),
                                       "ef3db887e": ("", "full")
                                   })
    scaling.plot(all[all["nx"] == 2048], by=["suite", "fipy_rev"],
                 xdata="tasks", ydata="speedup")

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("tasks")
    ax.set_ylabel("speedup / ($t_{PySparse} / t_N$)")
    ax.set_xlim(xmin=None, xmax=None)
    ax.set_ylim(ymin=1e-1, ymax=3e1)

    scaling.suite_legend(suites=all[all["nx"] == 2048]["suite"], linestyle="")

    fit_tasks = np.logspace(0, 2, 100)

    groups = all[(all["nx"] == 2048)].groupby(by=["suite", "fipy_rev"])
    
    for suite, rev, n, color in [["petsc", "a5f233aa7", 1, "blue"],
                                 ["trilinos", "a5f233aa7", 0, "red"],
                                 ["petsc", "371d28468", 0, "blue"],
                                 ["trilinos", "371d28468", 0, "red"]]:
        print(f"{suite} {rev}")
        
        d = groups.get_group((suite, rev))

        d = d[d["tasks"] >= 2**n]
        
        baseline = d.groupby("tasks")["speedup"].mean().iloc[0] / 2**n

#         ax.plot(d["tasks"], d["speedup"], marker="x", linestyle="", color=color)

        amdahl_bl = partial(amdahl, baseline=baseline)
        popt, pcov = curve_fit(amdahl_bl, d["tasks"], d["speedup"])
        ax.plot(fit_tasks, amdahl_bl(fit_tasks, *popt), color=color, linestyle="--")
        print(f"Amdahl")
        sigma = ufloat(popt[0], np.sqrt(pcov[0][0]))
        print(f"{sigma = :0.1u%S}")

        USL_bl = partial(USL, baseline=baseline)
        popt, pcov = curve_fit(USL_bl, d["tasks"], d["speedup"])
        ax.plot(fit_tasks, USL_bl(fit_tasks, *popt), color=color, linestyle="-")
        print(f"Gunther")
        sigma = ufloat(popt[0], np.sqrt(pcov[0][0]))
        kappa = ufloat(popt[1], np.sqrt(pcov[1][1]))
        print(f"{sigma = :0.1u%S}")
        print(f"{kappa = :0.1uS}")
        print()

    legend_elements = [Line2D([0], [0], color="black", marker="o", fillstyle="none", linestyle=""),
                       Line2D([0], [0], color="black", marker="o", linestyle=""),
                       Patch(facecolor="black", edgecolor=None, alpha=0.1),
                       Line2D([0], [0], color="black", marker="", linestyle="--"),
                       Line2D([0], [0], color="black", marker="", linestyle="-"),
                       Line2D([0], [0], color="black", linestyle="-", linewidth=1, alpha=0.1)]
    leg = Legend(ax, handles=legend_elements,
                 labels=["FiPy 3.4.4 (371d28468)", "FiPy 4.0 (a5f233aa7)",
                         "$\pm 1$ standard deviation",
                         "Amdahl's Law", "Gunther's Law", r"$\sim N$"],
                 loc=(0.7, 0.35)) #, frameon=False)
    ax.add_artist(leg);
                 
    for n in range(-3, 2):
        ax.axline((1, 10**n), (10, 10**(n+1)), color="black", linewidth=1, alpha=0.1)

    outer_grid = fig.add_gridspec(3, 3)
    inner_fig = fig.add_subfigure(outer_grid[2, 1])
    inner_grid = inner_fig.add_gridspec(10, 10)
    
    sparse = LogSparsePlot(coo=mmread("parallel_scaling.mtx"),
                           rhs=np.load("parallel_scaling.rhs.npz")["rhs"],
                           fig=inner_fig.add_subfigure(inner_grid[:6, :9]))
    sparse.plot()
    
    data = np.load("parallel_initial.npz")
    soln_ax = fig.add_subplot(outer_grid[2, 2],
                              position=[0.645, 0.11, 0.23, 0.23])
    soln_ax.imshow(data["eta"])
    soln_ax.set_title("initial condition")
    soln_ax.set(xticks=[], yticks=[])
    
    plt.show()
